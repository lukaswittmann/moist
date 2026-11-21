!> Reverse-mode nuclear gradient of the DROP cavity
!>
!> [[get_surface_gradient_drop]] contracts a surface adjoint straight into
!> `dE/dR_A` without ever forming the forward Jacobian. These tests pin it
!> against the forward `*_rA` arrays, which are independently validated by the
!> finite-difference suite in gradient.f90
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
!> drives them one at a time so a bug in one cannot hide behind another
module test_cavity_drop_nuclear_adjoint
!$ use omp_lib, only: omp_in_parallel, omp_get_num_threads
   use moist_cavity_drop_lsf_svdw_param, only: moist_cavity_drop_lsf_svdw_param_type
   use moist_cavity_drop_parameters, only: moist_cavity_drop_parameters_type, drop_branch_parameters_type
   use moist_model_continuum_component_pcm_type, only: moist_pcm_parameters_type
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use mctc_env_accuracy, only: wp
   use mctc_env_error, only: mctc_error => error_type
   use mctc_io, only: structure_type, new
   use testdrive, only: new_unittest, unittest_type, error_type, check, to_string, test_failed, &
      & skip_test
   use test_helpers, only: get_test_cross, fd6_offsets, get_test_cube, nudge_structure, get_test_forward_pair, &
                           get_test_sibling_gate
   use moist_cavity_drop, only: cavity_type_drop, new_cavity_drop
   use moist_cavity_drop_tile, only: surface_tile_width
   use moist_cavity_drop_lsf_svdw, only: moist_cavity_drop_lsf_svdw_type
   use moist_cavity_drop_lsf_base, only: moist_cavity_drop_lsf_type
   use moist_cavity_drop_gaussian, only: iswig_workspace_type
   use moist_cavity_drop_lsf_isodensity_callback, only: moist_cavity_drop_lsf_isodensity_callback_type
   use moist_cavity_drop_lsf_isodensity_internal, only: moist_cavity_drop_lsf_isodensity_internal_type
   use moist_cavity_drop_lsf_isodensity_param, only: moist_cavity_drop_lsf_isodensity_param_type
   use, intrinsic :: iso_c_binding, only: c_ptr, c_double, c_int, c_funloc, c_null_ptr, c_associated, c_f_pointer
   use moist_cavity_adjoint, only: cavity_adjoint_type
   use moist_radii, only: default_cpcm_radii, new_radii_custom_atoms, radius_type
   use moist_data_radii_legacy, only: get_radius_func
   use moist_context, only: moist_context_type, new_context
   use moist_channels_coupling, only: coupling_type, density_matrix_request_type, &
      & nuclear_charge_request_type
   use moist_channels_response, only: response_type, potential_adjoint_response_type
   use moist_model_continuum, only: model_continuum_type, new_continuum_model
   use moist_model_continuum_component_pcm_conductor, only: model_continuum_component_cpcm, new_component_cpcm, &
      & model_continuum_component_cosmo, new_component_cosmo, model_continuum_component_pcm_conductor
   use moist_model_continuum_component_pcm_type, only: solver_type
   use moist_model_continuum_component, only: model_continuum_component_pv, new_component_pv
   use test_helpers, only: stage_model_point_charge_energy, fill_missing_with_zeros, &
      & fill_point_charge_field, &
      & fill_legacy_radii, copy_potential_adjoint
   use moist_cavity_type, only: cavity_type
   implicit none(type, external)
   private

   public :: collect_cavity_drop_nuclear_adjoint

   integer, parameter :: ndim = 3

   !> Level-set smoothing of the shared fixture
   real(wp), parameter :: k = 2.5_wp
   real(wp), parameter :: gamma = 1.0_wp
   integer, parameter :: NUM_LEB = 50

   real(wp), parameter :: PROJ_TOL = 1E-14_wp
   integer, parameter :: PROJ_LEVEL = 1

   !> Absolute forward-versus-reverse bound across all channels
   real(wp), parameter :: EQ_ABS = 5.0E-10_wp
   !> Relative bound on each forward reference entry
   real(wp), parameter :: EQ_REL = 1.0E-12_wp

   !> Model paths sum the same analytic terms with round-off differences
   real(wp), parameter :: AB_ABS = 1.0E-12_wp
   !> Relative model bound on each forward reference entry
   real(wp), parameter :: AB_REL = 1.0E-12_wp

   !> Below this the reference gradient is too small to carry a relative test
   real(wp), parameter :: VACUITY_THR = 1.0E-6_wp

   !> Branch window for the branching test
   real(wp), parameter :: BRANCH_DPDIST_MAX = 3.7_wp

   !> Branching-width finite-difference step
   real(wp), parameter :: XI_FD_STEP = 1.0E-4_wp
   !> Absolute branching-width derivative bound
   real(wp), parameter :: XI_FD_ABS = 1.0E-10_wp
   !> Relative branching-width derivative bound
   real(wp), parameter :: XI_FD_REL = 5.0E-9_wp

   !> Symmetry-breaking displacement for branching fixture
   real(wp), parameter :: FIXTURE_NUDGE = 1.0E-4_wp

   !> Channel identifiers, in the order [[apply_channel]] understands
   integer, parameter :: CH_XI = 1, CH_F = 2, CH_A = 3, CH_W = 4
   integer, parameter :: CH_XYZ = 5, CH_N = 6, CH_K1 = 7, CH_K2 = 8
   integer, parameter :: NCHAN = 8

   !> SvdW fade onset of [[test_forward_screened_neighbour]] (e-folds); band end
   !> 4.5*0.4/5.5 = 0.327 Bohr past a sphere
   real(wp), parameter :: NB_FADE = 0.4_wp
   !> Reach gate of the same fixture (Bohr): closed before the fade band end, so
   !> a point stays inside its owner's reach
   real(wp), parameter :: NB_REACH_ON = 0.2_wp, NB_REACH_OFF = 0.32_wp
   !> Finite-difference bound of the same fixture
   real(wp), parameter :: NB_FD_TOL = 2.0e-4_wp
   !> Smallest area gradient the screened neighbours must carry, five times the bound
   real(wp), parameter :: NB_AREA_MIN = 5.0_wp*NB_FD_TOL

   !> Centers of the callback density of [[test_forward_isodensity_owner]] (Bohr)
   real(wp), parameter :: density_centers(3, 2) = reshape([0.0_wp, 0.0_wp, -1.4_wp, 0.3_wp, 0.0_wp, 1.5_wp], [3, 2])

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
                  new_unittest("branching_f_fd", test_branching_f_fd), &
                  new_unittest("reach_gate_channels", test_reach_gate_channels), &
                  new_unittest("forward_pair_channels", test_forward_pair_channels), &
                  new_unittest("sibling_gate_channels", test_sibling_gate_channels), &
                  new_unittest("forward_screened_neighbour", test_forward_screened_neighbour), &
                  new_unittest("forward_screened_owner", test_forward_screened_owner), &
                  new_unittest("forward_isodensity_owner", test_forward_isodensity_owner), &
                  new_unittest("forward_frozen_field", test_forward_frozen_field), &
                  new_unittest("model_forward_reverse", test_model_forward_reverse), &
                  new_unittest("cosmo_model_forward_reverse", test_cosmo_model_forward_reverse), &
                  new_unittest("shape_guard", test_shape_guard), &
                  new_unittest("thread_count_invariance", test_thread_count_invariance), &
                  new_unittest("isodensity_thread_count_invariance", &
                               test_isodensity_thread_count_invariance) &
                  ]

   end subroutine collect_cavity_drop_nuclear_adjoint

   !> Every surface-adjoint channel driven simultaneously
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
   !> branch-weight reverse pass and its explicit owner term completely untested
   !> A near-symmetric dimer under multistart projection does branch, and
   !> [[run_equivalence]] asserts that it actually did
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
   !> checked in isolation with its own non-vacuity guard
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
   !> The width carries no branch weight, but its rows still come from a branched
   !> projection, so it is pinned on the branching fixture as well
   subroutine test_branching_xi_fd(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      call run_branching_fd(CH_XI, error)

   end subroutine test_branching_xi_fd

   !> Finite-difference check of the switching channel under branching
   !>
   !> `f = f_anchor * wbranch` is where the branch weight enters, so this is the
   !> FD check of the branch-weight reverse pass and forward post-pass
   subroutine test_branching_f_fd(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      call run_branching_fd(CH_F, error)

   end subroutine test_branching_f_fd

   !> Finite-difference check of one per-point channel under branching
   !>
   !> `test_branching_channels` only proves the reverse path matches the
   !> forward one. This pins both against an independent 6-point stencil on
   !> the scalar functional
   !>
   !>   L(R) = sum_i w_i * q_i(R),   q = xi or f
   !>
   !> Grid points are filtered and reordered on every rebuild, so the weights
   !> are keyed on the persistent `cavity%numbering` and restricted to points
   !> that survive at every stencil geometry -- otherwise a point appearing or
   !> vanishing would put a step discontinuity into L
   !>
   !> @param[in]  channel `CH_XI` or `CH_F`
   !> @param[out] error   Error handle
   subroutine run_branching_fd(channel, error)
      !> `CH_XI` or `CH_F`
      integer, intent(in) :: channel
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      !> O(h^6) central-difference stencil: with the 79 branched anchors of the
      !> level-3 cross the O(h^4) truncation error exceeds XI_FD_REL (2.3e-8)
      integer, parameter :: NSTEP = 6
      !> Coordinate along which the branch weights respond
      integer, parameter :: FD_ATOM = 5, FD_AXIS = 3
      !> Smallest stencil derivative of a branch weight that counts as a response
      real(wp), parameter :: BRANCH_RESPONSE_THR = 1.0E-3_wp

      type(cavity_type_drop), allocatable :: cavity
      type(moist_context_type), target :: ctx
      type(moist_context_type), target :: ctx_fd
      type(cavity_adjoint_type) :: acc
      type(mctc_error), allocatable :: cav_error
      type(structure_type) :: mol, mol_fd

      real(wp) :: fd_delta(NSTEP)
      real(wp), allocatable :: grad_ana(:, :), wmap(:), q_store(:, :), wb_store(:, :), w_q(:)
      character(len=:), allocatable :: label
      logical, allocatable :: have(:, :), usable(:), in_ref(:)
      !> Branch count per persistent numbering at each stencil geometry
      integer, allocatable :: bc_store(:, :)
      real(wp) :: num_deriv, ana_deriv, diff, dwb, dwb_max
      integer :: max_num, istep, igrid, inum, ngrid, ncommon, nphantom

      select case (channel)
      case (CH_XI)
         label = "xi"
      case (CH_F)
         label = "f"
      case default
         call test_failed(error, "no FD check for channel "//to_string(channel))
         return
      end select
      fd_delta = fd6_offsets*XI_FD_STEP

      call fixture_geometry(.true., mol)
      call build_cavity(cavity, ctx, mol, .true., error)
      if (allocated(error)) return

      ngrid = cavity%ngrid
      if (.not. any(cavity%offset(2:cavity%anchor%ngrid + 1) - cavity%offset(1:cavity%anchor%ngrid) > 1)) then
         call test_failed(error, "FD fixture did not branch")
         return
      end if
      max_num = maxval(cavity%numbering(1:ngrid))

      allocate (wmap(max_num))
      do inum = 1, max_num
         wmap(inum) = ramp(inum, 0)
      end do
      allocate (q_store(max_num, NSTEP), source=0.0_wp)
      allocate (wb_store(max_num, NSTEP), source=0.0_wp)
      allocate (have(max_num, NSTEP), source=.false.)
      allocate (bc_store(max_num, NSTEP), source=0)

      ! Collect the channel at each stencil geometry, keyed on the persistent numbering
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
               if (channel == CH_XI) then
                  q_store(inum, istep) = cavity_fd%xi0(igrid)
               else
                  q_store(inum, istep) = cavity_fd%f(igrid)
               end if
               wb_store(inum, istep) = cavity_fd%wbranch(igrid)
               bc_store(inum, istep) = cavity_fd%offset(cavity_fd%anchor_of_row(igrid) + 1) &
                                        - cavity_fd%offset(cavity_fd%anchor_of_row(igrid))
            end do
         end block
      end do

      ! Only points present at every geometry carry weight
      allocate (usable(max_num))
      usable = all(have, dim=2)

      allocate (in_ref(max_num), source=.false.)
      do igrid = 1, ngrid
         inum = cavity%numbering(igrid)
         if (inum >= 1 .and. inum <= max_num) in_ref(inum) = .true.
      end do
      nphantom = count(usable .and. .not. in_ref)
      if (nphantom > 0) then
         call test_failed(error, "stencil point set does not match the reference cavity: "// &
                          to_string(nphantom)//" of "//to_string(count(usable))//" numbers common "// &
                          "to every displaced geometry are absent from the reference surface ("// &
                          to_string(ngrid)//" points), so the FD reference and the analytic sum "// &
                          "would cover different points")
         return
      end if

      do inum = 1, max_num
         if (.not. usable(inum)) cycle
         if (any(bc_store(inum, :) /= bc_store(inum, 1))) then
            call test_failed(error, "branch topology is not constant across the stencil for "// &
                             "point "//to_string(inum)//": branch counts "// &
                             to_string(bc_store(inum, 1))//", "//to_string(bc_store(inum, 2))// &
                             ", "//to_string(bc_store(inum, 3))//", "//to_string(bc_store(inum, 4))// &
                             ", "//to_string(bc_store(inum, 5))//", "//to_string(bc_store(inum, 6))// &
                             ", so its "//label//" is compared across different branch sets")
            return
         end if
      end do

      ncommon = count(usable)
      if (ncommon < 10) then
         call test_failed(error, "too few grid points survive the stencil ("// &
                          to_string(ncommon)//" of "//to_string(max_num)//")")
         return
      end if

      ! The branch weights must actually move along the coordinate
      dwb_max = 0.0_wp
      do inum = 1, max_num
         if (.not. usable(inum)) cycle
         dwb = (45.0_wp*(wb_store(inum, 3) - wb_store(inum, 4)) - &
                9.0_wp*(wb_store(inum, 2) - wb_store(inum, 5)) + &
                (wb_store(inum, 1) - wb_store(inum, 6)))/(60.0_wp*XI_FD_STEP)
         dwb_max = max(dwb_max, abs(dwb))
      end do
      if (dwb_max <= BRANCH_RESPONSE_THR) then
         call test_failed(error, "no branch weight responds along the FD coordinate (max |dw| "// &
                          to_string(dwb_max)//")")
         return
      end if

      ! Pair nearby values before summing to avoid subtracting large functionals
      num_deriv = 0.0_wp
      do inum = 1, max_num
         if (.not. usable(inum)) cycle
         num_deriv = num_deriv + wmap(inum)* &
                     (45.0_wp*(q_store(inum, 3) - q_store(inum, 4)) - &
                      9.0_wp*(q_store(inum, 2) - q_store(inum, 5)) + &
                      (q_store(inum, 1) - q_store(inum, 6)))/(60.0_wp*XI_FD_STEP)
      end do
      if (.not. ieee_is_finite(num_deriv)) then
         call test_failed(error, "branching "//label//" FD reference is not finite")
         return
      end if

      ! Analytic side: same restricted weights through the reverse path
      allocate (w_q(ngrid), source=0.0_wp)
      do igrid = 1, ngrid
         inum = cavity%numbering(igrid)
         if (inum >= 1 .and. inum <= max_num) then
            if (usable(inum)) w_q(igrid) = wmap(inum)
         end if
      end do
      call acc%init(ngrid)
      if (channel == CH_XI) then
         call acc%add_surface_weights(cav_error, w_xi=w_q)
      else
         call acc%add_surface_weights(cav_error, w_f=w_q)
      end if
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

      if (.not. ieee_is_finite(ana_deriv)) then
         call test_failed(error, "branching "//label//" analytic derivative is not finite")
         return
      end if
      diff = abs(ana_deriv - num_deriv)
      if (diff > XI_FD_ABS .and. diff > XI_FD_REL*abs(num_deriv)) then
         call test_failed(error, "branching "//label//" gradient disagrees with finite differences: "// &
                          "analytic "//to_string(ana_deriv)//" numeric "//to_string(num_deriv)// &
                          " (common points "//to_string(ncommon)//")")
         return
      end if

      ! Same stencil against the forward path, whose branch post-pass adds the
      ! branch-weight response to f1_rA
      call cavity%get_gradient(cav_error)
      if (allocated(cav_error)) then
         call test_failed(error, "forward gradient failed: "//cav_error%message)
         return
      end if
      ana_deriv = 0.0_wp
      do igrid = 1, ngrid
         if (channel == CH_XI) then
            ana_deriv = ana_deriv + w_q(igrid)*cavity%xi1_rA(FD_AXIS, FD_ATOM, igrid)
         else
            ana_deriv = ana_deriv + w_q(igrid)*cavity%f1_rA(FD_AXIS, FD_ATOM, igrid)
         end if
      end do

      if (.not. ieee_is_finite(ana_deriv)) then
         call test_failed(error, "forward "//label//" derivative is not finite under branching")
         return
      end if
      diff = abs(ana_deriv - num_deriv)
      if (diff > XI_FD_ABS .and. diff > XI_FD_REL*abs(num_deriv)) then
         call test_failed(error, "forward "//label//" derivative disagrees with finite differences under "// &
                          "branching: analytic "//to_string(ana_deriv)// &
                          " numeric "//to_string(num_deriv))
         return
      end if

   end subroutine run_branching_fd

   !> CPCM and PV gradients through the forward and reverse model paths
   subroutine test_model_forward_reverse(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      call check_model_forward_reverse(.false., error)
   end subroutine test_model_forward_reverse

   !> COSMO and PV gradients through the forward and reverse model paths
   subroutine test_cosmo_model_forward_reverse(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      call check_model_forward_reverse(.true., error)
   end subroutine test_cosmo_model_forward_reverse

   !> Compare model gradients and host responses from independent cavity paths
   !>
   !> @param[in] use_cosmo Select COSMO instead of CPCM
   !> @param[out] error Test failure
   subroutine check_model_forward_reverse(use_cosmo, error)
      !> Dielectric scaling choice
      logical, intent(in) :: use_cosmo
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      !> Dielectric constant and pressure of the probe model
      real(wp), parameter :: epsilon_r = 32.0_wp
      real(wp), parameter :: pressure = 0.75_wp
      !> Atomic charges driving the electrostatics
      real(wp), parameter :: qat_vals(*) = [0.20_wp, -0.15_wp, -0.05_wp]

      type(cavity_type_drop), allocatable :: cavity
      type(moist_context_type), target :: ctx
      type(model_continuum_type), target :: model_rev, model_fwd
      class(model_continuum_component_pcm_conductor), allocatable :: pcm_component
      type(model_continuum_component_pv) :: pv_component
      type(coupling_type), pointer :: coupling
      !> Host part of the gradient phase from each path
      type(response_type) :: response_rev, response_fwd
      !> Copies of the potential adjoint items published by each path
      type(potential_adjoint_response_type), allocatable :: charge_rev, charge_fwd
      type(structure_type) :: mol
      type(mctc_error), allocatable :: err
      !> Borrowed model cavities
      class(cavity_type), pointer :: cav_rev, cav_fwd

      real(wp), allocatable :: grad_rev(:, :), grad_fwd(:, :)
      real(wp) :: diff, scale
      !> Energy of the reverse model, solved to stage the gradient phase
      real(wp) :: energy
      integer :: nat, iatom, iaxis

      call fixture_geometry(.false., mol)
      call build_cavity(cavity, ctx, mol, .false., error)
      if (allocated(error)) return
      nat = mol%nat

      if (use_cosmo) then
         allocate (model_continuum_component_cosmo :: pcm_component)
      else
         allocate (model_continuum_component_cpcm :: pcm_component)
      end if
      select type (pcm_component)
      type is (model_continuum_component_cosmo)
         call new_component_cosmo(pcm_component, epsilon=epsilon_r, error=err, &
                                  param=moist_pcm_parameters_type(solver=solver_type%dense), ctx=ctx)
      type is (model_continuum_component_cpcm)
         call new_component_cpcm(pcm_component, epsilon=epsilon_r, error=err, &
                                 param=moist_pcm_parameters_type(solver=solver_type%dense), ctx=ctx)
      end select
      if (allocated(err)) then
         call test_failed(error, "CPCM construction failed: "//err%message)
         return
      end if
      call new_component_pv(pv_component, pressure)

      call build_model(model_rev, cavity, ctx, pcm_component, pv_component, mol, error)
      if (allocated(error)) return
      call build_model(model_fwd, cavity, ctx, pcm_component, pv_component, mol, error)
      if (allocated(error)) return
      call model_fwd%use_forward_gradient(.true.)
      cav_rev => model_rev%cavity
      cav_fwd => model_fwd%cavity

      ! Point-charge potential trace plus its total position weight
      call stage_model_point_charge_energy(error, model_rev, qat_vals, mol, coupling)
      if (allocated(error)) return
      energy = 0.0_wp
      call model_rev%get_energy(coupling, energy, err)
      if (allocated(err)) then
         call test_failed(error, "reverse-mode model energy failed: "//err%message)
         return
      end if
      call model_rev%prepare_gradient(coupling, err)
      if (allocated(err)) then
         call test_failed(error, "gradient staging failed: "//err%message)
         return
      end if
      call fill_missing_with_zeros(cav_rev, coupling)
      call fill_point_charge_field(cav_rev, coupling, qat_vals, mol)

      allocate (grad_rev(3, nat), source=0.0_wp)
      allocate (grad_fwd(3, nat), source=0.0_wp)

      call model_rev%get_gradient(coupling, response_rev, grad_rev, err)
      if (allocated(err)) then
         call test_failed(error, "reverse-mode model gradient failed: "//err%message)
         return
      end if
      call check(error,.not. allocated(cav_rev%xyz1_rA), &
                 more="reverse model built a surface nuclear Jacobian")
      if (allocated(error)) return
      call stage_model_point_charge_energy(error, model_fwd, qat_vals, mol, coupling)
      if (allocated(error)) return
      call model_fwd%prepare_gradient(coupling, err)
      call fill_point_charge_field(cav_fwd, coupling, qat_vals, mol)
      call model_fwd%get_gradient(coupling, response_fwd, grad_fwd, err)
      if (allocated(err)) then
         call test_failed(error, "forward model gradient failed: "//err%message)
         return
      end if

      call check(error, allocated(cav_fwd%xyz1_rA), &
                 more="forward model did not build a surface nuclear Jacobian")
      if (allocated(error)) return

      ! Both paths publish the host part of the phase
      call copy_potential_adjoint(response_rev, charge_rev)
      call copy_potential_adjoint(response_fwd, charge_fwd)
      if (.not. allocated(charge_rev) .or. .not. allocated(charge_fwd)) then
         call test_failed(error, "model gradient published no potential adjoint")
         return
      end if
      call check(error, maxval(abs(charge_rev%w_phi - charge_fwd%w_phi)), 0.0_wp, thr=EQ_ABS, &
         & more="reverse and forward gradient paths published different potential adjoints")
      if (allocated(error)) return
      call check(error, maxval(abs(charge_rev%w_phi)) > 0.0_wp, &
         & more="model gradient published a vanishing potential adjoint")
      if (allocated(error)) return

      if (maxval(abs(grad_fwd)) <= VACUITY_THR) then
         call test_failed(error, "model gradient is vacuous")
         return
      end if

      do iatom = 1, nat
         do iaxis = 1, 3
            diff = abs(grad_rev(iaxis, iatom) - grad_fwd(iaxis, iatom))
            scale = abs(grad_fwd(iaxis, iatom))
            ! Negated pass test, so a NaN on either path fails
            if (.not. (diff <= AB_ABS .or. diff <= AB_REL*scale)) then
               call test_failed(error, "model gradient differs between paths at atom "// &
                                to_string(iatom)//" axis "//to_string(iaxis)// &
                                ": reverse "//to_string(grad_rev(iaxis, iatom))// &
                                " forward "//to_string(grad_fwd(iaxis, iatom)))
               return
            end if
         end do
      end do

   end subroutine check_model_forward_reverse

   !> Assemble a conductor PCM + PV model on a DROP cavity
   !>
   !> @param[out]   model  Model to build
   !> @param[in]    cavity Cavity template copied into the model
   !> @param[in]    ctx    Run context owned by the caller
   !> @param[in]    pcmc   Conductor PCM component template
   !> @param[in]    pvc    PV component template
   !> @param[in]    mol    Molecular structure
   !> @param[out]   error  Error handle
   subroutine build_model(model, cavity, ctx, pcmc, pvc, mol, error)
      !> Model to build
      type(model_continuum_type), intent(out) :: model
      !> Cavity template copied into the model
      type(cavity_type_drop), intent(in) :: cavity
      !> Run context owned by the caller
      type(moist_context_type), intent(in), target :: ctx
      !> Component templates
      class(model_continuum_component_pcm_conductor), intent(in) :: pcmc
      type(model_continuum_component_pv), intent(in) :: pvc
      !> Molecular structure
      type(structure_type), intent(in) :: mol
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(mctc_error), allocatable :: err

      call new_continuum_model(model, ctx, cavity, err)
      if (.not. allocated(err)) call model%add_component(pcmc, err)
      if (.not. allocated(err)) call model%add_component(pvc, err)
      if (.not. allocated(err)) call model%update(mol, err)
      if (allocated(err)) then
         call test_failed(error, "model setup failed: "//err%message)
         return
      end if

   end subroutine build_model

   !> The reverse gradient is bit-identical on one and on four threads
   subroutine test_thread_count_invariance(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      !> Team size compared against one thread
      integer, parameter :: TEAM = 4
      type(cavity_type_drop), allocatable :: cavity
      type(moist_context_type), target :: ctx, ctx_serial, ctx_team
      type(cavity_adjoint_type) :: acc
      type(mctc_error), allocatable :: cav_error
      real(wp), allocatable :: grad_serial(:, :), grad_team(:, :)
      integer :: ich
      logical :: in_parallel

      ! Inside test-drive's suite-level team the library region would get one thread
      in_parallel = .false.
!$    in_parallel = omp_in_parallel()
      call check(error, .not. in_parallel, &
         & "thread-count check must run as a selected test, not inside the suite-level parallel run")
      if (allocated(error)) return

      ! The branching cross: several tiles per thread and multi-row anchors
      call build_fixture(cavity, ctx, .true., error)
      if (allocated(error)) return
      if (cavity%ngrid <= 2*TEAM*surface_tile_width) then
         call test_failed(error, "fixture too small to spread over "//to_string(TEAM)// &
                          " threads: "//to_string(cavity%ngrid)//" points")
         return
      end if

      call acc%init(cavity%ngrid)
      do ich = 1, NCHAN
         call apply_channel(acc, ich, cavity, error)
         if (allocated(error)) return
      end do

      call new_context(ctx_serial, nthreads=1, verbosity=0)
      call new_context(ctx_team, nthreads=TEAM, verbosity=0)
      allocate (grad_serial(3, cavity%nsph), source=0.0_wp)
      allocate (grad_team(3, cavity%nsph), source=0.0_wp)

      cavity%ctx => ctx_serial
      call cavity%get_surface_gradient(acc, grad_serial, cav_error)
      if (.not. allocated(cav_error)) then
         cavity%ctx => ctx_team
         call cavity%get_surface_gradient(acc, grad_team, cav_error)
      end if
      cavity%ctx => ctx
      if (allocated(cav_error)) then
         call test_failed(error, "reverse gradient failed: "//cav_error%message)
         return
      end if

      if (maxval(abs(grad_serial)) <= VACUITY_THR) then
         call test_failed(error, "reverse gradient is vacuous")
         return
      end if
      if (any(grad_serial /= grad_team)) then
         call test_failed(error, "reverse gradient depends on the thread count: max deviation "// &
                          to_string(maxval(abs(grad_serial - grad_team)))// &
                          " between 1 and "//to_string(TEAM)//" threads")
      end if

   end subroutine test_thread_count_invariance

   !> The internal isodensity cavity, its surface LSF weights and its reverse
   !> gradient are bit-identical on one and on several threads, and on a repeated call
   !>
   !> - the per-thread workers borrow the model's basis and density, in the
   !>   projection, both reverse loops, the forward Jacobian, the curvature and
   !>   the atom areas
   !> - the atom areas, a sum over the points in thread order, to round-off only
   !> - must run as a selected test, see [[test_thread_count_invariance]]
   subroutine test_isodensity_thread_count_invariance(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: TEAM = 4
      type(cavity_type_drop), allocatable :: cav_serial, cav_team
      type(moist_context_type), target :: ctx_serial, ctx_team
      type(cavity_adjoint_type) :: acc
      type(mctc_error), allocatable :: cav_error
      real(wp), allocatable :: grad(:, :, :), w0(:, :), w1(:, :, :), w2(:, :, :, :)
      integer :: ich, irun, nteam
      logical :: in_parallel

      in_parallel = .false.
!$    in_parallel = omp_in_parallel()
      call check(error, .not. in_parallel, &
         & "thread-count check must run as a selected test, not inside the suite-level parallel run")
      if (allocated(error)) return

      ! Without a full team the comparison below would pass without testing anything
      nteam = 1
      !$omp parallel num_threads(TEAM) default(shared)
      !$omp single
!$    nteam = omp_get_num_threads()
      !$omp end single
      !$omp end parallel
      if (nteam /= TEAM) then
         call skip_test(error, "a team of "//to_string(TEAM)//" threads is not available, got "// &
                        to_string(nteam))
         return
      end if

      call new_context(ctx_serial, nthreads=1, verbosity=0)
      call new_context(ctx_team, nthreads=TEAM, verbosity=0)
      call build_internal_iso_cavity(cav_serial, ctx_serial, error)
      if (allocated(error)) return
      call build_internal_iso_cavity(cav_team, ctx_team, error)
      if (allocated(error)) return
      if (cav_serial%ngrid <= 2*TEAM*surface_tile_width) then
         call test_failed(error, "fixture too small to spread over "//to_string(TEAM)// &
                          " threads: "//to_string(cav_serial%ngrid)//" points")
         return
      end if
      if (cav_team%ngrid /= cav_serial%ngrid) then
         call test_failed(error, "the isodensity grid depends on the thread count")
         return
      end if
      if (any(cav_team%xyz /= cav_serial%xyz) .or. any(cav_team%a /= cav_serial%a)) then
         call test_failed(error, "the isodensity surface depends on the thread count")
         return
      end if

      call acc%init(cav_serial%ngrid)
      do ich = 1, NCHAN
         call apply_channel(acc, ich, cav_serial, error)
         if (allocated(error)) return
      end do

      ! Run 1 serial, runs 2 and 3 on the team: thread count and repeated call
      allocate (grad(3, cav_serial%nsph, 3), source=0.0_wp)
      allocate (w0(cav_serial%ngrid, 3), w1(3, cav_serial%ngrid, 3), w2(3, 3, cav_serial%ngrid, 3))
      do irun = 1, 3
         if (irun == 1) then
            call cav_serial%contract_surface_lsf_weights(acc, w0(:, irun), w1(:, :, irun), &
                                                         w2(:, :, :, irun), cav_error)
            if (.not. allocated(cav_error)) &
               call cav_serial%get_surface_gradient(acc, grad(:, :, irun), cav_error)
         else
            call cav_team%contract_surface_lsf_weights(acc, w0(:, irun), w1(:, :, irun), &
                                                       w2(:, :, :, irun), cav_error)
            if (.not. allocated(cav_error)) &
               call cav_team%get_surface_gradient(acc, grad(:, :, irun), cav_error)
         end if
         if (allocated(cav_error)) then
            call test_failed(error, "isodensity reverse pass failed: "//cav_error%message)
            return
         end if
      end do

      if (maxval(abs(grad(:, :, 1))) <= VACUITY_THR .or. maxval(abs(w0(:, 1))) <= VACUITY_THR) then
         call test_failed(error, "isodensity reverse pass is vacuous")
         return
      end if
      do irun = 2, 3
         call check(error, all(w0(:, irun) == w0(:, 1)) .and. all(w1(:, :, irun) == w1(:, :, 1)) &
                    .and. all(w2(:, :, :, irun) == w2(:, :, :, 1)), &
                    "isodensity LSF weights depend on the thread count or the call, run "// &
                    to_string(irun))
         if (allocated(error)) return
         call check(error, all(grad(:, :, irun) == grad(:, :, 1)), &
                    "isodensity reverse gradient depends on the thread count or the call, run "// &
                    to_string(irun))
         if (allocated(error)) return
      end do

      ! The curvature and atom-area loops of the build have their own workers
      if (.not. (allocated(cav_serial%KM) .and. allocated(cav_serial%KG) &
                 .and. allocated(cav_serial%asph_partition))) then
         call test_failed(error, "isodensity curvature or atom areas missing")
         return
      end if
      if (maxval(abs(cav_serial%KM)) <= VACUITY_THR .or. &
          maxval(abs(cav_serial%asph_partition)) <= VACUITY_THR) then
         call test_failed(error, "isodensity curvature or atom areas are vacuous")
         return
      end if
      call check(error, all(cav_team%KM == cav_serial%KM) .and. all(cav_team%KG == cav_serial%KG), &
                 "isodensity curvature depends on the thread count")
      if (allocated(error)) return
      ! The per-atom sums are reduced in thread order, so they agree to round-off only
      call check(error, maxval(abs(cav_team%asph_partition - cav_serial%asph_partition)) &
                 <= 1.0e-12_wp*maxval(abs(cav_serial%asph_partition)), &
                 "isodensity atom areas depend on the thread count")
      if (allocated(error)) return

      ! The forward Jacobian loop as well
      call cav_serial%compute_gradient_drop(cav_error)
      if (.not. allocated(cav_error)) call cav_team%compute_gradient_drop(cav_error)
      if (allocated(cav_error)) then
         call test_failed(error, "isodensity forward Jacobian failed: "//cav_error%message)
         return
      end if
      if (maxval(abs(cav_serial%xyz1_rA)) <= VACUITY_THR) then
         call test_failed(error, "isodensity forward Jacobian is vacuous")
         return
      end if
      call check(error, all(cav_team%xyz1_rA == cav_serial%xyz1_rA) &
                 .and. all(cav_team%a_i1_rA == cav_serial%a_i1_rA), &
                 "isodensity forward Jacobian depends on the thread count")
      if (allocated(error)) return

   end subroutine test_isodensity_thread_count_invariance

   !> Internal isodensity cavity on two d-shell atoms, built through the coupling
   !>
   !> - with the curvature and atom-area diagnostics
   !>
   !> @param[out]   cavity Built cavity
   !> @param[inout] ctx    Run context borrowed by the cavity, team size set
   !> @param[out]   error  Test failure
   subroutine build_internal_iso_cavity(cavity, ctx, error)
      !> Built cavity
      type(cavity_type_drop), allocatable, intent(out) :: cavity
      !> Run context borrowed by the cavity
      type(moist_context_type), target, intent(inout) :: ctx
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      type(moist_cavity_drop_lsf_isodensity_internal_type) :: lsf
      type(moist_cavity_drop_lsf_isodensity_param_type) :: param
      type(structure_type) :: mol
      class(radius_type), allocatable :: radii
      type(coupling_type), target :: coupling
      type(mctc_error), allocatable :: err
      real(wp) :: dcart(12, 12)
      integer :: i

      call new(mol, [54, 54], reshape([0.0_wp, 0.0_wp, 0.0_wp, 0.5_wp, 0.0_wp, 5.0_wp], [3, 2]))
      param = moist_cavity_drop_lsf_isodensity_param_type(core_fill=.false.)
      param%rho_iso = 0.4_wp
      param%scale = 1.0_wp/param%rho_iso
      dcart = 0.0_wp
      do i = 1, 12
         dcart(i, i) = 1.0_wp
      end do
      call lsf%new(sh_atom=[1, 2], sh_l=[2, 2], sh_nprim=[1, 1], exps=[0.5_wp, 0.4_wp], &
                   coeffs=[1.0_wp, 1.0_wp], error=err, param=param)
      if (.not. allocated(err)) call new_radii_custom_atoms([4.0_wp, 4.0_wp], radii, err)
      if (.not. allocated(err)) then
         allocate (cavity)
         call new_cavity_drop(cavity, radius_model=radii, lsf_model=lsf, error=err, &
                              param=moist_cavity_drop_parameters_type(num_leb=302, do_curvature=.true., &
                                                                      do_atomarea=.true.), ctx=ctx)
      end if
      if (.not. allocated(err)) call cavity%update(mol, err)
      if (.not. allocated(err)) call cavity%prepare_build(coupling, err)
      if (.not. allocated(err)) then
         do while (coupling%next())
            select type (request => coupling%request())
            type is (density_matrix_request_type)
               call coupling%answer("dcart", dcart, err)
            type is (nuclear_charge_request_type)
               call coupling%answer("z", [54.0_wp, 54.0_wp], err)
            end select
            if (allocated(err)) exit
         end do
      end if
      if (.not. allocated(err)) call cavity%build(coupling, err)
      if (allocated(err)) then
         call test_failed(error, "failed to build isodensity cavity: "//err%message)
         return
      end if
      if (cavity%ngrid == 0) call test_failed(error, "empty isodensity grid")

   end subroutine build_internal_iso_cavity

   !> A mis-shaped gradient accumulator must be rejected, not silently ignored
   subroutine test_shape_guard(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(cavity_type_drop), allocatable :: cavity
      type(moist_context_type), target :: ctx
      type(cavity_adjoint_type) :: acc
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
      !> Grid extent
      integer :: ngrid

      call build_fixture(cavity, ctx, branching, error)
      if (allocated(error)) return

      ngrid = cavity%ngrid

      ! Without real branching the branch-weight reverse pass and its explicit owner
      ! term are dead code here, so the branching variant must prove it branched
      if (branching) then
         if (.not. allocated(cavity%offset)) then
            call test_failed(error, "branching fixture has no row layout")
            return
         end if
         if (.not. any(cavity%offset(2:cavity%anchor%ngrid + 1) - cavity%offset(1:cavity%anchor%ngrid) > 1)) then
            call test_failed(error, "branching fixture did not branch; "// &
                             "the branch-weight channel would go untested (ngrid "// &
                             to_string(ngrid)//", max branch_count "// &
                             to_string(maxval(cavity%offset(2:cavity%anchor%ngrid + 1) &
                                              - cavity%offset(1:cavity%anchor%ngrid)))//")")
            return
         end if
      end if

      call compare_paths(cavity, channels, label, error)

   end subroutine run_equivalence

   !> Largest response of a row weight to the gate of a partly open sibling row
   !>
   !> @param[in] cavity Built cavity
   !> @return    live   `max |d w_m / d g_l|`, `l /= m` rows of one anchor, `0 < g_l < 1`
   function sibling_gate_response(cavity) result(live)
      type(cavity_type_drop), intent(in) :: cavity
      real(wp) :: live

      integer :: ia, r0, n, m, l
      real(wp), allocatable :: w(:), jac_phi(:, :), jac_gate(:, :)

      live = 0.0_wp
      do ia = 1, cavity%anchor%ngrid
         r0 = cavity%offset(ia)
         n = cavity%offset(ia + 1) - r0
         if (n < 2) cycle
         allocate (w(n), jac_phi(n, n), jac_gate(n, n))
         call cavity%branch_weight%weights_jac(cavity%phi0(r0 + 1:r0 + n), cavity%w_f0(r0 + 1:r0 + n), &
                                               w, jac_phi, jac_gate)
         do l = 1, n
            if (.not. (cavity%w_f0(r0 + l) > 0.0_wp .and. cavity%w_f0(r0 + l) < 1.0_wp)) cycle
            do m = 1, n
               if (m /= l) live = max(live, abs(jac_gate(m, l)))
            end do
         end do
         deallocate (w, jac_phi, jac_gate)
      end do
   end function sibling_gate_response

   !> Every channel with many rows inside the reach band, whose direct response
   !> to the anchor (the owner) the reverse path adds
   subroutine test_reach_gate_channels(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(structure_type) :: mol
      type(cavity_type_drop), allocatable :: cavity
      type(moist_context_type), target :: ctx
      integer :: i, nband
      real(wp) :: f

      call get_test_cube(2.5093779394_wp, mol)
      call build_switch_fixture(cavity, ctx, mol, 5.5_wp, 0.2_wp, error, reach_on=0.2_wp, reach_off=1.2_wp)
      if (allocated(error)) return
      nband = 0
      do i = 1, cavity%ngrid
         f = cavity%f_reach%f0(cavity%pdist(i))
         if (f > 0.0_wp .and. f < 1.0_wp) nband = nband + 1
      end do
      if (nband == 0) then
         call test_failed(error, "fixture has no row inside the reach band")
         return
      end if
      call compare_paths(cavity, [CH_XI, CH_F, CH_A, CH_W, CH_XYZ, CH_N, CH_K1, CH_K2], "reach gate", error)

   end subroutine test_reach_gate_channels

   !> Every channel on anchors that share between two forward-facing points under
   !> the default switches, on and off the mirror geometry
   subroutine test_forward_pair_channels(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(structure_type) :: mol
      type(cavity_type_drop), allocatable :: cavity
      type(moist_context_type), target :: ctx
      real(wp) :: dz
      integer :: i

      do i = 1, 2
         dz = merge(0.0_wp, 0.06_wp, i == 1)
         call get_test_forward_pair(dz, mol)
         call build_switch_fixture(cavity, ctx, mol, 2.0_wp, 0.2_wp, error)
         if (allocated(error)) return
         if (.not. any(cavity%offset(2:cavity%anchor%ngrid + 1) - cavity%offset(1:cavity%anchor%ngrid) > 1)) then
            call test_failed(error, "forward pair did not branch")
            return
         end if
         call compare_paths(cavity, [CH_XI, CH_F, CH_A, CH_W, CH_XYZ, CH_N, CH_K1, CH_K2], "forward pair", error)
         if (allocated(error)) return
      end do

   end subroutine test_forward_pair_channels

   !> Every channel on the five-atom fixture whose anchors hold a partly open
   !> point next to a forward-facing row, under the default switches
   subroutine test_sibling_gate_channels(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(structure_type) :: mol
      type(cavity_type_drop), allocatable :: cavity
      type(moist_context_type), target :: ctx

      call get_test_sibling_gate(mol)
      call build_switch_fixture(cavity, ctx, mol, 1.5_wp, 0.2_wp, error)
      if (allocated(error)) return
      if (.not. sibling_gate_response(cavity) > 1.0e-2_wp) then
         call test_failed(error, "fixture has no partly open row that moves a sibling weight")
         return
      end if
      call compare_paths(cavity, [CH_XI, CH_F, CH_A, CH_W, CH_XYZ, CH_N, CH_K1, CH_K2], "sibling gate", error)

   end subroutine test_sibling_gate_channels

   !> Every channel and the total area gradient where a switching neighbour of an
   !> anchor is not LSF-active at the projected point
   !>
   !> - SvdW fade band of 0.33 Bohr past the sphere, shorter than the reach of
   !>   the switching function: the neighbour has no level set partial at the
   !>   point, yet the switching value of the row moves with it
   !> - reach gate closed at 0.32 Bohr, inside the fade band: every owner stays
   !>   LSF-active, so only the neighbour path is under test
   !> - finite differences: the short band makes the level set stiff, the 4-point
   !>   stencil at 1e-4 Bohr is off by 3e-5 (falling as h^4); the screened
   !>   neighbours carry 1.9e-3, which a gradient without them is off by
   subroutine test_forward_screened_neighbour(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(structure_type) :: mol
      type(cavity_type_drop), allocatable :: cavity
      type(moist_context_type), target :: ctx
      integer :: n_owner, n_neighbour
      real(wp), allocatable :: area_nb(:, :)

      call fixture_geometry(.false., mol)
      call build_screened_fixture(cavity, ctx, mol, NB_FADE, error, NB_REACH_ON, NB_REACH_OFF)
      if (allocated(error)) return
      call count_screened(cavity, n_owner, n_neighbour, error)
      if (allocated(error)) return
      if (n_neighbour == 0 .or. n_owner /= 0) then
         call test_failed(error, "fixture needs rows with a screened switching neighbour and no screened owner: "// &
                          to_string(n_neighbour)//", "//to_string(n_owner))
         return
      end if
      call compare_paths(cavity, [CH_XI, CH_F, CH_A, CH_W, CH_XYZ, CH_N, CH_K1, CH_K2], "screened neighbour", error)
      if (allocated(error)) return
      ! The screened neighbours must carry more area gradient than the bound hides
      call count_screened(cavity, n_owner, n_neighbour, error, area_nb)
      if (allocated(error)) return
      if (.not. maxval(abs(area_nb)) > NB_AREA_MIN) then
         call test_failed(error, "screened neighbours carry too little area gradient for the finite-difference "// &
                          "bound: "//to_string(maxval(abs(area_nb))))
         return
      end if
      call check_area_fd(cavity, mol, NB_FADE, NB_FD_TOL, error, NB_REACH_ON, NB_REACH_OFF)

   end subroutine test_forward_screened_neighbour

   !> Every channel and the total area gradient where the owner of an anchor is
   !> not LSF-active at the projected point
   !>
   !> - SvdW fade band of 0.16 Bohr past the sphere: an anchor just inside a
   !>   neighbour sphere projects onto that sphere, out of reach of its owner,
   !>   whose motion still carries the anchor
   !> - finite differences: off by 2e-6 at 1e-4 Bohr; a gradient without the
   !>   owner seed is off by 0.4
   subroutine test_forward_screened_owner(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(structure_type) :: mol
      type(cavity_type_drop), allocatable :: cavity
      type(moist_context_type), target :: ctx
      integer :: n_owner, n_neighbour

      call fixture_geometry(.false., mol)
      call build_screened_fixture(cavity, ctx, mol, 0.2_wp, error)
      if (allocated(error)) return
      call count_screened(cavity, n_owner, n_neighbour, error)
      if (allocated(error)) return
      if (n_owner == 0) then
         call test_failed(error, "fixture has no row whose owner is screened at the projected point")
         return
      end if
      call compare_paths(cavity, [CH_XI, CH_F, CH_A, CH_W, CH_XYZ, CH_N, CH_K1, CH_K2], "screened owner", error)
      if (allocated(error)) return
      call check_area_fd(cavity, mol, 0.2_wp, 1.0e-5_wp, error)

   end subroutine test_forward_screened_owner

   !> Every channel on a callback isodensity level set, which has no active atom
   !> and no nuclear partial: the anchors still move with their owners
   !>
   !> - the full forward gradient equals the anchor-only one
   subroutine test_forward_isodensity_owner(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(structure_type) :: mol
      type(cavity_type_drop), allocatable :: cavity
      type(moist_context_type), target :: ctx
      type(moist_cavity_drop_lsf_isodensity_callback_type) :: lsf
      class(radius_type), allocatable :: radii
      type(mctc_error), allocatable :: cav_error
      real(wp), allocatable :: xyz1(:, :, :, :), a1(:, :, :), area1(:, :)

      call new(mol, [10, 10], density_centers)
      call new_context(ctx, nthreads=0, verbosity=0)
      call new_radii_custom_atoms([4.0_wp, 4.0_wp], radii, cav_error)
      if (.not. allocated(cav_error)) then
         call lsf%new(callback_ptr=c_funloc(two_center_density), context=c_null_ptr, &
                      param=moist_cavity_drop_lsf_isodensity_param_type(core_fill=.false.))
         allocate (cavity)
         call new_cavity_drop(cavity, radius_model=radii, lsf_model=lsf, error=cav_error, &
                              param=moist_cavity_drop_parameters_type(num_leb=302, do_fine=.true.), ctx=ctx)
      end if
      if (.not. allocated(cav_error)) call cavity%update(mol, error=cav_error)
      if (allocated(cav_error)) then
         call test_failed(error, "failed to build cavity: "//cav_error%message)
         return
      end if
      if (cavity%ngrid == 0) then
         call test_failed(error, "empty grid")
         return
      end if

      call compare_paths(cavity, [CH_XI, CH_F, CH_A, CH_W, CH_XYZ, CH_N, CH_K1, CH_K2], "callback isodensity", error)
      if (allocated(error)) return
      xyz1 = cavity%xyz1_rA
      a1 = cavity%a_i1_rA
      area1 = cavity%A_tot1_rA
      if (.not. maxval(abs(xyz1)) > 0.1_wp) then
         call test_failed(error, "the projected points do not move with the atoms")
         return
      end if
      call cavity%compute_anchor_gradient(cav_error)
      if (allocated(cav_error)) then
         call test_failed(error, "anchor gradient failed: "//cav_error%message)
         return
      end if
      call check(error, maxval(abs(cavity%xyz1_rA - xyz1)) <= 1.0e-14_wp &
                 .and. maxval(abs(cavity%a_i1_rA - a1)) <= 1.0e-14_wp &
                 .and. maxval(abs(cavity%A_tot1_rA - area1)) <= 1.0e-12_wp, &
                 "full and anchor-only forward gradients differ on a level set without nuclear partials")

   end subroutine test_forward_isodensity_owner

   !> Anchor-only gradient on a level set with nuclear partials: the field is
   !> frozen, so every point moves with its owner alone
   !>
   !> - control: the full gradient moves points with non-owner atoms
   subroutine test_forward_frozen_field(error)
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(cavity_type_drop), allocatable :: cavity
      type(moist_context_type), target :: ctx
      type(mctc_error), allocatable :: cav_error
      real(wp) :: full_other, full_owner, frozen_other, frozen_owner

      call build_fixture(cavity, ctx, .false., error)
      if (allocated(error)) return
      call cavity%get_gradient(cav_error)
      if (.not. allocated(cav_error)) then
         call point_motion(cavity, full_other, full_owner)
         call cavity%compute_anchor_gradient(cav_error)
      end if
      if (allocated(cav_error)) then
         call test_failed(error, "gradient failed: "//cav_error%message)
         return
      end if
      call point_motion(cavity, frozen_other, frozen_owner)

      call check(error, full_other > 1.0e-3_wp, &
                 "control: the full gradient moves no point with a non-owner atom")
      if (allocated(error)) return
      call check(error, frozen_owner > 1.0e-3_wp, "the anchor-only points do not move with their owners")
      if (allocated(error)) return
      call check(error, frozen_other == 0.0_wp, &
                 "the anchor-only gradient moves a point with a non-owner atom: "//to_string(frozen_other))

   end subroutine test_forward_frozen_field

   !> Largest point response to the owner and to any other atom
   !>
   !> @param[in]  cavity Cavity holding `xyz1_rA`
   !> @param[out] other  `max |d r_i / d R_A|`, `A` not the owner of `i`
   !> @param[out] owner  `max |d r_i / d R_A|`, `A` the owner of `i`
   subroutine point_motion(cavity, other, owner)
      type(cavity_type_drop), intent(in) :: cavity
      real(wp), intent(out) :: other
      real(wp), intent(out) :: owner

      integer :: igrid, iatom
      real(wp) :: m

      other = 0.0_wp
      owner = 0.0_wp
      do igrid = 1, cavity%ngrid
         do iatom = 1, cavity%nsph
            m = maxval(abs(cavity%xyz1_rA(:, :, iatom, igrid)))
            if (iatom == cavity%owner(igrid)) then
               owner = max(owner, m)
            else
               other = max(other, m)
            end if
         end do
      end do
   end subroutine point_motion

   !> Density of two Gaussians `exp(-|r - c|^2/2)` at `density_centers`, fixed in space
   !>
   !> @param[in]  context Unused
   !> @param[in]  point   Evaluation point (Bohr)
   !> @param[out] rho     Density
   !> @param[out] drho    Gradient
   !> @param[in]  d2rho   Hessian (3, 3), or NULL
   !> @param[in]  d3rho   Third derivative (3, 3, 3), or NULL
   !> @return     status  Zero
   function two_center_density(context, point, rho, drho, d2rho, d3rho) result(status) bind(C)
      type(c_ptr), value :: context
      real(c_double), intent(in) :: point(3)
      real(c_double), intent(out) :: rho
      real(c_double), intent(out) :: drho(3)
      type(c_ptr), value :: d2rho
      type(c_ptr), value :: d3rho
      integer(c_int) :: status

      real(c_double), pointer :: hptr(:, :), tptr(:, :, :)
      real(wp) :: g, x(3), delta(3, 3)
      integer :: ic, i, j, l

      status = 0_c_int
      if (c_associated(context)) status = 1_c_int
      delta = 0.0_wp
      do i = 1, 3
         delta(i, i) = 1.0_wp
      end do
      rho = 0.0_wp
      drho = 0.0_wp
      if (c_associated(d2rho)) then
         call c_f_pointer(d2rho, hptr, [3, 3])
         hptr = 0.0_wp
      end if
      if (c_associated(d3rho)) then
         call c_f_pointer(d3rho, tptr, [3, 3, 3])
         tptr = 0.0_wp
      end if
      do ic = 1, size(density_centers, 2)
         x = point - density_centers(:, ic)
         g = exp(-0.5_wp*sum(x**2))
         rho = rho + g
         drho = drho - x*g
         if (c_associated(d2rho)) then
            do j = 1, 3
               do i = 1, 3
                  hptr(i, j) = hptr(i, j) + (x(i)*x(j) - delta(i, j))*g
               end do
            end do
         end if
         if (c_associated(d3rho)) then
            do l = 1, 3
               do j = 1, 3
                  do i = 1, 3
                     tptr(i, j, l) = tptr(i, j, l) + (-x(i)*x(j)*x(l) + x(i)*delta(j, l) + x(j)*delta(i, l) &
                                                      + x(l)*delta(i, j))*g
                  end do
               end do
            end do
         end if
      end do
   end function two_center_density

   !> Cavity of the screening fixtures: 302 points, level 4, SvdW with a short fade band
   !>
   !> @param[out]   cavity     Built cavity
   !> @param[inout] ctx        Run context borrowed by the cavity
   !> @param[in]    mol        Structure
   !> @param[in]    fade_start Fade onset of the level set in e-folds; the band is half as wide
   !> @param[out]   error      Error handle
   !> @param[in]    reach_on   Open edge of the reach gate (default: the library's)
   !> @param[in]    reach_off  Closed edge of the reach gate (default: the library's)
   subroutine build_screened_fixture(cavity, ctx, mol, fade_start, error, reach_on, reach_off)
      type(cavity_type_drop), allocatable, intent(out) :: cavity
      type(moist_context_type), target, intent(inout) :: ctx
      type(structure_type), intent(in) :: mol
      real(wp), intent(in) :: fade_start
      type(error_type), allocatable, intent(out) :: error
      real(wp), intent(in), optional :: reach_on, reach_off

      type(moist_cavity_drop_lsf_svdw_type) :: svdw_template
      type(moist_cavity_drop_parameters_type) :: param
      type(mctc_error), allocatable :: cav_error

      call new_context(ctx, nthreads=0, verbosity=0)
      param = moist_cavity_drop_parameters_type(num_leb=302, tolerance=PROJ_TOL, &
                                                proj_level=4, do_fine=.true.)
      if (present(reach_on)) param%branch%reach_on = reach_on
      if (present(reach_off)) param%branch%reach_off = reach_off
      call svdw_template%new(param=moist_cavity_drop_lsf_svdw_param_type(blend_k=5.5_wp, blend_3b=3.0_wp, &
                                                                         fade_start=fade_start, &
                                                                         fade_width=0.5_wp*fade_start))
      allocate (cavity)
      call new_cavity_drop(cavity, radius_model=default_cpcm_radii(), lsf_model=svdw_template, &
                           error=cav_error, param=param, ctx=ctx)
      if (.not. allocated(cav_error)) call cavity%update(mol, error=cav_error)
      if (allocated(cav_error)) call test_failed(error, "failed to build cavity: "//cav_error%message)
   end subroutine build_screened_fixture

   !> Forward total area gradient of a screening fixture against 4-point central differences
   !>
   !> @param[inout] cavity     Built cavity with its forward gradient
   !> @param[in]    mol        Reference structure
   !> @param[in]    fade_start Fade onset of the fixture
   !> @param[in]    tol        Absolute bound
   !> @param[out]   error      Error handle
   !> @param[in]    reach_on   Open edge of the fixture's reach gate (optional)
   !> @param[in]    reach_off  Closed edge of the fixture's reach gate (optional)
   subroutine check_area_fd(cavity, mol, fade_start, tol, error, reach_on, reach_off)
      type(cavity_type_drop), intent(inout) :: cavity
      type(structure_type), intent(in) :: mol
      real(wp), intent(in) :: fade_start, tol
      type(error_type), allocatable, intent(out) :: error
      real(wp), intent(in), optional :: reach_on, reach_off

      real(wp), parameter :: h = 1.0e-4_wp
      real(wp), parameter :: steps(4) = [-2.0_wp, -1.0_wp, 1.0_wp, 2.0_wp]
      real(wp), parameter :: weights(4) = [1.0_wp, -8.0_wp, 8.0_wp, -1.0_wp]
      type(cavity_type_drop), allocatable :: moved
      type(moist_context_type), target :: ctx
      type(structure_type) :: mol_moved
      real(wp), allocatable :: dA(:, :)
      real(wp) :: fd
      integer :: iat, idir, m

      dA = cavity%A_tot1_rA
      do iat = 1, mol%nat
         do idir = 1, 3
            fd = 0.0_wp
            do m = 1, 4
               mol_moved = mol
               mol_moved%xyz(idir, iat) = mol_moved%xyz(idir, iat) + steps(m)*h
               call build_screened_fixture(moved, ctx, mol_moved, fade_start, error, reach_on, reach_off)
               if (allocated(error)) return
               fd = fd + weights(m)*sum(moved%a(1:moved%ngrid))
            end do
            fd = fd/(12.0_wp*h)
            call check(error, dA(idir, iat), fd, thr=tol, &
                       more="total area gradient, atom "//to_string(iat)//" axis "//to_string(idir))
            if (allocated(error)) return
         end do
      end do
   end subroutine check_area_fd

   !> Rows whose owner, and rows one of whose switching neighbours, is not
   !> LSF-active at the projected point
   !>
   !> @param[in]  cavity      Built cavity
   !> @param[out] n_owner     Rows whose owner is not LSF-active at the point
   !> @param[out] n_neighbour Rows with a switching response to an atom that is not
   !> @param[out] error       Error handle
   !> @param[out] area_nb     Forward total area gradient restricted to those atoms
   !>                         (needs the forward gradient)
   subroutine count_screened(cavity, n_owner, n_neighbour, error, area_nb)
      type(cavity_type_drop), intent(in) :: cavity
      integer, intent(out) :: n_owner, n_neighbour
      type(error_type), allocatable, intent(out) :: error
      real(wp), allocatable, intent(out), optional :: area_nb(:, :)

      class(moist_cavity_drop_lsf_type), allocatable :: lsf
      type(iswig_workspace_type) :: work
      type(mctc_error), allocatable :: lsf_error
      real(wp), allocatable :: grad(:, :)
      logical, allocatable :: active(:)
      integer :: i, ia, k, owner
      logical :: counted

      n_owner = 0
      n_neighbour = 0
      if (present(area_nb)) then
         if (.not. allocated(cavity%a_i1_rA)) then
            call test_failed(error, "screened neighbour area gradient needs the forward gradient")
            return
         end if
         allocate (area_nb(3, cavity%nsph), source=0.0_wp)
      end if
      allocate (lsf, source=cavity%lsf_model)
      call lsf%set_max_deriv(1)
      call cavity%iswig%new_workspace(work)
      allocate (grad(3, cavity%nsph), active(cavity%nsph))
      do i = 1, cavity%ngrid
         call lsf%prepare(cavity%xyz(:, i), lsf_error)
         if (allocated(lsf_error)) then
            call test_failed(error, "level set failed at a grid point: "//lsf_error%message)
            return
         end if
         active = .false.
         do k = 1, lsf%active_count()
            active(lsf%active_atom(k)) = .true.
         end do
         ia = cavity%anchor_of_row(i)
         owner = cavity%owner(i)
         if (.not. active(owner)) n_owner = n_owner + 1
         grad = 0.0_wp
         call cavity%iswig%add_f1_rA(cavity%anchor%xyz(:, ia), owner, cavity%anchor%xi0(ia), 1.0_wp, work, grad)
         counted = .false.
         do k = 1, cavity%nsph
            if (k == owner .or. active(k)) cycle
            if (all(grad(:, k) == 0.0_wp)) cycle
            if (.not. counted) n_neighbour = n_neighbour + 1
            counted = .true.
            if (present(area_nb)) area_nb(:, k) = area_nb(:, k) + cavity%a_i1_rA(:, k, i)
         end do
      end do
   end subroutine count_screened

   !> Cavity of the switch fixtures: 302 points, level 3, SvdW with three-body weight 3
   !>
   !> @param[out]   cavity  Built cavity
   !> @param[inout] ctx     Run context borrowed by the cavity
   !> @param[in]    mol     Structure
   !> @param[in]    blend_k SvdW blending sharpness
   !> @param[in]    dpdist  Branch window (bohr)
   !> @param[out]   error   Error handle
   !> @param[in]    reach_on  Open edge of the reach gate (default: the library's)
   !> @param[in]    reach_off Closed edge of the reach gate (default: the library's)
   subroutine build_switch_fixture(cavity, ctx, mol, blend_k, dpdist, error, reach_on, reach_off)
      type(cavity_type_drop), allocatable, intent(out) :: cavity
      type(moist_context_type), target, intent(inout) :: ctx
      type(structure_type), intent(in) :: mol
      real(wp), intent(in) :: blend_k, dpdist
      type(error_type), allocatable, intent(out) :: error
      real(wp), intent(in), optional :: reach_on, reach_off

      type(moist_cavity_drop_lsf_svdw_type) :: svdw_template
      type(moist_cavity_drop_parameters_type) :: param
      type(mctc_error), allocatable :: cav_error

      call new_context(ctx, nthreads=0, verbosity=0)
      param = moist_cavity_drop_parameters_type(num_leb=302, tolerance=PROJ_TOL, &
                                                proj_level=3, branch=drop_branch_parameters_type(dpdist_max=dpdist), &
                                                do_fine=.true.)
      if (present(reach_on)) param%branch%reach_on = reach_on
      if (present(reach_off)) param%branch%reach_off = reach_off
      call svdw_template%new(param=moist_cavity_drop_lsf_svdw_param_type(blend_k=blend_k, blend_3b=3.0_wp))
      allocate (cavity)
      call new_cavity_drop(cavity, radius_model=default_cpcm_radii(), lsf_model=svdw_template, &
                           error=cav_error, param=param, ctx=ctx)
      if (.not. allocated(cav_error)) call cavity%update(mol, error=cav_error)
      if (allocated(cav_error)) call test_failed(error, "failed to build cavity: "//cav_error%message)
   end subroutine build_switch_fixture

   !> Compare the reverse contraction of a built cavity against its forward arrays
   !>
   !> @param[inout] cavity   Built cavity; its forward arrays are filled here
   !> @param[in]    channels Channel identifiers to populate
   !> @param[in]    label    Human-readable channel description
   !> @param[out]   error    Error handle
   subroutine compare_paths(cavity, channels, label, error)
      !> Built cavity
      type(cavity_type_drop), intent(inout) :: cavity
      !> Channel identifiers to populate
      integer, intent(in) :: channels(:)
      !> Human-readable channel description
      character(len=*), intent(in) :: label
      !> Error handle
      type(error_type), allocatable, intent(out) :: error

      type(cavity_adjoint_type) :: acc
      type(mctc_error), allocatable :: cav_error
      !> Reverse-mode and forward-mode gradients
      real(wp), allocatable :: grad_rev(:, :), grad_fwd(:, :)
      !> Grid extents and loop indices
      integer :: ngrid, nsph, igrid, iatom, iaxis, ich
      !> Comparison bookkeeping
      real(wp) :: diff, scale, worst
      integer :: worst_atom, worst_axis

      ngrid = cavity%ngrid
      nsph = cavity%nsph

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
            ! Negated pass test, so a NaN on either path fails
            if (.not. (diff <= EQ_ABS .or. diff <= EQ_REL*abs(grad_fwd(iaxis, iatom)))) then
               call test_failed(error, "reverse/forward gradient mismatch for "//label// &
                                " at atom "//to_string(iatom)//" axis "//to_string(iaxis)// &
                                ": reverse "//to_string(grad_rev(iaxis, iatom))// &
                                " forward "//to_string(grad_fwd(iaxis, iatom)))
               return
            end if
         end do
      end do

   end subroutine compare_paths

   !> Populate one surface-adjoint channel with a reproducible weight pattern
   !>
   !> The pattern is deterministic and varies across the grid so that a bug
   !> that happens to cancel for uniform weights still shows up
   !>
   !> @param[inout] acc      Surface-adjoint accumulator
   !> @param[in]    channel  Channel identifier
   !> @param[in]    cavity   Cavity supplying the grid extent
   !> @param[out]   error    Error handle
   subroutine apply_channel(acc, channel, cavity, error)
      !> Surface-adjoint accumulator
      type(cavity_adjoint_type), intent(inout) :: acc
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
   !> The two heavy centers give an elongated cavity with well-separated
   !> principal curvatures and the off-axis hydrogen removes the residual
   !> rotational symmetry, so no channel is accidentally degenerate
   !>
   !> The `branching` variant instead uses the five-carbon cross with concave
   !> seams. Branching needs the projector to find sibling minima, which takes
   !> a concave seam and multistart projection -- a plain dimer never branches
   !> at any Lebedev order or branch window -- *and* a window wide enough that
   !> the siblings keep weight. The caller asserts `branch_count > 1`
   !> actually occurred rather than trusting the configuration
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

      !> Displacement axis and atom indices for the symmetry-breaking nudge
      integer :: iat, iaxis
      !> Unit direction of the nudge for one atom
      real(wp) :: dvec(3)

      if (branching) then
         call get_test_cross(mol)
         ! Add tiny perturbation so we find all branches on all systems
         do iat = 1, mol%nat
            do iaxis = 1, 3
               dvec(iaxis) = sin(2.7_wp*real(iat, wp) + 1.1_wp*real(iaxis, wp))
            end do
            mol%xyz(:, iat) = mol%xyz(:, iat) + FIXTURE_NUDGE*dvec/norm2(dvec)
         end do
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
      integer :: proj_level_loc, nleb_loc, prune_loc
      real(wp) :: blend_k_loc, gamma_loc
      !> Fixture settings
      type(moist_cavity_drop_parameters_type) :: param

      nleb_loc = NUM_LEB
      prune_loc = 4
      gamma_loc = gamma
      blend_k_loc = k
      proj_level_loc = PROJ_LEVEL
      if (branching) then
         proj_level_loc = 3
         blend_k_loc = 1.0_wp
         gamma_loc = 1.0_wp
         nleb_loc = 110
         prune_loc = 0
      end if

      call fill_legacy_radii(mol, radii, error)
      if (allocated(error)) return

      ! do_fine turns on curvature, normals, r_iI and pdist so that every
      ! optional *_rA array the reference contraction reads is allocated
      allocate (cavity)
      block
         type(moist_cavity_drop_lsf_svdw_type) :: svdw_template
         call svdw_template%new(param=moist_cavity_drop_lsf_svdw_param_type(blend_k=blend_k_loc, &
                                                                            blend_3b=gamma_loc))
         call new_context(ctx, nthreads=0, verbosity=0)
         param = moist_cavity_drop_parameters_type(num_leb=nleb_loc, tolerance=PROJ_TOL, &
                                                   proj_level=proj_level_loc, &
                                                   wleb_prune_level=prune_loc, &
                                                   do_fine=.true.)
         ! The cross's siblings are backward projections (normal into the anchor
         ! sphere, c < 0), which the signed horizon switch closes under defaults
         if (branching) param%w_0hor_on = 0.0_wp
         call new_cavity_drop(cavity, radius_model=default_cpcm_radii(), lsf_model=svdw_template, &
                              error=cav_error, param=param, ctx=ctx)
      end block
      if (allocated(cav_error)) then
         call test_failed(error, "failed to initialize cavity: "//cav_error%message)
         return
      end if

      ! The projection finds sibling branches, but none inside the default
      ! branch window, so every anchor keeps a single row. Widening the window
      ! lets siblings share, which is what puts the branch channel under test
      if (branching) then
         cavity%param%branch%dpdist_max = BRANCH_DPDIST_MAX
         ! The search window follows the branch window, so recompute the
         ! derived parameters rather than leaving a stale admissible set
         call cavity%param%compute_derived(cav_error)
         if (allocated(cav_error)) then
            call test_failed(error, "failed to recompute derived parameters: "//cav_error%message)
            return
         end if
         call cavity%branch_weight%init(BRANCH_DPDIST_MAX, cavity%param%phi_alpha)
      end if

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
