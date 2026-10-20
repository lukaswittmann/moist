!> Complex and response phases of the generated PCM kernels
!>
!> Tested:
!> - Every model on some iSwiG surface of a trimer (SvdW-DROP for SLIC)
!> - Gaussian and point kernels
!> - Three solvers: dense LAPACK, GMRES dense and matrix-free GMRES
!>
!> Tests:
!> - `response` against six-point differences of the v (and vn) gradient along
!>   the first direction, against `tangent` with only the trace directions set,
!>   and symmetric between two directions
!> - Complex variant with eps + 0i and real traces and directions against the
!>   real phases: energy, charges, gradients, weights, tangent, response and
!>   parameter derivatives; the real phases repeated afterwards are unchanged
!> - Complex variant at a lossy eps: complex response against differences of
!>   the complex v gradient, against the complex trace tangent, and symmetric
!>   under the plain transpose (holomorphic, never conjugated)
!> - Dense: complex-step derivative Im E(eps + ih)/Im(eps) against the real
!>   parameter derivative (direct solves only, the imaginary part sits far
!>   below any GMRES tolerance)
!> - GMRES: every complex phase at the lossy eps against the dense solver
!> - Phase order: real chain, complex chain on the same state, real solve
!>   repeated, then the stencils
module test_model_component_pcm_kernel_variants
   use mctc_env, only: wp, fatal_error
   use mctc_io, only: structure_type
   use mctc_env_error, only: moist_error_type => error_type
   use testdrive, only: new_unittest, unittest_type, error_type, check, test_failed
   use test_helpers, only: fd6_scalar, fd6_offsets, check_moist_error, point_charge_potential, &
      & point_charge_field
   use test_model_component_pcm_fixtures, only: pcm_cpcm, pcm_iefpcm, pcm_ssvpe, pcm_dpcm, &
      & pcm_dpcm_original, pcm_lpb, pcm_slic, pcm_nonlocal, pcm_physics, first_seed, fd_cavity, &
      & get_seed, get_solute, mol_triatomic, host_charges
   use moist_model_continuum_component_pcm_solvers, only: pcm_solver, pcm_dense_solver, &
      & pcm_gmres_solver, pcm_block_jacobi, pcm_storage_dense, pcm_storage_matrix_free
   use pcmgen_gauss_conductor, only: conductor_complex_state, conductor_assemble, &
      & conductor_medium, conductor_solve, conductor_weights, conductor_tangent, &
      & conductor_response, conductor_param_weights, conductor_complex_medium, &
      & conductor_complex_solve, conductor_complex_weights, conductor_complex_tangent, &
      & conductor_complex_response, conductor_complex_param_weights
   use pcmgen_point_conductor, only: conductor_point_complex_state, conductor_point_assemble, &
      & conductor_point_medium, conductor_point_solve, conductor_point_weights, &
      & conductor_point_tangent, conductor_point_response, conductor_point_param_weights, &
      & conductor_point_complex_medium, conductor_point_complex_solve, &
      & conductor_point_complex_weights, conductor_point_complex_tangent, &
      & conductor_point_complex_response, conductor_point_complex_param_weights
   use pcmgen_gauss_iefpcm, only: iefpcm_complex_state, iefpcm_assemble, iefpcm_medium, &
      & iefpcm_solve, iefpcm_weights, iefpcm_tangent, iefpcm_response, iefpcm_param_weights, &
      & iefpcm_complex_medium, iefpcm_complex_solve, iefpcm_complex_weights, &
      & iefpcm_complex_tangent, iefpcm_complex_response, iefpcm_complex_param_weights
   use pcmgen_point_iefpcm, only: iefpcm_point_complex_state, iefpcm_point_assemble, &
      & iefpcm_point_medium, iefpcm_point_solve, iefpcm_point_weights, iefpcm_point_tangent, &
      & iefpcm_point_response, iefpcm_point_param_weights, iefpcm_point_complex_medium, &
      & iefpcm_point_complex_solve, iefpcm_point_complex_weights, iefpcm_point_complex_tangent, &
      & iefpcm_point_complex_response, iefpcm_point_complex_param_weights
   use pcmgen_gauss_ssvpe, only: ssvpe_complex_state, ssvpe_assemble, ssvpe_medium, ssvpe_solve, &
      & ssvpe_weights, ssvpe_tangent, ssvpe_response, ssvpe_param_weights, ssvpe_complex_medium, &
      & ssvpe_complex_solve, ssvpe_complex_weights, ssvpe_complex_tangent, ssvpe_complex_response, &
      & ssvpe_complex_param_weights
   use pcmgen_point_ssvpe, only: ssvpe_point_complex_state, ssvpe_point_assemble, &
      & ssvpe_point_medium, ssvpe_point_solve, ssvpe_point_weights, ssvpe_point_tangent, &
      & ssvpe_point_response, ssvpe_point_param_weights, ssvpe_point_complex_medium, &
      & ssvpe_point_complex_solve, ssvpe_point_complex_weights, ssvpe_point_complex_tangent, &
      & ssvpe_point_complex_response, ssvpe_point_complex_param_weights
   use pcmgen_gauss_dpcm, only: dpcm_complex_state, dpcm_assemble, dpcm_medium, dpcm_solve, &
      & dpcm_weights, dpcm_tangent, dpcm_response, dpcm_param_weights, dpcm_complex_medium, &
      & dpcm_complex_solve, dpcm_complex_weights, dpcm_complex_tangent, dpcm_complex_response, &
      & dpcm_complex_param_weights
   use pcmgen_point_dpcm, only: dpcm_point_complex_state, dpcm_point_assemble, dpcm_point_medium, &
      & dpcm_point_solve, dpcm_point_weights, dpcm_point_tangent, dpcm_point_response, &
      & dpcm_point_param_weights, dpcm_point_complex_medium, dpcm_point_complex_solve, &
      & dpcm_point_complex_weights, dpcm_point_complex_tangent, dpcm_point_complex_response, &
      & dpcm_point_complex_param_weights
   use pcmgen_gauss_dpcm_original, only: dpcm_original_complex_state, dpcm_original_assemble, &
      & dpcm_original_medium, dpcm_original_solve, dpcm_original_weights, dpcm_original_tangent, &
      & dpcm_original_response, dpcm_original_param_weights, dpcm_original_complex_medium, &
      & dpcm_original_complex_solve, dpcm_original_complex_weights, dpcm_original_complex_tangent, &
      & dpcm_original_complex_response, dpcm_original_complex_param_weights
   use pcmgen_point_dpcm_original, only: dpcm_original_point_complex_state, &
      & dpcm_original_point_assemble, dpcm_original_point_medium, dpcm_original_point_solve, &
      & dpcm_original_point_weights, dpcm_original_point_tangent, dpcm_original_point_response, &
      & dpcm_original_point_param_weights, dpcm_original_point_complex_medium, &
      & dpcm_original_point_complex_solve, dpcm_original_point_complex_weights, &
      & dpcm_original_point_complex_tangent, dpcm_original_point_complex_response, &
      & dpcm_original_point_complex_param_weights
   use pcmgen_gauss_lpb, only: lpb_complex_state, lpb_assemble, lpb_medium, lpb_solve, &
      & lpb_weights, lpb_tangent, lpb_response, lpb_param_weights, lpb_complex_medium, &
      & lpb_complex_solve, lpb_complex_weights, lpb_complex_tangent, lpb_complex_response, &
      & lpb_complex_param_weights
   use pcmgen_point_lpb, only: lpb_point_complex_state, lpb_point_assemble, lpb_point_medium, &
      & lpb_point_solve, lpb_point_weights, lpb_point_tangent, lpb_point_response, &
      & lpb_point_param_weights, lpb_point_complex_medium, lpb_point_complex_solve, &
      & lpb_point_complex_weights, lpb_point_complex_tangent, lpb_point_complex_response, &
      & lpb_point_complex_param_weights
   use pcmgen_gauss_slic, only: slic_state, slic_assemble, slic_medium, slic_solve, slic_weights, &
      & slic_tangent, slic_response, slic_param_weights
   use pcmgen_point_slic, only: slic_point_state, slic_point_assemble, slic_point_medium, &
      & slic_point_solve, slic_point_weights, slic_point_tangent, slic_point_response, &
      & slic_point_param_weights
   use pcmgen_gauss_nonlocal, only: nonlocal_pcm_state, nonlocal_pcm_assemble, &
      & nonlocal_pcm_medium, nonlocal_pcm_solve, nonlocal_pcm_weights, nonlocal_pcm_tangent, &
      & nonlocal_pcm_response, nonlocal_pcm_param_weights
   use pcmgen_point_nonlocal, only: nonlocal_pcm_point_state, nonlocal_pcm_point_assemble, &
      & nonlocal_pcm_point_medium, nonlocal_pcm_point_solve, nonlocal_pcm_point_weights, &
      & nonlocal_pcm_point_tangent, nonlocal_pcm_point_response, nonlocal_pcm_point_param_weights
   implicit none(type, external)
   private

   public :: collect_model_component_pcm_kernel_variants

   !> Solvers: dense LAPACK, GMRES on stored matrices, GMRES on matrix-free kernels
   integer, parameter :: solver_dense = 1, solver_gmres = 2, solver_mf = 3
   !> Solver names for messages
   character(len=5), parameter :: solver_name(3) = ["dense", "gmres", "mf   "]
   !> Model names for messages, indexed by the fixture model identifiers
   character(len=13), parameter :: model_name(9) = [character(len=13) :: "conductor", "cosmo", &
      & "iefpcm", "ssvpe", "dpcm", "dpcm_original", "lpb", "slic", "nonlocal"]

   !> Number of direction blocks
   integer, parameter :: nblk = 2
   !> Relative GMRES tolerance
   real(wp), parameter :: gmres_tol = 1.0e-13_wp
   !> Points per block-Jacobi cluster, well below the surface size
   integer, parameter :: jacobi_block = 8
   !> Newton tolerance of SLIC, far below the difference tolerances
   real(wp), parameter :: newton_tol = 1.0e-13_wp
   !> Difference step along the trace directions, linear models
   !>
   !> Trace gradient linear in v and vn: stencil exact at any step
   real(wp), parameter :: step_linear = 0.25_wp
   !> Difference step along the trace directions, SLIC
   !>
   !> Balances the truncation against the Newton and GMRES noise, both at 1e-13
   real(wp), parameter :: step_slic = 1.0e-2_wp
   !> Imaginary part of the complex-step permittivity
   real(wp), parameter :: step_complex = 1.0e-20_wp
   !> Imaginary part of the lossy permittivity
   real(wp), parameter :: eps_loss = 9.0_wp

   !* --------------- Tolerances relative to the largest reference entry -------------- *!

   !> Response against the six-point differences, linear models, real and lossy
   real(wp), parameter :: fd_rtol = 5.0e-12_wp
   !> Response against the differences, SLIC with dense solves
   real(wp), parameter :: fd_slic_rtol = 2.5e-11_wp
   !> Response against the differences, SLIC with GMRES
   real(wp), parameter :: fd_slic_gmres_rtol = 1.0e-10_wp
   !> Response against the trace tangent and between two directions
   real(wp), parameter :: resp_rtol = 5.0e-14_wp
   !> Complex variant at eps + 0i against the real phases
   real(wp), parameter :: cx_rtol = 3.0e-14_wp
   !> Complex-step parameter derivative against the real one
   real(wp), parameter :: cs_rtol = 5.0e-15_wp
   !> GMRES against dense at the lossy eps
   real(wp), parameter :: gmres_rtol = 2.5e-12_wp

   !> Inputs of a kernel run on one surface
   type :: kernel_inputs
      !> Gaussian surface (widths and switching factors set)
      logical :: gaussian = .false.
      !> Number of surface points
      integer :: npt = 0
      !> Physical parameters
      type(pcm_physics) :: phys
      !> Difference step along the trace direction
      real(wp) :: h = 0.0_wp
      !> Positions (3, npt), bohr
      real(wp), allocatable :: r(:, :)
      !> Unit normals (3, npt)
      real(wp), allocatable :: n(:, :)
      !> Principal curvatures (npt), 1/bohr
      real(wp), allocatable :: k1(:), k2(:)
      !> Gaussian widths (npt), Gaussian surface only
      real(wp), allocatable :: xi(:)
      !> Switching factors (npt), Gaussian surface only
      real(wp), allocatable :: f(:)
      !> Areas (npt), bohr^2
      real(wp), allocatable :: a(:)
      !> Potential and normal field traces (npt)
      real(wp), allocatable :: v(:), vn(:)
      !> Position and normal directions (3, npt, nblk)
      real(wp), allocatable :: dr(:, :, :), dn(:, :, :)
      !> Curvature, width, switching and area directions (npt, nblk)
      real(wp), allocatable :: dk1(:, :), dk2(:, :), dxi(:, :), df(:, :), da(:, :)
      !> Trace directions (npt, nblk)
      real(wp), allocatable :: dv(:, :), dvn(:, :)
      !> Zero point-field and vector directions, (3, npt, nblk) and (npt, nblk)
      real(wp), allocatable :: zr(:, :, :), zs(:, :)
   end type kernel_inputs

   !> Real outputs of the phases
   type :: real_work
      !> Charges and trace gradients (npt)
      real(wp), allocatable :: q(:), gv(:), gvn(:)
      !> Interior potential and exterior auxiliary traces of the nonlocal model (npt)
      real(wp), allocatable :: u(:), wext(:)
      !> Weights of positions and normals (3, npt)
      real(wp), allocatable :: gr(:, :), gn(:, :)
      !> Weights of curvatures, widths, switching and areas (npt)
      real(wp), allocatable :: gk1(:), gk2(:), gxi(:), gf(:), ga(:)
      !> Tangents of positions and normals (3, npt, nblk)
      real(wp), allocatable :: hr(:, :, :), hn(:, :, :)
      !> Tangents of the other channels and traces (npt, nblk)
      real(wp), allocatable :: hk1(:, :), hk2(:, :), hxi(:, :), hf(:, :), ha(:, :), hv(:, :), &
         & hvn(:, :)
   end type real_work

   !> Complex outputs of the phases
   type :: complex_work
      !> Charges and trace gradients (npt)
      complex(wp), allocatable :: q(:), gv(:), gvn(:)
      !> Weights of positions and normals (3, npt)
      complex(wp), allocatable :: gr(:, :), gn(:, :)
      !> Weights of curvatures, widths, switching and areas (npt)
      complex(wp), allocatable :: gk1(:), gk2(:), gxi(:), gf(:), ga(:)
      !> Tangents of positions and normals (3, npt, nblk)
      complex(wp), allocatable :: hr(:, :, :), hn(:, :, :)
      !> Tangents of the other channels and traces (npt, nblk)
      complex(wp), allocatable :: hk1(:, :), hk2(:, :), hxi(:, :), hf(:, :), ha(:, :), hv(:, :), &
         & hvn(:, :)
   end type complex_work

   !> Results of a kernel run; trace quantities stack v over vn
   type :: kernel_run
      !> The kernel has a complex variant
      logical :: has_complex = .false.
      !> Energy
      real(wp) :: energy = 0.0_wp
      !> Charges (npt)
      real(wp), allocatable :: q(:)
      !> Trace gradient (2 npt)
      real(wp), allocatable :: grad(:)
      !> Weights, all channels flattened
      real(wp), allocatable :: weights(:)
      !> Tangent along the full directions, all channels flattened
      real(wp), allocatable :: tangent(:)
      !> Response (2 npt, nblk)
      real(wp), allocatable :: response(:, :)
      !> Trace block of the tangent along the trace directions only (2 npt, nblk)
      real(wp), allocatable :: response_tangent(:, :)
      !> Parameter derivatives, first the permittivity (or feps)
      real(wp), allocatable :: params(:)
      !> Trace gradient on the stencil along the first trace direction (2 npt, 6)
      real(wp), allocatable :: stencil(:, :)
      !> Energy of the real solve repeated last
      real(wp) :: energy_again = 0.0_wp
      !> Charges and trace gradient of the real solve repeated last
      real(wp), allocatable :: q_again(:), grad_again(:)
      !> Permittivity (or feps) of the complex phases
      complex(wp) :: cx_param = (0.0_wp, 0.0_wp)
      !> Complex energy
      complex(wp) :: cx_energy = (0.0_wp, 0.0_wp)
      !> Complex counterparts of the real results
      complex(wp), allocatable :: cx_q(:), cx_grad(:), cx_weights(:), cx_tangent(:), &
         & cx_response(:, :), cx_response_tangent(:, :), cx_params(:), cx_stencil(:, :)
   end type kernel_run

contains

   !> Collect all exported unit tests
   subroutine collect_model_component_pcm_kernel_variants(testsuite)
      !> Collection of tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
         & new_unittest("conductor_gauss_dense", test_conductor_gauss_dense), &
         & new_unittest("conductor_gauss_gmres", test_conductor_gauss_gmres), &
         & new_unittest("conductor_gauss_mf", test_conductor_gauss_mf), &
         & new_unittest("conductor_point_dense", test_conductor_point_dense), &
         & new_unittest("conductor_point_gmres", test_conductor_point_gmres), &
         & new_unittest("conductor_point_mf", test_conductor_point_mf), &
         & new_unittest("iefpcm_gauss_dense", test_iefpcm_gauss_dense), &
         & new_unittest("iefpcm_gauss_gmres", test_iefpcm_gauss_gmres), &
         & new_unittest("iefpcm_gauss_mf", test_iefpcm_gauss_mf), &
         & new_unittest("iefpcm_point_dense", test_iefpcm_point_dense), &
         & new_unittest("iefpcm_point_gmres", test_iefpcm_point_gmres), &
         & new_unittest("iefpcm_point_mf", test_iefpcm_point_mf), &
         & new_unittest("ssvpe_gauss_dense", test_ssvpe_gauss_dense), &
         & new_unittest("ssvpe_gauss_gmres", test_ssvpe_gauss_gmres), &
         & new_unittest("ssvpe_gauss_mf", test_ssvpe_gauss_mf), &
         & new_unittest("ssvpe_point_dense", test_ssvpe_point_dense), &
         & new_unittest("ssvpe_point_gmres", test_ssvpe_point_gmres), &
         & new_unittest("ssvpe_point_mf", test_ssvpe_point_mf), &
         & new_unittest("dpcm_gauss_dense", test_dpcm_gauss_dense), &
         & new_unittest("dpcm_gauss_gmres", test_dpcm_gauss_gmres), &
         & new_unittest("dpcm_gauss_mf", test_dpcm_gauss_mf), &
         & new_unittest("dpcm_point_dense", test_dpcm_point_dense), &
         & new_unittest("dpcm_point_gmres", test_dpcm_point_gmres), &
         & new_unittest("dpcm_point_mf", test_dpcm_point_mf), &
         & new_unittest("dpcm_original_gauss_dense", test_dpcm_original_gauss_dense), &
         & new_unittest("dpcm_original_gauss_gmres", test_dpcm_original_gauss_gmres), &
         & new_unittest("dpcm_original_gauss_mf", test_dpcm_original_gauss_mf), &
         & new_unittest("dpcm_original_point_dense", test_dpcm_original_point_dense), &
         & new_unittest("dpcm_original_point_gmres", test_dpcm_original_point_gmres), &
         & new_unittest("dpcm_original_point_mf", test_dpcm_original_point_mf), &
         & new_unittest("lpb_gauss_dense", test_lpb_gauss_dense), &
         & new_unittest("lpb_gauss_gmres", test_lpb_gauss_gmres), &
         & new_unittest("lpb_gauss_mf", test_lpb_gauss_mf), &
         & new_unittest("lpb_point_dense", test_lpb_point_dense), &
         & new_unittest("lpb_point_gmres", test_lpb_point_gmres), &
         & new_unittest("lpb_point_mf", test_lpb_point_mf), &
         & new_unittest("slic_gauss_dense", test_slic_gauss_dense), &
         & new_unittest("slic_gauss_gmres", test_slic_gauss_gmres), &
         & new_unittest("slic_gauss_mf", test_slic_gauss_mf), &
         & new_unittest("slic_point_dense", test_slic_point_dense), &
         & new_unittest("slic_point_gmres", test_slic_point_gmres), &
         & new_unittest("slic_point_mf", test_slic_point_mf), &
         & new_unittest("nonlocal_gauss_dense", test_nonlocal_gauss_dense), &
         & new_unittest("nonlocal_gauss_gmres", test_nonlocal_gauss_gmres), &
         & new_unittest("nonlocal_gauss_mf", test_nonlocal_gauss_mf), &
         & new_unittest("nonlocal_point_dense", test_nonlocal_point_dense), &
         & new_unittest("nonlocal_point_gmres", test_nonlocal_point_gmres), &
         & new_unittest("nonlocal_point_mf", test_nonlocal_point_mf) &
         & ]
   end subroutine collect_model_component_pcm_kernel_variants

   !> Conductor, Gaussian surface, dense solver
   subroutine test_conductor_gauss_dense(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call check_variants(error, pcm_cpcm, .true., solver_dense)
   end subroutine test_conductor_gauss_dense

   !> Conductor, Gaussian surface, GMRES on stored matrices
   subroutine test_conductor_gauss_gmres(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call check_variants(error, pcm_cpcm, .true., solver_gmres)
   end subroutine test_conductor_gauss_gmres

   !> Conductor, Gaussian surface, matrix-free GMRES
   subroutine test_conductor_gauss_mf(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call check_variants(error, pcm_cpcm, .true., solver_mf)
   end subroutine test_conductor_gauss_mf

   !> Conductor, point surface, dense solver
   subroutine test_conductor_point_dense(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call check_variants(error, pcm_cpcm, .false., solver_dense)
   end subroutine test_conductor_point_dense

   !> Conductor, point surface, GMRES on stored matrices
   subroutine test_conductor_point_gmres(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call check_variants(error, pcm_cpcm, .false., solver_gmres)
   end subroutine test_conductor_point_gmres

   !> Conductor, point surface, matrix-free GMRES
   subroutine test_conductor_point_mf(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call check_variants(error, pcm_cpcm, .false., solver_mf)
   end subroutine test_conductor_point_mf

   !> IEF-PCM, Gaussian surface, dense solver
   subroutine test_iefpcm_gauss_dense(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call check_variants(error, pcm_iefpcm, .true., solver_dense)
   end subroutine test_iefpcm_gauss_dense

   !> IEF-PCM, Gaussian surface, GMRES on stored matrices
   subroutine test_iefpcm_gauss_gmres(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call check_variants(error, pcm_iefpcm, .true., solver_gmres)
   end subroutine test_iefpcm_gauss_gmres

   !> IEF-PCM, Gaussian surface, matrix-free GMRES
   subroutine test_iefpcm_gauss_mf(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call check_variants(error, pcm_iefpcm, .true., solver_mf)
   end subroutine test_iefpcm_gauss_mf

   !> IEF-PCM, point surface, dense solver
   subroutine test_iefpcm_point_dense(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call check_variants(error, pcm_iefpcm, .false., solver_dense)
   end subroutine test_iefpcm_point_dense

   !> IEF-PCM, point surface, GMRES on stored matrices
   subroutine test_iefpcm_point_gmres(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call check_variants(error, pcm_iefpcm, .false., solver_gmres)
   end subroutine test_iefpcm_point_gmres

   !> IEF-PCM, point surface, matrix-free GMRES
   subroutine test_iefpcm_point_mf(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call check_variants(error, pcm_iefpcm, .false., solver_mf)
   end subroutine test_iefpcm_point_mf

   !> SS(V)PE, Gaussian surface, dense solver
   subroutine test_ssvpe_gauss_dense(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call check_variants(error, pcm_ssvpe, .true., solver_dense)
   end subroutine test_ssvpe_gauss_dense

   !> SS(V)PE, Gaussian surface, GMRES on stored matrices
   subroutine test_ssvpe_gauss_gmres(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call check_variants(error, pcm_ssvpe, .true., solver_gmres)
   end subroutine test_ssvpe_gauss_gmres

   !> SS(V)PE, Gaussian surface, matrix-free GMRES
   subroutine test_ssvpe_gauss_mf(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call check_variants(error, pcm_ssvpe, .true., solver_mf)
   end subroutine test_ssvpe_gauss_mf

   !> SS(V)PE, point surface, dense solver
   subroutine test_ssvpe_point_dense(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call check_variants(error, pcm_ssvpe, .false., solver_dense)
   end subroutine test_ssvpe_point_dense

   !> SS(V)PE, point surface, GMRES on stored matrices
   subroutine test_ssvpe_point_gmres(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call check_variants(error, pcm_ssvpe, .false., solver_gmres)
   end subroutine test_ssvpe_point_gmres

   !> SS(V)PE, point surface, matrix-free GMRES
   subroutine test_ssvpe_point_mf(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call check_variants(error, pcm_ssvpe, .false., solver_mf)
   end subroutine test_ssvpe_point_mf

   !> DPCM, Gaussian surface, dense solver
   subroutine test_dpcm_gauss_dense(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call check_variants(error, pcm_dpcm, .true., solver_dense)
   end subroutine test_dpcm_gauss_dense

   !> DPCM, Gaussian surface, GMRES on stored matrices
   subroutine test_dpcm_gauss_gmres(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call check_variants(error, pcm_dpcm, .true., solver_gmres)
   end subroutine test_dpcm_gauss_gmres

   !> DPCM, Gaussian surface, matrix-free GMRES
   subroutine test_dpcm_gauss_mf(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call check_variants(error, pcm_dpcm, .true., solver_mf)
   end subroutine test_dpcm_gauss_mf

   !> DPCM, point surface, dense solver
   subroutine test_dpcm_point_dense(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call check_variants(error, pcm_dpcm, .false., solver_dense)
   end subroutine test_dpcm_point_dense

   !> DPCM, point surface, GMRES on stored matrices
   subroutine test_dpcm_point_gmres(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call check_variants(error, pcm_dpcm, .false., solver_gmres)
   end subroutine test_dpcm_point_gmres

   !> DPCM, point surface, matrix-free GMRES
   subroutine test_dpcm_point_mf(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call check_variants(error, pcm_dpcm, .false., solver_mf)
   end subroutine test_dpcm_point_mf

   !> Original DPCM, Gaussian surface, dense solver
   subroutine test_dpcm_original_gauss_dense(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call check_variants(error, pcm_dpcm_original, .true., solver_dense)
   end subroutine test_dpcm_original_gauss_dense

   !> Original DPCM, Gaussian surface, GMRES on stored matrices
   subroutine test_dpcm_original_gauss_gmres(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call check_variants(error, pcm_dpcm_original, .true., solver_gmres)
   end subroutine test_dpcm_original_gauss_gmres

   !> Original DPCM, Gaussian surface, matrix-free GMRES
   subroutine test_dpcm_original_gauss_mf(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call check_variants(error, pcm_dpcm_original, .true., solver_mf)
   end subroutine test_dpcm_original_gauss_mf

   !> Original DPCM, point surface, dense solver
   subroutine test_dpcm_original_point_dense(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call check_variants(error, pcm_dpcm_original, .false., solver_dense)
   end subroutine test_dpcm_original_point_dense

   !> Original DPCM, point surface, GMRES on stored matrices
   subroutine test_dpcm_original_point_gmres(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call check_variants(error, pcm_dpcm_original, .false., solver_gmres)
   end subroutine test_dpcm_original_point_gmres

   !> Original DPCM, point surface, matrix-free GMRES
   subroutine test_dpcm_original_point_mf(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call check_variants(error, pcm_dpcm_original, .false., solver_mf)
   end subroutine test_dpcm_original_point_mf

   !> LPB, Gaussian surface, dense solver
   subroutine test_lpb_gauss_dense(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call check_variants(error, pcm_lpb, .true., solver_dense)
   end subroutine test_lpb_gauss_dense

   !> LPB, Gaussian surface, GMRES on stored matrices
   subroutine test_lpb_gauss_gmres(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call check_variants(error, pcm_lpb, .true., solver_gmres)
   end subroutine test_lpb_gauss_gmres

   !> LPB, Gaussian surface, matrix-free GMRES
   subroutine test_lpb_gauss_mf(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call check_variants(error, pcm_lpb, .true., solver_mf)
   end subroutine test_lpb_gauss_mf

   !> LPB, point surface, dense solver
   subroutine test_lpb_point_dense(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call check_variants(error, pcm_lpb, .false., solver_dense)
   end subroutine test_lpb_point_dense

   !> LPB, point surface, GMRES on stored matrices
   subroutine test_lpb_point_gmres(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call check_variants(error, pcm_lpb, .false., solver_gmres)
   end subroutine test_lpb_point_gmres

   !> LPB, point surface, matrix-free GMRES
   subroutine test_lpb_point_mf(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call check_variants(error, pcm_lpb, .false., solver_mf)
   end subroutine test_lpb_point_mf

   !> SLIC, Gaussian surface, dense solver
   subroutine test_slic_gauss_dense(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call check_variants(error, pcm_slic, .true., solver_dense)
   end subroutine test_slic_gauss_dense

   !> SLIC, Gaussian surface, GMRES on stored matrices
   subroutine test_slic_gauss_gmres(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call check_variants(error, pcm_slic, .true., solver_gmres)
   end subroutine test_slic_gauss_gmres

   !> SLIC, Gaussian surface, matrix-free GMRES
   subroutine test_slic_gauss_mf(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call check_variants(error, pcm_slic, .true., solver_mf)
   end subroutine test_slic_gauss_mf

   !> SLIC, point surface, dense solver
   subroutine test_slic_point_dense(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call check_variants(error, pcm_slic, .false., solver_dense)
   end subroutine test_slic_point_dense

   !> SLIC, point surface, GMRES on stored matrices
   subroutine test_slic_point_gmres(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call check_variants(error, pcm_slic, .false., solver_gmres)
   end subroutine test_slic_point_gmres

   !> SLIC, point surface, matrix-free GMRES
   subroutine test_slic_point_mf(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call check_variants(error, pcm_slic, .false., solver_mf)
   end subroutine test_slic_point_mf

   !> Nonlocal, Gaussian surface, dense solver
   subroutine test_nonlocal_gauss_dense(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call check_variants(error, pcm_nonlocal, .true., solver_dense)
   end subroutine test_nonlocal_gauss_dense

   !> Nonlocal, Gaussian surface, GMRES on stored matrices
   subroutine test_nonlocal_gauss_gmres(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call check_variants(error, pcm_nonlocal, .true., solver_gmres)
   end subroutine test_nonlocal_gauss_gmres

   !> Nonlocal, Gaussian surface, matrix-free GMRES
   subroutine test_nonlocal_gauss_mf(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call check_variants(error, pcm_nonlocal, .true., solver_mf)
   end subroutine test_nonlocal_gauss_mf

   !> Nonlocal, point surface, dense solver
   subroutine test_nonlocal_point_dense(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call check_variants(error, pcm_nonlocal, .false., solver_dense)
   end subroutine test_nonlocal_point_dense

   !> Nonlocal, point surface, GMRES on stored matrices
   subroutine test_nonlocal_point_gmres(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call check_variants(error, pcm_nonlocal, .false., solver_gmres)
   end subroutine test_nonlocal_point_gmres

   !> Nonlocal, point surface, matrix-free GMRES
   subroutine test_nonlocal_point_mf(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      call check_variants(error, pcm_nonlocal, .false., solver_mf)
   end subroutine test_nonlocal_point_mf

   !> Run every check of one model, surface and solver
   !>
   !> @param[out] error    Test failure
   !> @param[in]  model    Model identifier of the fixtures
   !> @param[in]  gaussian Gaussian surface
   !> @param[in]  kind     Solver kind
   subroutine check_variants(error, model, gaussian, kind)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Model identifier of the fixtures
      integer, intent(in) :: model
      !> Gaussian surface
      logical, intent(in) :: gaussian
      !> Solver kind
      integer, intent(in) :: kind

      !> Library error handling
      type(moist_error_type), allocatable :: err
      !> Kernel inputs
      type(kernel_inputs) :: inp
      !> Solver under test and dense reference solver
      class(pcm_solver), allocatable :: solver, dense
      !> Runs at eps + 0i, complex-step eps, lossy eps, lossy eps with dense solves
      type(kernel_run) :: run, step, lossy, ref
      !> Failure context prefix
      character(len=:), allocatable :: prefix
      !> Complex-step derivative
      real(wp) :: dstep
      !> Tolerance of the differences
      real(wp) :: fd_tol

      prefix = trim(model_name(model))//" "//trim(merge("gauss", "point", gaussian))//" " &
         & //trim(solver_name(kind))
      call new_inputs(model, gaussian, inp, error)
      if (allocated(error)) return
      call new_solver(kind, solver)
      call run_kernel(model, inp, solver, cmplx(inp%phys%eps, 0.0_wp, wp), run, err)
      call check_moist_error(error, err, prefix//": run")
      if (allocated(error)) return

      fd_tol = fd_rtol
      if (model == pcm_slic) fd_tol = merge(fd_slic_rtol, fd_slic_gmres_rtol, kind == solver_dense)
      call check_response(error, prefix, "response", inp, cmplx(run%response, 0.0_wp, wp), &
         & cmplx(run%response_tangent, 0.0_wp, wp), cmplx(run%stencil, 0.0_wp, wp), fd_tol)
      if (allocated(error)) return
      call compare(error, prefix, "repeated energy", [run%energy_again], [run%energy], cx_rtol)
      if (allocated(error)) return
      call compare(error, prefix, "repeated charges", run%q_again, run%q, cx_rtol)
      if (allocated(error)) return
      call compare(error, prefix, "repeated trace gradient", run%grad_again, run%grad, cx_rtol)
      if (allocated(error)) return
      if (.not. run%has_complex) return

      call check_complex(error, prefix, run)
      if (allocated(error)) return
      call run_kernel(model, inp, solver, cmplx(inp%phys%eps, eps_loss, wp), lossy, err)
      call check_moist_error(error, err, prefix//": lossy run")
      if (allocated(error)) return
      call check_response(error, prefix, "lossy response", inp, lossy%cx_response, &
         & lossy%cx_response_tangent, lossy%cx_stencil, fd_rtol)
      if (allocated(error)) return
      if (kind == solver_dense) then
         call run_kernel(model, inp, solver, cmplx(inp%phys%eps, step_complex, wp), step, err)
         call check_moist_error(error, err, prefix//": complex-step run")
         if (allocated(error)) return
         dstep = aimag(step%cx_energy)/aimag(step%cx_param)
         call compare(error, prefix, "complex-step dE/deps", [dstep], run%params(1:1), cs_rtol)
      else
         call new_solver(solver_dense, dense)
         call run_kernel(model, inp, dense, cmplx(inp%phys%eps, eps_loss, wp), ref, err)
         call check_moist_error(error, err, prefix//": lossy dense run")
         if (allocated(error)) return
         call check_lossy(error, prefix, lossy, ref)
      end if
   end subroutine check_variants

   !> Response against differences, the trace tangent and its plain transpose
   !>
   !> Complex responses are holomorphic: the symmetry uses dot products without
   !> conjugation
   !>
   !> @param[out] error    Test failure
   !> @param[in]  prefix   Failure context prefix
   !> @param[in]  name     Quantity name
   !> @param[in]  inp      Kernel inputs
   !> @param[in]  response Response (2 npt, nblk)
   !> @param[in]  trace    Trace block of the tangent along the trace directions (2 npt, nblk)
   !> @param[in]  stencil  Trace gradient on the stencil (2 npt, 6), in `fd6_offsets` order
   !> @param[in]  fd_tol   Tolerance of the differences
   subroutine check_response(error, prefix, name, inp, response, trace, stencil, fd_tol)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Failure context prefix
      character(len=*), intent(in) :: prefix
      !> Quantity name
      character(len=*), intent(in) :: name
      !> Kernel inputs
      type(kernel_inputs), intent(in) :: inp
      !> Response
      complex(wp), intent(in) :: response(:, :)
      !> Trace block of the tangent along the trace directions
      complex(wp), intent(in) :: trace(:, :)
      !> Trace gradient on the stencil
      complex(wp), intent(in) :: stencil(:, :)
      !> Tolerance of the differences
      real(wp), intent(in) :: fd_tol

      !> Differences of the trace gradient
      complex(wp), allocatable :: fd(:)
      !> Differences of the real and imaginary parts
      real(wp) :: re, im
      !> Stacked trace directions
      real(wp) :: d1(2*inp%npt), d2(2*inp%npt)
      !> Entry index
      integer :: i

      allocate (fd(size(stencil, 1)))
      do i = 1, size(fd)
         call fd6_scalar(real(stencil(i, 1)), real(stencil(i, 2)), real(stencil(i, 3)), &
            & real(stencil(i, 4)), real(stencil(i, 5)), real(stencil(i, 6)), inp%h, re, error)
         if (allocated(error)) return
         call fd6_scalar(aimag(stencil(i, 1)), aimag(stencil(i, 2)), aimag(stencil(i, 3)), &
            & aimag(stencil(i, 4)), aimag(stencil(i, 5)), aimag(stencil(i, 6)), inp%h, im, error)
         if (allocated(error)) return
         fd(i) = cmplx(re, im, wp)
      end do
      call compare(error, prefix, name//" against differences", response(:, 1), fd, fd_tol)
      if (allocated(error)) return
      call compare(error, prefix, name//" against tangent", pack(response, .true.), &
         & pack(trace, .true.), resp_rtol)
      if (allocated(error)) return
      d1 = [inp%dv(:, 1), inp%dvn(:, 1)]
      d2 = [inp%dv(:, 2), inp%dvn(:, 2)]
      call compare(error, prefix, name//" symmetry", [sum(d1*response(:, 2))], &
         & [sum(d2*response(:, 1))], resp_rtol)
   end subroutine check_response

   !> Complex variant at eps + 0i against the real phases
   !>
   !> @param[out] error  Test failure
   !> @param[in]  prefix Failure context prefix
   !> @param[in]  run    Kernel results
   subroutine check_complex(error, prefix, run)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Failure context prefix
      character(len=*), intent(in) :: prefix
      !> Kernel results
      type(kernel_run), intent(in) :: run

      call compare(error, prefix, "complex energy", [run%cx_energy], [run%energy], cx_rtol)
      if (allocated(error)) return
      call compare(error, prefix, "complex charges", run%cx_q, run%q, cx_rtol)
      if (allocated(error)) return
      call compare(error, prefix, "complex trace gradient", run%cx_grad, run%grad, cx_rtol)
      if (allocated(error)) return
      call compare(error, prefix, "complex weights", run%cx_weights, run%weights, cx_rtol)
      if (allocated(error)) return
      call compare(error, prefix, "complex tangent", run%cx_tangent, run%tangent, cx_rtol)
      if (allocated(error)) return
      call compare(error, prefix, "complex response", pack(run%cx_response, .true.), &
         & pack(run%response, .true.), cx_rtol)
      if (allocated(error)) return
      call compare(error, prefix, "complex parameter derivatives", run%cx_params, run%params, &
         & cx_rtol)
   end subroutine check_complex

   !> Complex phases at the lossy eps against the dense solver
   !>
   !> @param[out] error  Test failure
   !> @param[in]  prefix Failure context prefix
   !> @param[in]  run    Results of the solver under test
   !> @param[in]  ref    Results of the dense solver
   subroutine check_lossy(error, prefix, run, ref)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Failure context prefix
      character(len=*), intent(in) :: prefix
      !> Results of the solver under test
      type(kernel_run), intent(in) :: run
      !> Results of the dense solver
      type(kernel_run), intent(in) :: ref

      call compare(error, prefix, "lossy energy", [run%cx_energy], [ref%cx_energy], gmres_rtol)
      if (allocated(error)) return
      call compare(error, prefix, "lossy charges", run%cx_q, ref%cx_q, gmres_rtol)
      if (allocated(error)) return
      call compare(error, prefix, "lossy trace gradient", run%cx_grad, ref%cx_grad, gmres_rtol)
      if (allocated(error)) return
      call compare(error, prefix, "lossy weights", run%cx_weights, ref%cx_weights, gmres_rtol)
      if (allocated(error)) return
      call compare(error, prefix, "lossy tangent", run%cx_tangent, ref%cx_tangent, gmres_rtol)
      if (allocated(error)) return
      call compare(error, prefix, "lossy response", pack(run%cx_response, .true.), &
         & pack(ref%cx_response, .true.), gmres_rtol)
      if (allocated(error)) return
      call compare(error, prefix, "lossy parameter derivatives", run%cx_params, ref%cx_params, &
         & gmres_rtol)
   end subroutine check_lossy

   !> Compare two arrays: max |got - ref| <= rtol max |ref|, ref not vanishing
   !>
   !> @param[out] error  Test failure
   !> @param[in]  prefix Failure context prefix
   !> @param[in]  name   Quantity name
   !> @param[in]  got    Values under test
   !> @param[in]  ref    Reference values
   !> @param[in]  rtol   Tolerance relative to the largest reference entry
   subroutine compare(error, prefix, name, got, ref, rtol)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Failure context prefix
      character(len=*), intent(in) :: prefix
      !> Quantity name
      character(len=*), intent(in) :: name
      !> Values under test, real or complex
      class(*), intent(in) :: got(:)
      !> Reference values, real or complex
      class(*), intent(in) :: ref(:)
      !> Tolerance relative to the largest reference entry
      real(wp), intent(in) :: rtol

      !> Values as complex numbers
      complex(wp), allocatable :: x(:), y(:)
      !> Largest deviation and largest reference entry
      real(wp) :: dev, scale
      !> Failure context
      character(len=200) :: context

      call as_complex(got, x)
      call as_complex(ref, y)
      if (size(x) /= size(y)) then
         call test_failed(error, prefix//", "//name//": size mismatch")
         return
      end if
      dev = maxval(abs(x - y))
      scale = maxval(abs(y))
      write (context, '(a,", ",a,": max dev ",es10.2," scale ",es10.2)') prefix, name, dev, scale
      call check(error, scale > 0.0_wp, message="vacuous reference", more=trim(context))
      if (allocated(error)) return
      call check(error, dev <= rtol*scale, message="deviation above tolerance", &
         & more=trim(context))
   end subroutine compare

   !> Copy a real or complex array into a complex one
   !>
   !> @param[in]  values Real or complex values
   !> @param[out] x      Complex copy
   subroutine as_complex(values, x)
      !> Real or complex values
      class(*), intent(in) :: values(:)
      !> Complex copy
      complex(wp), allocatable, intent(out) :: x(:)

      select type (values)
      type is (real(wp))
         x = cmplx(values, 0.0_wp, wp)
      type is (complex(wp))
         x = values
      class default
         allocate (x(0))
      end select
   end subroutine as_complex

   !> Solver of a kind; GMRES with small block-Jacobi clusters
   !>
   !> @param[in]  kind   Solver kind
   !> @param[out] solver Solver
   subroutine new_solver(kind, solver)
      !> Solver kind
      integer, intent(in) :: kind
      !> Solver
      class(pcm_solver), allocatable, intent(out) :: solver

      if (kind == solver_dense) then
         allocate (pcm_dense_solver :: solver)
      else
         allocate (pcm_gmres_solver :: solver)
         select type (solver)
         type is (pcm_gmres_solver)
            solver%tol = gmres_tol
            solver%storage = merge(pcm_storage_matrix_free, pcm_storage_dense, kind == solver_mf)
            allocate (pcm_block_jacobi :: solver%preconditioner)
            select type (preconditioner => solver%preconditioner)
            type is (pcm_block_jacobi)
               preconditioner%block_size = jacobi_block
            end select
         end select
      end if
      solver%nonlinear_tol = newton_tol
   end subroutine new_solver

   !> Kernel inputs on the seed surface of the triatomic
   !>
   !> - Traces of the fixture host charges: potential and normal field
   !> - Smooth deterministic directions of every channel, nblk blocks
   !>
   !> @param[in]  model    Model identifier of the fixtures
   !> @param[in]  gaussian Gaussian surface
   !> @param[out] inp      Kernel inputs
   !> @param[out] error    Test failure
   subroutine new_inputs(model, gaussian, inp, error)
      !> Model identifier of the fixtures
      integer, intent(in) :: model
      !> Gaussian surface
      logical, intent(in) :: gaussian
      !> Kernel inputs
      type(kernel_inputs), intent(out) :: inp
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      !> Solute
      type(structure_type) :: mol
      !> Synthetic cavity
      type(fd_cavity) :: cavity
      !> Host charges
      real(wp), allocatable :: qsrc(:)
      !> Potential derivatives of the host charges
      real(wp), allocatable :: g(:, :), gxi(:), hrr(:, :, :), hrxi(:, :)
      !> Trace scales
      real(wp) :: vscale, vnscale
      !> Point, block and component indices
      integer :: i, ib, ic
      !> Phase of the direction patterns
      real(wp) :: t

      call get_solute(mol_triatomic, mol)
      call get_seed(first_seed(model), gaussian, mol, cavity, error)
      if (allocated(error)) return
      inp%gaussian = gaussian
      inp%npt = cavity%ngrid
      inp%h = merge(step_slic, step_linear, model == pcm_slic)
      associate (npt => inp%npt)
         inp%r = cavity%xyz(:, :npt)
         inp%n = cavity%normal0(:, :npt)
         inp%k1 = cavity%k1(:npt)
         inp%k2 = cavity%k2(:npt)
         inp%a = cavity%a(:npt)
         if (gaussian) then
            inp%xi = cavity%xi0(:npt)
            inp%f = cavity%f(:npt)
         end if
         qsrc = host_charges(mol)
         inp%v = point_charge_potential(cavity, qsrc, mol, gaussian)
         call point_charge_field(cavity, qsrc, mol, gaussian, g, gxi, hrr, hrxi)
         allocate (inp%vn(npt))
         do i = 1, npt
            inp%vn(i) = dot_product(inp%n(:, i), g(:, i))
         end do
         vscale = maxval(abs(inp%v))
         vnscale = maxval(abs(inp%vn))

         allocate (inp%dr(3, npt, nblk), inp%dn(3, npt, nblk), inp%dk1(npt, nblk), &
            & inp%dk2(npt, nblk), inp%dxi(npt, nblk), inp%df(npt, nblk), inp%da(npt, nblk), &
            & inp%dv(npt, nblk), inp%dvn(npt, nblk))
         do ib = 1, nblk
            do i = 1, npt
               t = 0.37_wp*real(i, wp) + 1.1_wp*real(ib, wp)
               do ic = 1, 3
                  inp%dr(ic, i, ib) = 0.05_wp*sin(t + 0.9_wp*real(ic, wp))
                  inp%dn(ic, i, ib) = 0.05_wp*cos(1.3_wp*t + 0.7_wp*real(ic, wp))
               end do
               inp%dk1(i, ib) = 0.02_wp*sin(1.7_wp*t)
               inp%dk2(i, ib) = 0.02_wp*cos(0.6_wp*t)
               inp%dxi(i, ib) = 0.05_wp*sin(2.3_wp*t)
               inp%df(i, ib) = 0.05_wp*cos(1.9_wp*t)
               inp%da(i, ib) = 0.01_wp*sin(0.8_wp*t)
               inp%dv(i, ib) = 0.5_wp*vscale*sin(2.9_wp*t)
               inp%dvn(i, ib) = 0.5_wp*vnscale*cos(3.1_wp*t)
            end do
         end do
         allocate (inp%zr(3, npt, nblk), inp%zs(npt, nblk), source=0.0_wp)
      end associate
   end subroutine new_inputs

   !> Real output buffers of a surface, zero
   !>
   !> @param[in]  npt Number of surface points
   !> @param[out] w   Buffers
   subroutine new_real_work(npt, w)
      !> Number of surface points
      integer, intent(in) :: npt
      !> Buffers
      type(real_work), intent(out) :: w

      allocate (w%q(npt), w%gv(npt), w%gvn(npt), w%u(npt), w%wext(npt), w%gr(3, npt), &
         & w%gn(3, npt), w%gk1(npt), w%gk2(npt), w%gxi(npt), w%gf(npt), w%ga(npt), &
         & w%hr(3, npt, nblk), w%hn(3, npt, nblk), w%hk1(npt, nblk), w%hk2(npt, nblk), &
         & w%hxi(npt, nblk), w%hf(npt, nblk), w%ha(npt, nblk), w%hv(npt, nblk), &
         & w%hvn(npt, nblk), source=0.0_wp)
   end subroutine new_real_work

   !> Complex output buffers of a surface, zero
   !>
   !> @param[in]  npt Number of surface points
   !> @param[out] w   Buffers
   subroutine new_complex_work(npt, w)
      !> Number of surface points
      integer, intent(in) :: npt
      !> Buffers
      type(complex_work), intent(out) :: w

      allocate (w%q(npt), w%gv(npt), w%gvn(npt), w%gr(3, npt), w%gn(3, npt), w%gk1(npt), &
         & w%gk2(npt), w%gxi(npt), w%gf(npt), w%ga(npt), w%hr(3, npt, nblk), &
         & w%hn(3, npt, nblk), w%hk1(npt, nblk), w%hk2(npt, nblk), w%hxi(npt, nblk), &
         & w%hf(npt, nblk), w%ha(npt, nblk), w%hv(npt, nblk), w%hvn(npt, nblk), &
         & source=(0.0_wp, 0.0_wp))
   end subroutine new_complex_work

   !> Keep the real solve: energy, charges and stacked trace gradient
   !>
   !> @param[in,out] run    Kernel results
   !> @param[in]     energy Energy
   !> @param[in]     w      Buffers holding the solve outputs
   subroutine keep_solve(run, energy, w)
      !> Kernel results
      type(kernel_run), intent(inout) :: run
      !> Energy
      real(wp), intent(in) :: energy
      !> Buffers holding the solve outputs
      type(real_work), intent(in) :: w

      run%energy = energy
      run%q = w%q
      run%grad = [w%gv, w%gvn]
   end subroutine keep_solve

   !> Keep the real weights, all channels flattened
   !>
   !> @param[in,out] run Kernel results
   !> @param[in]     w   Buffers holding the weights
   subroutine keep_weights(run, w)
      !> Kernel results
      type(kernel_run), intent(inout) :: run
      !> Buffers holding the weights
      type(real_work), intent(in) :: w

      run%weights = [pack(w%gr, .true.), pack(w%gn, .true.), w%gk1, w%gk2, w%gxi, w%gf, w%ga]
   end subroutine keep_weights

   !> Keep the real tangent, all channels flattened
   !>
   !> @param[in,out] run Kernel results
   !> @param[in]     w   Buffers holding the tangent
   subroutine keep_tangent(run, w)
      !> Kernel results
      type(kernel_run), intent(inout) :: run
      !> Buffers holding the tangent
      type(real_work), intent(in) :: w

      run%tangent = [pack(w%hr, .true.), pack(w%hn, .true.), pack(w%hk1, .true.), &
         & pack(w%hk2, .true.), pack(w%hxi, .true.), pack(w%hf, .true.), pack(w%ha, .true.), &
         & pack(w%hv, .true.), pack(w%hvn, .true.)]
   end subroutine keep_tangent

   !> Trace blocks stacked, v over vn
   !>
   !> @param[in] hv  v block (npt, nblk)
   !> @param[in] hvn vn block (npt, nblk)
   pure function stack(hv, hvn) result(h)
      !> v block
      real(wp), intent(in) :: hv(:, :)
      !> vn block
      real(wp), intent(in) :: hvn(:, :)
      !> Stacked blocks (2 npt, nblk)
      real(wp) :: h(2*size(hv, 1), size(hv, 2))

      h(:size(hv, 1), :) = hv
      h(size(hv, 1) + 1:, :) = hvn
   end function stack

   !> Keep the real solve repeated after the complex phases
   !>
   !> @param[in,out] run    Kernel results
   !> @param[in]     energy Energy
   !> @param[in]     w      Buffers holding the solve outputs
   subroutine keep_again(run, energy, w)
      !> Kernel results
      type(kernel_run), intent(inout) :: run
      !> Energy
      real(wp), intent(in) :: energy
      !> Buffers holding the solve outputs
      type(real_work), intent(in) :: w

      run%energy_again = energy
      run%q_again = w%q
      run%grad_again = [w%gv, w%gvn]
   end subroutine keep_again

   !> Keep the complex solve: energy, charges and stacked trace gradient
   !>
   !> @param[in,out] run    Kernel results
   !> @param[in]     energy Energy
   !> @param[in]     z      Buffers holding the solve outputs
   subroutine keep_complex_solve(run, energy, z)
      !> Kernel results
      type(kernel_run), intent(inout) :: run
      !> Energy
      complex(wp), intent(in) :: energy
      !> Buffers holding the solve outputs
      type(complex_work), intent(in) :: z

      run%cx_energy = energy
      run%cx_q = z%q
      run%cx_grad = [z%gv, z%gvn]
   end subroutine keep_complex_solve

   !> Keep the complex weights and tangent, all channels flattened
   !>
   !> @param[in,out] run Kernel results
   !> @param[in]     z   Buffers holding the weights and the tangent
   subroutine keep_complex_derivatives(run, z)
      !> Kernel results
      type(kernel_run), intent(inout) :: run
      !> Buffers holding the weights and the tangent
      type(complex_work), intent(in) :: z

      run%cx_weights = [pack(z%gr, .true.), pack(z%gn, .true.), z%gk1, z%gk2, z%gxi, z%gf, z%ga]
      run%cx_tangent = [pack(z%hr, .true.), pack(z%hn, .true.), pack(z%hk1, .true.), &
         & pack(z%hk2, .true.), pack(z%hxi, .true.), pack(z%hf, .true.), pack(z%ha, .true.), &
         & pack(z%hv, .true.), pack(z%hvn, .true.)]
   end subroutine keep_complex_derivatives

   !> Trace blocks stacked, v over vn, complex
   !>
   !> @param[in] hv  v block (npt, nblk)
   !> @param[in] hvn vn block (npt, nblk)
   pure function stack_complex(hv, hvn) result(h)
      !> v block
      complex(wp), intent(in) :: hv(:, :)
      !> vn block
      complex(wp), intent(in) :: hvn(:, :)
      !> Stacked blocks (2 npt, nblk)
      complex(wp) :: h(2*size(hv, 1), size(hv, 2))

      h(:size(hv, 1), :) = hv
      h(size(hv, 1) + 1:, :) = hvn
   end function stack_complex

   !> Run one model on its kernel of the surface
   !>
   !> @param[in]  model  Model identifier of the fixtures
   !> @param[in]  inp    Kernel inputs
   !> @param[in]  solver Linear solver
   !> @param[in]  eps    Permittivity of the complex phases
   !> @param[out] run    Kernel results
   !> @param[out] error  Library error
   subroutine run_kernel(model, inp, solver, eps, run, error)
      !> Model identifier of the fixtures
      integer, intent(in) :: model
      !> Kernel inputs
      type(kernel_inputs), intent(in) :: inp
      !> Linear solver
      class(pcm_solver), intent(in) :: solver
      !> Permittivity of the complex phases
      complex(wp), intent(in) :: eps
      !> Kernel results
      type(kernel_run), intent(out) :: run
      !> Library error
      type(moist_error_type), allocatable, intent(out) :: error

      allocate (run%stencil(2*inp%npt, size(fd6_offsets)), &
         & run%cx_stencil(2*inp%npt, size(fd6_offsets)))
      select case (model)
      case (pcm_cpcm)
         if (inp%gaussian) then
            call run_conductor_gauss(inp, solver, eps, run, error)
         else
            call run_conductor_point(inp, solver, eps, run, error)
         end if
      case (pcm_iefpcm)
         if (inp%gaussian) then
            call run_iefpcm_gauss(inp, solver, eps, run, error)
         else
            call run_iefpcm_point(inp, solver, eps, run, error)
         end if
      case (pcm_ssvpe)
         if (inp%gaussian) then
            call run_ssvpe_gauss(inp, solver, eps, run, error)
         else
            call run_ssvpe_point(inp, solver, eps, run, error)
         end if
      case (pcm_dpcm)
         if (inp%gaussian) then
            call run_dpcm_gauss(inp, solver, eps, run, error)
         else
            call run_dpcm_point(inp, solver, eps, run, error)
         end if
      case (pcm_dpcm_original)
         if (inp%gaussian) then
            call run_dpcm_original_gauss(inp, solver, eps, run, error)
         else
            call run_dpcm_original_point(inp, solver, eps, run, error)
         end if
      case (pcm_lpb)
         if (inp%gaussian) then
            call run_lpb_gauss(inp, solver, eps, run, error)
         else
            call run_lpb_point(inp, solver, eps, run, error)
         end if
      case (pcm_slic)
         if (inp%gaussian) then
            call run_slic_gauss(inp, solver, eps, run, error)
         else
            call run_slic_point(inp, solver, eps, run, error)
         end if
      case (pcm_nonlocal)
         if (inp%gaussian) then
            call run_nonlocal_gauss(inp, solver, eps, run, error)
         else
            call run_nonlocal_point(inp, solver, eps, run, error)
         end if
      case default
         call fatal_error(error, "unknown PCM model")
      end select
   end subroutine run_kernel

   !> Conductor, Gaussian kernel: real and complex phases and stencils
   !>
   !> @param[in]     inp    Kernel inputs
   !> @param[in]     solver Linear solver
   !> @param[in]     eps    Permittivity of the complex phases
   !> @param[in,out] run    Kernel results
   !> @param[out]    error  Library error
   subroutine run_conductor_gauss(inp, solver, eps, run, error)
      !> Kernel inputs
      type(kernel_inputs), intent(in) :: inp
      !> Linear solver
      class(pcm_solver), intent(in) :: solver
      !> Permittivity of the complex phases
      complex(wp), intent(in) :: eps
      !> Kernel results
      type(kernel_run), intent(inout) :: run
      !> Library error
      type(moist_error_type), allocatable, intent(out) :: error

      !> Kernel state of both variants
      type(conductor_complex_state), target :: st
      !> Output buffers
      type(real_work) :: w
      type(complex_work) :: z
      !> Energies
      real(wp) :: e
      complex(wp) :: ez
      !> Stencil index
      integer :: k

      call new_real_work(inp%npt, w)
      call new_complex_work(inp%npt, z)
      run%has_complex = .true.
      allocate (run%params(1), run%cx_params(1))
      call conductor_assemble(st%conductor_state, inp%r, inp%xi, inp%f, error, solver)
      if (allocated(error)) return
      call conductor_medium(st%conductor_state, feps_of(inp%phys%eps), error)
      if (allocated(error)) return
      call conductor_solve(st%conductor_state, inp%v, e, w%q, w%gv, error)
      if (allocated(error)) return
      call keep_solve(run, e, w)
      call conductor_response(st%conductor_state, inp%dv, w%hv, error)
      if (allocated(error)) return
      run%response = stack(w%hv, w%hvn)
      call conductor_param_weights(st%conductor_state, run%params(1), error)
      if (allocated(error)) return
      call conductor_weights(st%conductor_state, w%gr, w%gxi, w%gf, error)
      if (allocated(error)) return
      call keep_weights(run, w)
      call conductor_tangent(st%conductor_state, inp%dr, inp%dxi, inp%df, inp%dv, w%hr, w%hxi, &
         & w%hf, w%hv, error)
      if (allocated(error)) return
      call keep_tangent(run, w)
      call conductor_tangent(st%conductor_state, inp%zr, inp%zs, inp%zs, inp%dv, w%hr, w%hxi, &
         & w%hf, w%hv, error)
      if (allocated(error)) return
      run%response_tangent = stack(w%hv, w%hvn)

      run%cx_param = cfeps_of(eps)
      call conductor_complex_medium(st, run%cx_param, error)
      if (allocated(error)) return
      call conductor_complex_solve(st, cmplx(inp%v, 0.0_wp, wp), ez, z%q, z%gv, error)
      if (allocated(error)) return
      call keep_complex_solve(run, ez, z)
      call conductor_complex_response(st, cmplx(inp%dv, 0.0_wp, wp), z%hv, error)
      if (allocated(error)) return
      run%cx_response = stack_complex(z%hv, z%hvn)
      call conductor_complex_param_weights(st, run%cx_params(1), error)
      if (allocated(error)) return
      call conductor_complex_weights(st, z%gr, z%gxi, z%gf, error)
      if (allocated(error)) return
      call conductor_complex_tangent(st, inp%dr, inp%dxi, inp%df, cmplx(inp%dv, 0.0_wp, wp), z%hr, &
         & z%hxi, z%hf, z%hv, error)
      if (allocated(error)) return
      call keep_complex_derivatives(run, z)
      call conductor_complex_tangent(st, inp%zr, inp%zs, inp%zs, cmplx(inp%dv, 0.0_wp, wp), z%hr, &
         & z%hxi, z%hf, z%hv, error)
      if (allocated(error)) return
      run%cx_response_tangent = stack_complex(z%hv, z%hvn)

      call conductor_solve(st%conductor_state, inp%v, e, w%q, w%gv, error)
      if (allocated(error)) return
      call keep_again(run, e, w)
      do k = 1, size(fd6_offsets)
         call conductor_solve(st%conductor_state, inp%v + fd6_offsets(k)*inp%h*inp%dv(:, 1), e, &
            & w%q, w%gv, error)
         if (allocated(error)) return
         run%stencil(:, k) = [w%gv, w%gvn]
      end do
      do k = 1, size(fd6_offsets)
         call conductor_complex_solve(st, &
            & cmplx(inp%v + fd6_offsets(k)*inp%h*inp%dv(:, 1), 0.0_wp, wp), ez, z%q, z%gv, error)
         if (allocated(error)) return
         run%cx_stencil(:, k) = [z%gv, z%gvn]
      end do
   end subroutine run_conductor_gauss

   !> Conductor, point kernel: real and complex phases and stencils
   !>
   !> @param[in]     inp    Kernel inputs
   !> @param[in]     solver Linear solver
   !> @param[in]     eps    Permittivity of the complex phases
   !> @param[in,out] run    Kernel results
   !> @param[out]    error  Library error
   subroutine run_conductor_point(inp, solver, eps, run, error)
      !> Kernel inputs
      type(kernel_inputs), intent(in) :: inp
      !> Linear solver
      class(pcm_solver), intent(in) :: solver
      !> Permittivity of the complex phases
      complex(wp), intent(in) :: eps
      !> Kernel results
      type(kernel_run), intent(inout) :: run
      !> Library error
      type(moist_error_type), allocatable, intent(out) :: error

      !> Kernel state of both variants
      type(conductor_point_complex_state), target :: st
      !> Output buffers
      type(real_work) :: w
      type(complex_work) :: z
      !> Energies
      real(wp) :: e
      complex(wp) :: ez
      !> Stencil index
      integer :: k

      call new_real_work(inp%npt, w)
      call new_complex_work(inp%npt, z)
      run%has_complex = .true.
      allocate (run%params(1), run%cx_params(1))
      call conductor_point_assemble(st%conductor_point_state, inp%r, inp%a, error, solver)
      if (allocated(error)) return
      call conductor_point_medium(st%conductor_point_state, feps_of(inp%phys%eps), error)
      if (allocated(error)) return
      call conductor_point_solve(st%conductor_point_state, inp%v, e, w%q, w%gv, error)
      if (allocated(error)) return
      call keep_solve(run, e, w)
      call conductor_point_response(st%conductor_point_state, inp%dv, w%hv, error)
      if (allocated(error)) return
      run%response = stack(w%hv, w%hvn)
      call conductor_point_param_weights(st%conductor_point_state, run%params(1), error)
      if (allocated(error)) return
      call conductor_point_weights(st%conductor_point_state, w%gr, w%ga, error)
      if (allocated(error)) return
      call keep_weights(run, w)
      call conductor_point_tangent(st%conductor_point_state, inp%dr, inp%da, inp%dv, w%hr, w%ha, &
         & w%hv, error)
      if (allocated(error)) return
      call keep_tangent(run, w)
      call conductor_point_tangent(st%conductor_point_state, inp%zr, inp%zs, inp%dv, w%hr, w%ha, &
         & w%hv, error)
      if (allocated(error)) return
      run%response_tangent = stack(w%hv, w%hvn)

      run%cx_param = cfeps_of(eps)
      call conductor_point_complex_medium(st, run%cx_param, error)
      if (allocated(error)) return
      call conductor_point_complex_solve(st, cmplx(inp%v, 0.0_wp, wp), ez, z%q, z%gv, error)
      if (allocated(error)) return
      call keep_complex_solve(run, ez, z)
      call conductor_point_complex_response(st, cmplx(inp%dv, 0.0_wp, wp), z%hv, error)
      if (allocated(error)) return
      run%cx_response = stack_complex(z%hv, z%hvn)
      call conductor_point_complex_param_weights(st, run%cx_params(1), error)
      if (allocated(error)) return
      call conductor_point_complex_weights(st, z%gr, z%ga, error)
      if (allocated(error)) return
      call conductor_point_complex_tangent(st, inp%dr, inp%da, cmplx(inp%dv, 0.0_wp, wp), z%hr, &
         & z%ha, z%hv, error)
      if (allocated(error)) return
      call keep_complex_derivatives(run, z)
      call conductor_point_complex_tangent(st, inp%zr, inp%zs, cmplx(inp%dv, 0.0_wp, wp), z%hr, &
         & z%ha, z%hv, error)
      if (allocated(error)) return
      run%cx_response_tangent = stack_complex(z%hv, z%hvn)

      call conductor_point_solve(st%conductor_point_state, inp%v, e, w%q, w%gv, error)
      if (allocated(error)) return
      call keep_again(run, e, w)
      do k = 1, size(fd6_offsets)
         call conductor_point_solve(st%conductor_point_state, &
            & inp%v + fd6_offsets(k)*inp%h*inp%dv(:, 1), e, w%q, w%gv, error)
         if (allocated(error)) return
         run%stencil(:, k) = [w%gv, w%gvn]
      end do
      do k = 1, size(fd6_offsets)
         call conductor_point_complex_solve(st, &
            & cmplx(inp%v + fd6_offsets(k)*inp%h*inp%dv(:, 1), 0.0_wp, wp), ez, z%q, z%gv, error)
         if (allocated(error)) return
         run%cx_stencil(:, k) = [z%gv, z%gvn]
      end do
   end subroutine run_conductor_point

   !> IEF-PCM, Gaussian kernel: real and complex phases and stencils
   !>
   !> @param[in]     inp    Kernel inputs
   !> @param[in]     solver Linear solver
   !> @param[in]     eps    Permittivity of the complex phases
   !> @param[in,out] run    Kernel results
   !> @param[out]    error  Library error
   subroutine run_iefpcm_gauss(inp, solver, eps, run, error)
      !> Kernel inputs
      type(kernel_inputs), intent(in) :: inp
      !> Linear solver
      class(pcm_solver), intent(in) :: solver
      !> Permittivity of the complex phases
      complex(wp), intent(in) :: eps
      !> Kernel results
      type(kernel_run), intent(inout) :: run
      !> Library error
      type(moist_error_type), allocatable, intent(out) :: error

      !> Kernel state of both variants
      type(iefpcm_complex_state), target :: st
      !> Output buffers
      type(real_work) :: w
      type(complex_work) :: z
      !> Energies
      real(wp) :: e
      complex(wp) :: ez
      !> Stencil index
      integer :: k

      call new_real_work(inp%npt, w)
      call new_complex_work(inp%npt, z)
      run%has_complex = .true.
      allocate (run%params(1), run%cx_params(1))
      call iefpcm_assemble(st%iefpcm_state, inp%r, inp%n, inp%k1, inp%k2, inp%xi, inp%f, inp%a, &
         & error, solver)
      if (allocated(error)) return
      call iefpcm_medium(st%iefpcm_state, inp%phys%eps, error)
      if (allocated(error)) return
      call iefpcm_solve(st%iefpcm_state, inp%v, e, w%q, w%gv, error)
      if (allocated(error)) return
      call keep_solve(run, e, w)
      call iefpcm_response(st%iefpcm_state, inp%dv, w%hv, error)
      if (allocated(error)) return
      run%response = stack(w%hv, w%hvn)
      call iefpcm_param_weights(st%iefpcm_state, run%params(1), error)
      if (allocated(error)) return
      call iefpcm_weights(st%iefpcm_state, w%gr, w%gn, w%gk1, w%gk2, w%gxi, w%gf, w%ga, error)
      if (allocated(error)) return
      call keep_weights(run, w)
      call iefpcm_tangent(st%iefpcm_state, inp%dr, inp%dn, inp%dk1, inp%dk2, inp%dxi, inp%df, &
         & inp%da, inp%dv, w%hr, w%hn, w%hk1, w%hk2, w%hxi, w%hf, w%ha, w%hv, error)
      if (allocated(error)) return
      call keep_tangent(run, w)
      call iefpcm_tangent(st%iefpcm_state, inp%zr, inp%zr, inp%zs, inp%zs, inp%zs, inp%zs, inp%zs, &
         & inp%dv, w%hr, w%hn, w%hk1, w%hk2, w%hxi, w%hf, w%ha, w%hv, error)
      if (allocated(error)) return
      run%response_tangent = stack(w%hv, w%hvn)

      run%cx_param = eps
      call iefpcm_complex_medium(st, run%cx_param, error)
      if (allocated(error)) return
      call iefpcm_complex_solve(st, cmplx(inp%v, 0.0_wp, wp), ez, z%q, z%gv, error)
      if (allocated(error)) return
      call keep_complex_solve(run, ez, z)
      call iefpcm_complex_response(st, cmplx(inp%dv, 0.0_wp, wp), z%hv, error)
      if (allocated(error)) return
      run%cx_response = stack_complex(z%hv, z%hvn)
      call iefpcm_complex_param_weights(st, run%cx_params(1), error)
      if (allocated(error)) return
      call iefpcm_complex_weights(st, z%gr, z%gn, z%gk1, z%gk2, z%gxi, z%gf, z%ga, error)
      if (allocated(error)) return
      call iefpcm_complex_tangent(st, inp%dr, inp%dn, inp%dk1, inp%dk2, inp%dxi, inp%df, inp%da, &
         & cmplx(inp%dv, 0.0_wp, wp), z%hr, z%hn, z%hk1, z%hk2, z%hxi, z%hf, z%ha, z%hv, error)
      if (allocated(error)) return
      call keep_complex_derivatives(run, z)
      call iefpcm_complex_tangent(st, inp%zr, inp%zr, inp%zs, inp%zs, inp%zs, inp%zs, inp%zs, &
         & cmplx(inp%dv, 0.0_wp, wp), z%hr, z%hn, z%hk1, z%hk2, z%hxi, z%hf, z%ha, z%hv, error)
      if (allocated(error)) return
      run%cx_response_tangent = stack_complex(z%hv, z%hvn)

      call iefpcm_solve(st%iefpcm_state, inp%v, e, w%q, w%gv, error)
      if (allocated(error)) return
      call keep_again(run, e, w)
      do k = 1, size(fd6_offsets)
         call iefpcm_solve(st%iefpcm_state, inp%v + fd6_offsets(k)*inp%h*inp%dv(:, 1), e, w%q, &
            & w%gv, error)
         if (allocated(error)) return
         run%stencil(:, k) = [w%gv, w%gvn]
      end do
      do k = 1, size(fd6_offsets)
         call iefpcm_complex_solve(st, &
            & cmplx(inp%v + fd6_offsets(k)*inp%h*inp%dv(:, 1), 0.0_wp, wp), ez, z%q, z%gv, error)
         if (allocated(error)) return
         run%cx_stencil(:, k) = [z%gv, z%gvn]
      end do
   end subroutine run_iefpcm_gauss

   !> IEF-PCM, point kernel: real and complex phases and stencils
   !>
   !> @param[in]     inp    Kernel inputs
   !> @param[in]     solver Linear solver
   !> @param[in]     eps    Permittivity of the complex phases
   !> @param[in,out] run    Kernel results
   !> @param[out]    error  Library error
   subroutine run_iefpcm_point(inp, solver, eps, run, error)
      !> Kernel inputs
      type(kernel_inputs), intent(in) :: inp
      !> Linear solver
      class(pcm_solver), intent(in) :: solver
      !> Permittivity of the complex phases
      complex(wp), intent(in) :: eps
      !> Kernel results
      type(kernel_run), intent(inout) :: run
      !> Library error
      type(moist_error_type), allocatable, intent(out) :: error

      !> Kernel state of both variants
      type(iefpcm_point_complex_state), target :: st
      !> Output buffers
      type(real_work) :: w
      type(complex_work) :: z
      !> Energies
      real(wp) :: e
      complex(wp) :: ez
      !> Stencil index
      integer :: k

      call new_real_work(inp%npt, w)
      call new_complex_work(inp%npt, z)
      run%has_complex = .true.
      allocate (run%params(1), run%cx_params(1))
      call iefpcm_point_assemble(st%iefpcm_point_state, inp%r, inp%n, inp%k1, inp%k2, inp%a, &
         & error, solver)
      if (allocated(error)) return
      call iefpcm_point_medium(st%iefpcm_point_state, inp%phys%eps, error)
      if (allocated(error)) return
      call iefpcm_point_solve(st%iefpcm_point_state, inp%v, e, w%q, w%gv, error)
      if (allocated(error)) return
      call keep_solve(run, e, w)
      call iefpcm_point_response(st%iefpcm_point_state, inp%dv, w%hv, error)
      if (allocated(error)) return
      run%response = stack(w%hv, w%hvn)
      call iefpcm_point_param_weights(st%iefpcm_point_state, run%params(1), error)
      if (allocated(error)) return
      call iefpcm_point_weights(st%iefpcm_point_state, w%gr, w%gn, w%gk1, w%gk2, w%ga, error)
      if (allocated(error)) return
      call keep_weights(run, w)
      call iefpcm_point_tangent(st%iefpcm_point_state, inp%dr, inp%dn, inp%dk1, inp%dk2, inp%da, &
         & inp%dv, w%hr, w%hn, w%hk1, w%hk2, w%ha, w%hv, error)
      if (allocated(error)) return
      call keep_tangent(run, w)
      call iefpcm_point_tangent(st%iefpcm_point_state, inp%zr, inp%zr, inp%zs, inp%zs, inp%zs, &
         & inp%dv, w%hr, w%hn, w%hk1, w%hk2, w%ha, w%hv, error)
      if (allocated(error)) return
      run%response_tangent = stack(w%hv, w%hvn)

      run%cx_param = eps
      call iefpcm_point_complex_medium(st, run%cx_param, error)
      if (allocated(error)) return
      call iefpcm_point_complex_solve(st, cmplx(inp%v, 0.0_wp, wp), ez, z%q, z%gv, error)
      if (allocated(error)) return
      call keep_complex_solve(run, ez, z)
      call iefpcm_point_complex_response(st, cmplx(inp%dv, 0.0_wp, wp), z%hv, error)
      if (allocated(error)) return
      run%cx_response = stack_complex(z%hv, z%hvn)
      call iefpcm_point_complex_param_weights(st, run%cx_params(1), error)
      if (allocated(error)) return
      call iefpcm_point_complex_weights(st, z%gr, z%gn, z%gk1, z%gk2, z%ga, error)
      if (allocated(error)) return
      call iefpcm_point_complex_tangent(st, inp%dr, inp%dn, inp%dk1, inp%dk2, inp%da, &
         & cmplx(inp%dv, 0.0_wp, wp), z%hr, z%hn, z%hk1, z%hk2, z%ha, z%hv, error)
      if (allocated(error)) return
      call keep_complex_derivatives(run, z)
      call iefpcm_point_complex_tangent(st, inp%zr, inp%zr, inp%zs, inp%zs, inp%zs, &
         & cmplx(inp%dv, 0.0_wp, wp), z%hr, z%hn, z%hk1, z%hk2, z%ha, z%hv, error)
      if (allocated(error)) return
      run%cx_response_tangent = stack_complex(z%hv, z%hvn)

      call iefpcm_point_solve(st%iefpcm_point_state, inp%v, e, w%q, w%gv, error)
      if (allocated(error)) return
      call keep_again(run, e, w)
      do k = 1, size(fd6_offsets)
         call iefpcm_point_solve(st%iefpcm_point_state, inp%v + fd6_offsets(k)*inp%h*inp%dv(:, 1), &
            & e, w%q, w%gv, error)
         if (allocated(error)) return
         run%stencil(:, k) = [w%gv, w%gvn]
      end do
      do k = 1, size(fd6_offsets)
         call iefpcm_point_complex_solve(st, &
            & cmplx(inp%v + fd6_offsets(k)*inp%h*inp%dv(:, 1), 0.0_wp, wp), ez, z%q, z%gv, error)
         if (allocated(error)) return
         run%cx_stencil(:, k) = [z%gv, z%gvn]
      end do
   end subroutine run_iefpcm_point

   !> SS(V)PE, Gaussian kernel: real and complex phases and stencils
   !>
   !> @param[in]     inp    Kernel inputs
   !> @param[in]     solver Linear solver
   !> @param[in]     eps    Permittivity of the complex phases
   !> @param[in,out] run    Kernel results
   !> @param[out]    error  Library error
   subroutine run_ssvpe_gauss(inp, solver, eps, run, error)
      !> Kernel inputs
      type(kernel_inputs), intent(in) :: inp
      !> Linear solver
      class(pcm_solver), intent(in) :: solver
      !> Permittivity of the complex phases
      complex(wp), intent(in) :: eps
      !> Kernel results
      type(kernel_run), intent(inout) :: run
      !> Library error
      type(moist_error_type), allocatable, intent(out) :: error

      !> Kernel state of both variants
      type(ssvpe_complex_state), target :: st
      !> Output buffers
      type(real_work) :: w
      type(complex_work) :: z
      !> Energies
      real(wp) :: e
      complex(wp) :: ez
      !> Stencil index
      integer :: k

      call new_real_work(inp%npt, w)
      call new_complex_work(inp%npt, z)
      run%has_complex = .true.
      allocate (run%params(1), run%cx_params(1))
      call ssvpe_assemble(st%ssvpe_state, inp%r, inp%n, inp%k1, inp%k2, inp%xi, inp%f, inp%a, &
         & error, solver)
      if (allocated(error)) return
      call ssvpe_medium(st%ssvpe_state, inp%phys%eps, error)
      if (allocated(error)) return
      call ssvpe_solve(st%ssvpe_state, inp%v, e, w%q, w%gv, error)
      if (allocated(error)) return
      call keep_solve(run, e, w)
      call ssvpe_response(st%ssvpe_state, inp%dv, w%hv, error)
      if (allocated(error)) return
      run%response = stack(w%hv, w%hvn)
      call ssvpe_param_weights(st%ssvpe_state, run%params(1), error)
      if (allocated(error)) return
      call ssvpe_weights(st%ssvpe_state, w%gr, w%gn, w%gk1, w%gk2, w%gxi, w%gf, w%ga, error)
      if (allocated(error)) return
      call keep_weights(run, w)
      call ssvpe_tangent(st%ssvpe_state, inp%dr, inp%dn, inp%dk1, inp%dk2, inp%dxi, inp%df, &
         & inp%da, inp%dv, w%hr, w%hn, w%hk1, w%hk2, w%hxi, w%hf, w%ha, w%hv, error)
      if (allocated(error)) return
      call keep_tangent(run, w)
      call ssvpe_tangent(st%ssvpe_state, inp%zr, inp%zr, inp%zs, inp%zs, inp%zs, inp%zs, inp%zs, &
         & inp%dv, w%hr, w%hn, w%hk1, w%hk2, w%hxi, w%hf, w%ha, w%hv, error)
      if (allocated(error)) return
      run%response_tangent = stack(w%hv, w%hvn)

      run%cx_param = eps
      call ssvpe_complex_medium(st, run%cx_param, error)
      if (allocated(error)) return
      call ssvpe_complex_solve(st, cmplx(inp%v, 0.0_wp, wp), ez, z%q, z%gv, error)
      if (allocated(error)) return
      call keep_complex_solve(run, ez, z)
      call ssvpe_complex_response(st, cmplx(inp%dv, 0.0_wp, wp), z%hv, error)
      if (allocated(error)) return
      run%cx_response = stack_complex(z%hv, z%hvn)
      call ssvpe_complex_param_weights(st, run%cx_params(1), error)
      if (allocated(error)) return
      call ssvpe_complex_weights(st, z%gr, z%gn, z%gk1, z%gk2, z%gxi, z%gf, z%ga, error)
      if (allocated(error)) return
      call ssvpe_complex_tangent(st, inp%dr, inp%dn, inp%dk1, inp%dk2, inp%dxi, inp%df, inp%da, &
         & cmplx(inp%dv, 0.0_wp, wp), z%hr, z%hn, z%hk1, z%hk2, z%hxi, z%hf, z%ha, z%hv, error)
      if (allocated(error)) return
      call keep_complex_derivatives(run, z)
      call ssvpe_complex_tangent(st, inp%zr, inp%zr, inp%zs, inp%zs, inp%zs, inp%zs, inp%zs, &
         & cmplx(inp%dv, 0.0_wp, wp), z%hr, z%hn, z%hk1, z%hk2, z%hxi, z%hf, z%ha, z%hv, error)
      if (allocated(error)) return
      run%cx_response_tangent = stack_complex(z%hv, z%hvn)

      call ssvpe_solve(st%ssvpe_state, inp%v, e, w%q, w%gv, error)
      if (allocated(error)) return
      call keep_again(run, e, w)
      do k = 1, size(fd6_offsets)
         call ssvpe_solve(st%ssvpe_state, inp%v + fd6_offsets(k)*inp%h*inp%dv(:, 1), e, w%q, w%gv, &
            & error)
         if (allocated(error)) return
         run%stencil(:, k) = [w%gv, w%gvn]
      end do
      do k = 1, size(fd6_offsets)
         call ssvpe_complex_solve(st, &
            & cmplx(inp%v + fd6_offsets(k)*inp%h*inp%dv(:, 1), 0.0_wp, wp), ez, z%q, z%gv, error)
         if (allocated(error)) return
         run%cx_stencil(:, k) = [z%gv, z%gvn]
      end do
   end subroutine run_ssvpe_gauss

   !> SS(V)PE, point kernel: real and complex phases and stencils
   !>
   !> @param[in]     inp    Kernel inputs
   !> @param[in]     solver Linear solver
   !> @param[in]     eps    Permittivity of the complex phases
   !> @param[in,out] run    Kernel results
   !> @param[out]    error  Library error
   subroutine run_ssvpe_point(inp, solver, eps, run, error)
      !> Kernel inputs
      type(kernel_inputs), intent(in) :: inp
      !> Linear solver
      class(pcm_solver), intent(in) :: solver
      !> Permittivity of the complex phases
      complex(wp), intent(in) :: eps
      !> Kernel results
      type(kernel_run), intent(inout) :: run
      !> Library error
      type(moist_error_type), allocatable, intent(out) :: error

      !> Kernel state of both variants
      type(ssvpe_point_complex_state), target :: st
      !> Output buffers
      type(real_work) :: w
      type(complex_work) :: z
      !> Energies
      real(wp) :: e
      complex(wp) :: ez
      !> Stencil index
      integer :: k

      call new_real_work(inp%npt, w)
      call new_complex_work(inp%npt, z)
      run%has_complex = .true.
      allocate (run%params(1), run%cx_params(1))
      call ssvpe_point_assemble(st%ssvpe_point_state, inp%r, inp%n, inp%k1, inp%k2, inp%a, error, &
         & solver)
      if (allocated(error)) return
      call ssvpe_point_medium(st%ssvpe_point_state, inp%phys%eps, error)
      if (allocated(error)) return
      call ssvpe_point_solve(st%ssvpe_point_state, inp%v, e, w%q, w%gv, error)
      if (allocated(error)) return
      call keep_solve(run, e, w)
      call ssvpe_point_response(st%ssvpe_point_state, inp%dv, w%hv, error)
      if (allocated(error)) return
      run%response = stack(w%hv, w%hvn)
      call ssvpe_point_param_weights(st%ssvpe_point_state, run%params(1), error)
      if (allocated(error)) return
      call ssvpe_point_weights(st%ssvpe_point_state, w%gr, w%gn, w%gk1, w%gk2, w%ga, error)
      if (allocated(error)) return
      call keep_weights(run, w)
      call ssvpe_point_tangent(st%ssvpe_point_state, inp%dr, inp%dn, inp%dk1, inp%dk2, inp%da, &
         & inp%dv, w%hr, w%hn, w%hk1, w%hk2, w%ha, w%hv, error)
      if (allocated(error)) return
      call keep_tangent(run, w)
      call ssvpe_point_tangent(st%ssvpe_point_state, inp%zr, inp%zr, inp%zs, inp%zs, inp%zs, &
         & inp%dv, w%hr, w%hn, w%hk1, w%hk2, w%ha, w%hv, error)
      if (allocated(error)) return
      run%response_tangent = stack(w%hv, w%hvn)

      run%cx_param = eps
      call ssvpe_point_complex_medium(st, run%cx_param, error)
      if (allocated(error)) return
      call ssvpe_point_complex_solve(st, cmplx(inp%v, 0.0_wp, wp), ez, z%q, z%gv, error)
      if (allocated(error)) return
      call keep_complex_solve(run, ez, z)
      call ssvpe_point_complex_response(st, cmplx(inp%dv, 0.0_wp, wp), z%hv, error)
      if (allocated(error)) return
      run%cx_response = stack_complex(z%hv, z%hvn)
      call ssvpe_point_complex_param_weights(st, run%cx_params(1), error)
      if (allocated(error)) return
      call ssvpe_point_complex_weights(st, z%gr, z%gn, z%gk1, z%gk2, z%ga, error)
      if (allocated(error)) return
      call ssvpe_point_complex_tangent(st, inp%dr, inp%dn, inp%dk1, inp%dk2, inp%da, &
         & cmplx(inp%dv, 0.0_wp, wp), z%hr, z%hn, z%hk1, z%hk2, z%ha, z%hv, error)
      if (allocated(error)) return
      call keep_complex_derivatives(run, z)
      call ssvpe_point_complex_tangent(st, inp%zr, inp%zr, inp%zs, inp%zs, inp%zs, &
         & cmplx(inp%dv, 0.0_wp, wp), z%hr, z%hn, z%hk1, z%hk2, z%ha, z%hv, error)
      if (allocated(error)) return
      run%cx_response_tangent = stack_complex(z%hv, z%hvn)

      call ssvpe_point_solve(st%ssvpe_point_state, inp%v, e, w%q, w%gv, error)
      if (allocated(error)) return
      call keep_again(run, e, w)
      do k = 1, size(fd6_offsets)
         call ssvpe_point_solve(st%ssvpe_point_state, inp%v + fd6_offsets(k)*inp%h*inp%dv(:, 1), &
            & e, w%q, w%gv, error)
         if (allocated(error)) return
         run%stencil(:, k) = [w%gv, w%gvn]
      end do
      do k = 1, size(fd6_offsets)
         call ssvpe_point_complex_solve(st, &
            & cmplx(inp%v + fd6_offsets(k)*inp%h*inp%dv(:, 1), 0.0_wp, wp), ez, z%q, z%gv, error)
         if (allocated(error)) return
         run%cx_stencil(:, k) = [z%gv, z%gvn]
      end do
   end subroutine run_ssvpe_point

   !> DPCM, Gaussian kernel: real and complex phases and stencils
   !>
   !> @param[in]     inp    Kernel inputs
   !> @param[in]     solver Linear solver
   !> @param[in]     eps    Permittivity of the complex phases
   !> @param[in,out] run    Kernel results
   !> @param[out]    error  Library error
   subroutine run_dpcm_gauss(inp, solver, eps, run, error)
      !> Kernel inputs
      type(kernel_inputs), intent(in) :: inp
      !> Linear solver
      class(pcm_solver), intent(in) :: solver
      !> Permittivity of the complex phases
      complex(wp), intent(in) :: eps
      !> Kernel results
      type(kernel_run), intent(inout) :: run
      !> Library error
      type(moist_error_type), allocatable, intent(out) :: error

      !> Kernel state of both variants
      type(dpcm_complex_state), target :: st
      !> Output buffers
      type(real_work) :: w
      type(complex_work) :: z
      !> Energies
      real(wp) :: e
      complex(wp) :: ez
      !> Stencil index
      integer :: k

      call new_real_work(inp%npt, w)
      call new_complex_work(inp%npt, z)
      run%has_complex = .true.
      allocate (run%params(1), run%cx_params(1))
      call dpcm_assemble(st%dpcm_state, inp%r, inp%n, inp%k1, inp%k2, inp%xi, inp%a, error, solver)
      if (allocated(error)) return
      call dpcm_medium(st%dpcm_state, inp%phys%eps, error)
      if (allocated(error)) return
      call dpcm_solve(st%dpcm_state, inp%v, inp%vn, e, w%q, w%gv, w%gvn, error)
      if (allocated(error)) return
      call keep_solve(run, e, w)
      call dpcm_response(st%dpcm_state, inp%dv, inp%dvn, w%hv, w%hvn, error)
      if (allocated(error)) return
      run%response = stack(w%hv, w%hvn)
      call dpcm_param_weights(st%dpcm_state, run%params(1), error)
      if (allocated(error)) return
      call dpcm_weights(st%dpcm_state, w%gr, w%gn, w%gk1, w%gk2, w%gxi, w%ga, error)
      if (allocated(error)) return
      call keep_weights(run, w)
      call dpcm_tangent(st%dpcm_state, inp%dr, inp%dn, inp%dk1, inp%dk2, inp%dxi, inp%da, inp%dv, &
         & inp%dvn, w%hr, w%hn, w%hk1, w%hk2, w%hxi, w%ha, w%hv, w%hvn, error)
      if (allocated(error)) return
      call keep_tangent(run, w)
      call dpcm_tangent(st%dpcm_state, inp%zr, inp%zr, inp%zs, inp%zs, inp%zs, inp%zs, inp%dv, &
         & inp%dvn, w%hr, w%hn, w%hk1, w%hk2, w%hxi, w%ha, w%hv, w%hvn, error)
      if (allocated(error)) return
      run%response_tangent = stack(w%hv, w%hvn)

      run%cx_param = eps
      call dpcm_complex_medium(st, run%cx_param, error)
      if (allocated(error)) return
      call dpcm_complex_solve(st, cmplx(inp%v, 0.0_wp, wp), cmplx(inp%vn, 0.0_wp, wp), ez, z%q, &
         & z%gv, z%gvn, error)
      if (allocated(error)) return
      call keep_complex_solve(run, ez, z)
      call dpcm_complex_response(st, cmplx(inp%dv, 0.0_wp, wp), cmplx(inp%dvn, 0.0_wp, wp), z%hv, &
         & z%hvn, error)
      if (allocated(error)) return
      run%cx_response = stack_complex(z%hv, z%hvn)
      call dpcm_complex_param_weights(st, run%cx_params(1), error)
      if (allocated(error)) return
      call dpcm_complex_weights(st, z%gr, z%gn, z%gk1, z%gk2, z%gxi, z%ga, error)
      if (allocated(error)) return
      call dpcm_complex_tangent(st, inp%dr, inp%dn, inp%dk1, inp%dk2, inp%dxi, inp%da, &
         & cmplx(inp%dv, 0.0_wp, wp), cmplx(inp%dvn, 0.0_wp, wp), z%hr, z%hn, z%hk1, z%hk2, z%hxi, &
         & z%ha, z%hv, z%hvn, error)
      if (allocated(error)) return
      call keep_complex_derivatives(run, z)
      call dpcm_complex_tangent(st, inp%zr, inp%zr, inp%zs, inp%zs, inp%zs, inp%zs, &
         & cmplx(inp%dv, 0.0_wp, wp), cmplx(inp%dvn, 0.0_wp, wp), z%hr, z%hn, z%hk1, z%hk2, z%hxi, &
         & z%ha, z%hv, z%hvn, error)
      if (allocated(error)) return
      run%cx_response_tangent = stack_complex(z%hv, z%hvn)

      call dpcm_solve(st%dpcm_state, inp%v, inp%vn, e, w%q, w%gv, w%gvn, error)
      if (allocated(error)) return
      call keep_again(run, e, w)
      do k = 1, size(fd6_offsets)
         call dpcm_solve(st%dpcm_state, inp%v + fd6_offsets(k)*inp%h*inp%dv(:, 1), &
            & inp%vn + fd6_offsets(k)*inp%h*inp%dvn(:, 1), e, w%q, w%gv, w%gvn, error)
         if (allocated(error)) return
         run%stencil(:, k) = [w%gv, w%gvn]
      end do
      do k = 1, size(fd6_offsets)
         call dpcm_complex_solve(st, cmplx(inp%v + fd6_offsets(k)*inp%h*inp%dv(:, 1), 0.0_wp, wp), &
            & cmplx(inp%vn + fd6_offsets(k)*inp%h*inp%dvn(:, 1), 0.0_wp, wp), ez, z%q, z%gv, &
            & z%gvn, error)
         if (allocated(error)) return
         run%cx_stencil(:, k) = [z%gv, z%gvn]
      end do
   end subroutine run_dpcm_gauss

   !> DPCM, point kernel: real and complex phases and stencils
   !>
   !> @param[in]     inp    Kernel inputs
   !> @param[in]     solver Linear solver
   !> @param[in]     eps    Permittivity of the complex phases
   !> @param[in,out] run    Kernel results
   !> @param[out]    error  Library error
   subroutine run_dpcm_point(inp, solver, eps, run, error)
      !> Kernel inputs
      type(kernel_inputs), intent(in) :: inp
      !> Linear solver
      class(pcm_solver), intent(in) :: solver
      !> Permittivity of the complex phases
      complex(wp), intent(in) :: eps
      !> Kernel results
      type(kernel_run), intent(inout) :: run
      !> Library error
      type(moist_error_type), allocatable, intent(out) :: error

      !> Kernel state of both variants
      type(dpcm_point_complex_state), target :: st
      !> Output buffers
      type(real_work) :: w
      type(complex_work) :: z
      !> Energies
      real(wp) :: e
      complex(wp) :: ez
      !> Stencil index
      integer :: k

      call new_real_work(inp%npt, w)
      call new_complex_work(inp%npt, z)
      run%has_complex = .true.
      allocate (run%params(1), run%cx_params(1))
      call dpcm_point_assemble(st%dpcm_point_state, inp%r, inp%n, inp%k1, inp%k2, inp%a, error, &
         & solver)
      if (allocated(error)) return
      call dpcm_point_medium(st%dpcm_point_state, inp%phys%eps, error)
      if (allocated(error)) return
      call dpcm_point_solve(st%dpcm_point_state, inp%v, inp%vn, e, w%q, w%gv, w%gvn, error)
      if (allocated(error)) return
      call keep_solve(run, e, w)
      call dpcm_point_response(st%dpcm_point_state, inp%dv, inp%dvn, w%hv, w%hvn, error)
      if (allocated(error)) return
      run%response = stack(w%hv, w%hvn)
      call dpcm_point_param_weights(st%dpcm_point_state, run%params(1), error)
      if (allocated(error)) return
      call dpcm_point_weights(st%dpcm_point_state, w%gr, w%gn, w%gk1, w%gk2, w%ga, error)
      if (allocated(error)) return
      call keep_weights(run, w)
      call dpcm_point_tangent(st%dpcm_point_state, inp%dr, inp%dn, inp%dk1, inp%dk2, inp%da, &
         & inp%dv, inp%dvn, w%hr, w%hn, w%hk1, w%hk2, w%ha, w%hv, w%hvn, error)
      if (allocated(error)) return
      call keep_tangent(run, w)
      call dpcm_point_tangent(st%dpcm_point_state, inp%zr, inp%zr, inp%zs, inp%zs, inp%zs, inp%dv, &
         & inp%dvn, w%hr, w%hn, w%hk1, w%hk2, w%ha, w%hv, w%hvn, error)
      if (allocated(error)) return
      run%response_tangent = stack(w%hv, w%hvn)

      run%cx_param = eps
      call dpcm_point_complex_medium(st, run%cx_param, error)
      if (allocated(error)) return
      call dpcm_point_complex_solve(st, cmplx(inp%v, 0.0_wp, wp), cmplx(inp%vn, 0.0_wp, wp), ez, &
         & z%q, z%gv, z%gvn, error)
      if (allocated(error)) return
      call keep_complex_solve(run, ez, z)
      call dpcm_point_complex_response(st, cmplx(inp%dv, 0.0_wp, wp), cmplx(inp%dvn, 0.0_wp, wp), &
         & z%hv, z%hvn, error)
      if (allocated(error)) return
      run%cx_response = stack_complex(z%hv, z%hvn)
      call dpcm_point_complex_param_weights(st, run%cx_params(1), error)
      if (allocated(error)) return
      call dpcm_point_complex_weights(st, z%gr, z%gn, z%gk1, z%gk2, z%ga, error)
      if (allocated(error)) return
      call dpcm_point_complex_tangent(st, inp%dr, inp%dn, inp%dk1, inp%dk2, inp%da, &
         & cmplx(inp%dv, 0.0_wp, wp), cmplx(inp%dvn, 0.0_wp, wp), z%hr, z%hn, z%hk1, z%hk2, z%ha, &
         & z%hv, z%hvn, error)
      if (allocated(error)) return
      call keep_complex_derivatives(run, z)
      call dpcm_point_complex_tangent(st, inp%zr, inp%zr, inp%zs, inp%zs, inp%zs, &
         & cmplx(inp%dv, 0.0_wp, wp), cmplx(inp%dvn, 0.0_wp, wp), z%hr, z%hn, z%hk1, z%hk2, z%ha, &
         & z%hv, z%hvn, error)
      if (allocated(error)) return
      run%cx_response_tangent = stack_complex(z%hv, z%hvn)

      call dpcm_point_solve(st%dpcm_point_state, inp%v, inp%vn, e, w%q, w%gv, w%gvn, error)
      if (allocated(error)) return
      call keep_again(run, e, w)
      do k = 1, size(fd6_offsets)
         call dpcm_point_solve(st%dpcm_point_state, inp%v + fd6_offsets(k)*inp%h*inp%dv(:, 1), &
            & inp%vn + fd6_offsets(k)*inp%h*inp%dvn(:, 1), e, w%q, w%gv, w%gvn, error)
         if (allocated(error)) return
         run%stencil(:, k) = [w%gv, w%gvn]
      end do
      do k = 1, size(fd6_offsets)
         call dpcm_point_complex_solve(st, &
            & cmplx(inp%v + fd6_offsets(k)*inp%h*inp%dv(:, 1), 0.0_wp, wp), &
            & cmplx(inp%vn + fd6_offsets(k)*inp%h*inp%dvn(:, 1), 0.0_wp, wp), ez, z%q, z%gv, &
            & z%gvn, error)
         if (allocated(error)) return
         run%cx_stencil(:, k) = [z%gv, z%gvn]
      end do
   end subroutine run_dpcm_point

   !> Original DPCM, Gaussian kernel: real and complex phases and stencils
   !>
   !> @param[in]     inp    Kernel inputs
   !> @param[in]     solver Linear solver
   !> @param[in]     eps    Permittivity of the complex phases
   !> @param[in,out] run    Kernel results
   !> @param[out]    error  Library error
   subroutine run_dpcm_original_gauss(inp, solver, eps, run, error)
      !> Kernel inputs
      type(kernel_inputs), intent(in) :: inp
      !> Linear solver
      class(pcm_solver), intent(in) :: solver
      !> Permittivity of the complex phases
      complex(wp), intent(in) :: eps
      !> Kernel results
      type(kernel_run), intent(inout) :: run
      !> Library error
      type(moist_error_type), allocatable, intent(out) :: error

      !> Kernel state of both variants
      type(dpcm_original_complex_state), target :: st
      !> Output buffers
      type(real_work) :: w
      type(complex_work) :: z
      !> Energies
      real(wp) :: e
      complex(wp) :: ez
      !> Stencil index
      integer :: k

      call new_real_work(inp%npt, w)
      call new_complex_work(inp%npt, z)
      run%has_complex = .true.
      allocate (run%params(1), run%cx_params(1))
      call dpcm_original_assemble(st%dpcm_original_state, inp%r, inp%n, inp%k1, inp%k2, inp%xi, &
         & inp%a, error, solver)
      if (allocated(error)) return
      call dpcm_original_medium(st%dpcm_original_state, inp%phys%eps, error)
      if (allocated(error)) return
      call dpcm_original_solve(st%dpcm_original_state, inp%v, inp%vn, e, w%q, w%gv, w%gvn, error)
      if (allocated(error)) return
      call keep_solve(run, e, w)
      call dpcm_original_response(st%dpcm_original_state, inp%dv, inp%dvn, w%hv, w%hvn, error)
      if (allocated(error)) return
      run%response = stack(w%hv, w%hvn)
      call dpcm_original_param_weights(st%dpcm_original_state, run%params(1), error)
      if (allocated(error)) return
      call dpcm_original_weights(st%dpcm_original_state, w%gr, w%gn, w%gk1, w%gk2, w%gxi, w%ga, &
         & error)
      if (allocated(error)) return
      call keep_weights(run, w)
      call dpcm_original_tangent(st%dpcm_original_state, inp%dr, inp%dn, inp%dk1, inp%dk2, &
         & inp%dxi, inp%da, inp%dv, inp%dvn, w%hr, w%hn, w%hk1, w%hk2, w%hxi, w%ha, w%hv, w%hvn, &
         & error)
      if (allocated(error)) return
      call keep_tangent(run, w)
      call dpcm_original_tangent(st%dpcm_original_state, inp%zr, inp%zr, inp%zs, inp%zs, inp%zs, &
         & inp%zs, inp%dv, inp%dvn, w%hr, w%hn, w%hk1, w%hk2, w%hxi, w%ha, w%hv, w%hvn, error)
      if (allocated(error)) return
      run%response_tangent = stack(w%hv, w%hvn)

      run%cx_param = eps
      call dpcm_original_complex_medium(st, run%cx_param, error)
      if (allocated(error)) return
      call dpcm_original_complex_solve(st, cmplx(inp%v, 0.0_wp, wp), cmplx(inp%vn, 0.0_wp, wp), &
         & ez, z%q, z%gv, z%gvn, error)
      if (allocated(error)) return
      call keep_complex_solve(run, ez, z)
      call dpcm_original_complex_response(st, cmplx(inp%dv, 0.0_wp, wp), &
         & cmplx(inp%dvn, 0.0_wp, wp), z%hv, z%hvn, error)
      if (allocated(error)) return
      run%cx_response = stack_complex(z%hv, z%hvn)
      call dpcm_original_complex_param_weights(st, run%cx_params(1), error)
      if (allocated(error)) return
      call dpcm_original_complex_weights(st, z%gr, z%gn, z%gk1, z%gk2, z%gxi, z%ga, error)
      if (allocated(error)) return
      call dpcm_original_complex_tangent(st, inp%dr, inp%dn, inp%dk1, inp%dk2, inp%dxi, inp%da, &
         & cmplx(inp%dv, 0.0_wp, wp), cmplx(inp%dvn, 0.0_wp, wp), z%hr, z%hn, z%hk1, z%hk2, z%hxi, &
         & z%ha, z%hv, z%hvn, error)
      if (allocated(error)) return
      call keep_complex_derivatives(run, z)
      call dpcm_original_complex_tangent(st, inp%zr, inp%zr, inp%zs, inp%zs, inp%zs, inp%zs, &
         & cmplx(inp%dv, 0.0_wp, wp), cmplx(inp%dvn, 0.0_wp, wp), z%hr, z%hn, z%hk1, z%hk2, z%hxi, &
         & z%ha, z%hv, z%hvn, error)
      if (allocated(error)) return
      run%cx_response_tangent = stack_complex(z%hv, z%hvn)

      call dpcm_original_solve(st%dpcm_original_state, inp%v, inp%vn, e, w%q, w%gv, w%gvn, error)
      if (allocated(error)) return
      call keep_again(run, e, w)
      do k = 1, size(fd6_offsets)
         call dpcm_original_solve(st%dpcm_original_state, &
            & inp%v + fd6_offsets(k)*inp%h*inp%dv(:, 1), &
            & inp%vn + fd6_offsets(k)*inp%h*inp%dvn(:, 1), e, w%q, w%gv, w%gvn, error)
         if (allocated(error)) return
         run%stencil(:, k) = [w%gv, w%gvn]
      end do
      do k = 1, size(fd6_offsets)
         call dpcm_original_complex_solve(st, &
            & cmplx(inp%v + fd6_offsets(k)*inp%h*inp%dv(:, 1), 0.0_wp, wp), &
            & cmplx(inp%vn + fd6_offsets(k)*inp%h*inp%dvn(:, 1), 0.0_wp, wp), ez, z%q, z%gv, &
            & z%gvn, error)
         if (allocated(error)) return
         run%cx_stencil(:, k) = [z%gv, z%gvn]
      end do
   end subroutine run_dpcm_original_gauss

   !> Original DPCM, point kernel: real and complex phases and stencils
   !>
   !> @param[in]     inp    Kernel inputs
   !> @param[in]     solver Linear solver
   !> @param[in]     eps    Permittivity of the complex phases
   !> @param[in,out] run    Kernel results
   !> @param[out]    error  Library error
   subroutine run_dpcm_original_point(inp, solver, eps, run, error)
      !> Kernel inputs
      type(kernel_inputs), intent(in) :: inp
      !> Linear solver
      class(pcm_solver), intent(in) :: solver
      !> Permittivity of the complex phases
      complex(wp), intent(in) :: eps
      !> Kernel results
      type(kernel_run), intent(inout) :: run
      !> Library error
      type(moist_error_type), allocatable, intent(out) :: error

      !> Kernel state of both variants
      type(dpcm_original_point_complex_state), target :: st
      !> Output buffers
      type(real_work) :: w
      type(complex_work) :: z
      !> Energies
      real(wp) :: e
      complex(wp) :: ez
      !> Stencil index
      integer :: k

      call new_real_work(inp%npt, w)
      call new_complex_work(inp%npt, z)
      run%has_complex = .true.
      allocate (run%params(1), run%cx_params(1))
      call dpcm_original_point_assemble(st%dpcm_original_point_state, inp%r, inp%n, inp%k1, &
         & inp%k2, inp%a, error, solver)
      if (allocated(error)) return
      call dpcm_original_point_medium(st%dpcm_original_point_state, inp%phys%eps, error)
      if (allocated(error)) return
      call dpcm_original_point_solve(st%dpcm_original_point_state, inp%v, inp%vn, e, w%q, w%gv, &
         & w%gvn, error)
      if (allocated(error)) return
      call keep_solve(run, e, w)
      call dpcm_original_point_response(st%dpcm_original_point_state, inp%dv, inp%dvn, w%hv, &
         & w%hvn, error)
      if (allocated(error)) return
      run%response = stack(w%hv, w%hvn)
      call dpcm_original_point_param_weights(st%dpcm_original_point_state, run%params(1), error)
      if (allocated(error)) return
      call dpcm_original_point_weights(st%dpcm_original_point_state, w%gr, w%gn, w%gk1, w%gk2, &
         & w%ga, error)
      if (allocated(error)) return
      call keep_weights(run, w)
      call dpcm_original_point_tangent(st%dpcm_original_point_state, inp%dr, inp%dn, inp%dk1, &
         & inp%dk2, inp%da, inp%dv, inp%dvn, w%hr, w%hn, w%hk1, w%hk2, w%ha, w%hv, w%hvn, error)
      if (allocated(error)) return
      call keep_tangent(run, w)
      call dpcm_original_point_tangent(st%dpcm_original_point_state, inp%zr, inp%zr, inp%zs, &
         & inp%zs, inp%zs, inp%dv, inp%dvn, w%hr, w%hn, w%hk1, w%hk2, w%ha, w%hv, w%hvn, error)
      if (allocated(error)) return
      run%response_tangent = stack(w%hv, w%hvn)

      run%cx_param = eps
      call dpcm_original_point_complex_medium(st, run%cx_param, error)
      if (allocated(error)) return
      call dpcm_original_point_complex_solve(st, cmplx(inp%v, 0.0_wp, wp), &
         & cmplx(inp%vn, 0.0_wp, wp), ez, z%q, z%gv, z%gvn, error)
      if (allocated(error)) return
      call keep_complex_solve(run, ez, z)
      call dpcm_original_point_complex_response(st, cmplx(inp%dv, 0.0_wp, wp), &
         & cmplx(inp%dvn, 0.0_wp, wp), z%hv, z%hvn, error)
      if (allocated(error)) return
      run%cx_response = stack_complex(z%hv, z%hvn)
      call dpcm_original_point_complex_param_weights(st, run%cx_params(1), error)
      if (allocated(error)) return
      call dpcm_original_point_complex_weights(st, z%gr, z%gn, z%gk1, z%gk2, z%ga, error)
      if (allocated(error)) return
      call dpcm_original_point_complex_tangent(st, inp%dr, inp%dn, inp%dk1, inp%dk2, inp%da, &
         & cmplx(inp%dv, 0.0_wp, wp), cmplx(inp%dvn, 0.0_wp, wp), z%hr, z%hn, z%hk1, z%hk2, z%ha, &
         & z%hv, z%hvn, error)
      if (allocated(error)) return
      call keep_complex_derivatives(run, z)
      call dpcm_original_point_complex_tangent(st, inp%zr, inp%zr, inp%zs, inp%zs, inp%zs, &
         & cmplx(inp%dv, 0.0_wp, wp), cmplx(inp%dvn, 0.0_wp, wp), z%hr, z%hn, z%hk1, z%hk2, z%ha, &
         & z%hv, z%hvn, error)
      if (allocated(error)) return
      run%cx_response_tangent = stack_complex(z%hv, z%hvn)

      call dpcm_original_point_solve(st%dpcm_original_point_state, inp%v, inp%vn, e, w%q, w%gv, &
         & w%gvn, error)
      if (allocated(error)) return
      call keep_again(run, e, w)
      do k = 1, size(fd6_offsets)
         call dpcm_original_point_solve(st%dpcm_original_point_state, &
            & inp%v + fd6_offsets(k)*inp%h*inp%dv(:, 1), &
            & inp%vn + fd6_offsets(k)*inp%h*inp%dvn(:, 1), e, w%q, w%gv, w%gvn, error)
         if (allocated(error)) return
         run%stencil(:, k) = [w%gv, w%gvn]
      end do
      do k = 1, size(fd6_offsets)
         call dpcm_original_point_complex_solve(st, &
            & cmplx(inp%v + fd6_offsets(k)*inp%h*inp%dv(:, 1), 0.0_wp, wp), &
            & cmplx(inp%vn + fd6_offsets(k)*inp%h*inp%dvn(:, 1), 0.0_wp, wp), ez, z%q, z%gv, &
            & z%gvn, error)
         if (allocated(error)) return
         run%cx_stencil(:, k) = [z%gv, z%gvn]
      end do
   end subroutine run_dpcm_original_point

   !> LPB, Gaussian kernel: real and complex phases and stencils
   !>
   !> @param[in]     inp    Kernel inputs
   !> @param[in]     solver Linear solver
   !> @param[in]     eps    Permittivity of the complex phases
   !> @param[in,out] run    Kernel results
   !> @param[out]    error  Library error
   subroutine run_lpb_gauss(inp, solver, eps, run, error)
      !> Kernel inputs
      type(kernel_inputs), intent(in) :: inp
      !> Linear solver
      class(pcm_solver), intent(in) :: solver
      !> Permittivity of the complex phases
      complex(wp), intent(in) :: eps
      !> Kernel results
      type(kernel_run), intent(inout) :: run
      !> Library error
      type(moist_error_type), allocatable, intent(out) :: error

      !> Kernel state of both variants
      type(lpb_complex_state), target :: st
      !> Output buffers
      type(real_work) :: w
      type(complex_work) :: z
      !> Energies
      real(wp) :: e
      complex(wp) :: ez
      !> Stencil index
      integer :: k

      call new_real_work(inp%npt, w)
      call new_complex_work(inp%npt, z)
      run%has_complex = .true.
      allocate (run%params(2), run%cx_params(2))
      call lpb_assemble(st%lpb_state, inp%r, inp%n, inp%k1, inp%k2, inp%xi, inp%f, inp%a, error, &
         & solver)
      if (allocated(error)) return
      call lpb_medium(st%lpb_state, inp%phys%eps, inp%phys%kappa, error)
      if (allocated(error)) return
      call lpb_solve(st%lpb_state, inp%v, e, w%q, w%gv, error)
      if (allocated(error)) return
      call keep_solve(run, e, w)
      call lpb_response(st%lpb_state, inp%dv, w%hv, error)
      if (allocated(error)) return
      run%response = stack(w%hv, w%hvn)
      call lpb_param_weights(st%lpb_state, run%params(1), run%params(2), error)
      if (allocated(error)) return
      call lpb_weights(st%lpb_state, w%gr, w%gn, w%gk1, w%gk2, w%gxi, w%gf, w%ga, error)
      if (allocated(error)) return
      call keep_weights(run, w)
      call lpb_tangent(st%lpb_state, inp%dr, inp%dn, inp%dk1, inp%dk2, inp%dxi, inp%df, inp%da, &
         & inp%dv, w%hr, w%hn, w%hk1, w%hk2, w%hxi, w%hf, w%ha, w%hv, error)
      if (allocated(error)) return
      call keep_tangent(run, w)
      call lpb_tangent(st%lpb_state, inp%zr, inp%zr, inp%zs, inp%zs, inp%zs, inp%zs, inp%zs, &
         & inp%dv, w%hr, w%hn, w%hk1, w%hk2, w%hxi, w%hf, w%ha, w%hv, error)
      if (allocated(error)) return
      run%response_tangent = stack(w%hv, w%hvn)

      run%cx_param = eps
      call lpb_complex_medium(st, run%cx_param, inp%phys%kappa, error)
      if (allocated(error)) return
      call lpb_complex_solve(st, cmplx(inp%v, 0.0_wp, wp), ez, z%q, z%gv, error)
      if (allocated(error)) return
      call keep_complex_solve(run, ez, z)
      call lpb_complex_response(st, cmplx(inp%dv, 0.0_wp, wp), z%hv, error)
      if (allocated(error)) return
      run%cx_response = stack_complex(z%hv, z%hvn)
      call lpb_complex_param_weights(st, run%cx_params(1), run%cx_params(2), error)
      if (allocated(error)) return
      call lpb_complex_weights(st, z%gr, z%gn, z%gk1, z%gk2, z%gxi, z%gf, z%ga, error)
      if (allocated(error)) return
      call lpb_complex_tangent(st, inp%dr, inp%dn, inp%dk1, inp%dk2, inp%dxi, inp%df, inp%da, &
         & cmplx(inp%dv, 0.0_wp, wp), z%hr, z%hn, z%hk1, z%hk2, z%hxi, z%hf, z%ha, z%hv, error)
      if (allocated(error)) return
      call keep_complex_derivatives(run, z)
      call lpb_complex_tangent(st, inp%zr, inp%zr, inp%zs, inp%zs, inp%zs, inp%zs, inp%zs, &
         & cmplx(inp%dv, 0.0_wp, wp), z%hr, z%hn, z%hk1, z%hk2, z%hxi, z%hf, z%ha, z%hv, error)
      if (allocated(error)) return
      run%cx_response_tangent = stack_complex(z%hv, z%hvn)

      call lpb_solve(st%lpb_state, inp%v, e, w%q, w%gv, error)
      if (allocated(error)) return
      call keep_again(run, e, w)
      do k = 1, size(fd6_offsets)
         call lpb_solve(st%lpb_state, inp%v + fd6_offsets(k)*inp%h*inp%dv(:, 1), e, w%q, w%gv, &
            & error)
         if (allocated(error)) return
         run%stencil(:, k) = [w%gv, w%gvn]
      end do
      do k = 1, size(fd6_offsets)
         call lpb_complex_solve(st, cmplx(inp%v + fd6_offsets(k)*inp%h*inp%dv(:, 1), 0.0_wp, wp), &
            & ez, z%q, z%gv, error)
         if (allocated(error)) return
         run%cx_stencil(:, k) = [z%gv, z%gvn]
      end do
   end subroutine run_lpb_gauss

   !> LPB, point kernel: real and complex phases and stencils
   !>
   !> @param[in]     inp    Kernel inputs
   !> @param[in]     solver Linear solver
   !> @param[in]     eps    Permittivity of the complex phases
   !> @param[in,out] run    Kernel results
   !> @param[out]    error  Library error
   subroutine run_lpb_point(inp, solver, eps, run, error)
      !> Kernel inputs
      type(kernel_inputs), intent(in) :: inp
      !> Linear solver
      class(pcm_solver), intent(in) :: solver
      !> Permittivity of the complex phases
      complex(wp), intent(in) :: eps
      !> Kernel results
      type(kernel_run), intent(inout) :: run
      !> Library error
      type(moist_error_type), allocatable, intent(out) :: error

      !> Kernel state of both variants
      type(lpb_point_complex_state), target :: st
      !> Output buffers
      type(real_work) :: w
      type(complex_work) :: z
      !> Energies
      real(wp) :: e
      complex(wp) :: ez
      !> Stencil index
      integer :: k

      call new_real_work(inp%npt, w)
      call new_complex_work(inp%npt, z)
      run%has_complex = .true.
      allocate (run%params(2), run%cx_params(2))
      call lpb_point_assemble(st%lpb_point_state, inp%r, inp%n, inp%k1, inp%k2, inp%a, error, &
         & solver)
      if (allocated(error)) return
      call lpb_point_medium(st%lpb_point_state, inp%phys%eps, inp%phys%kappa, error)
      if (allocated(error)) return
      call lpb_point_solve(st%lpb_point_state, inp%v, e, w%q, w%gv, error)
      if (allocated(error)) return
      call keep_solve(run, e, w)
      call lpb_point_response(st%lpb_point_state, inp%dv, w%hv, error)
      if (allocated(error)) return
      run%response = stack(w%hv, w%hvn)
      call lpb_point_param_weights(st%lpb_point_state, run%params(1), run%params(2), error)
      if (allocated(error)) return
      call lpb_point_weights(st%lpb_point_state, w%gr, w%gn, w%gk1, w%gk2, w%ga, error)
      if (allocated(error)) return
      call keep_weights(run, w)
      call lpb_point_tangent(st%lpb_point_state, inp%dr, inp%dn, inp%dk1, inp%dk2, inp%da, inp%dv, &
         & w%hr, w%hn, w%hk1, w%hk2, w%ha, w%hv, error)
      if (allocated(error)) return
      call keep_tangent(run, w)
      call lpb_point_tangent(st%lpb_point_state, inp%zr, inp%zr, inp%zs, inp%zs, inp%zs, inp%dv, &
         & w%hr, w%hn, w%hk1, w%hk2, w%ha, w%hv, error)
      if (allocated(error)) return
      run%response_tangent = stack(w%hv, w%hvn)

      run%cx_param = eps
      call lpb_point_complex_medium(st, run%cx_param, inp%phys%kappa, error)
      if (allocated(error)) return
      call lpb_point_complex_solve(st, cmplx(inp%v, 0.0_wp, wp), ez, z%q, z%gv, error)
      if (allocated(error)) return
      call keep_complex_solve(run, ez, z)
      call lpb_point_complex_response(st, cmplx(inp%dv, 0.0_wp, wp), z%hv, error)
      if (allocated(error)) return
      run%cx_response = stack_complex(z%hv, z%hvn)
      call lpb_point_complex_param_weights(st, run%cx_params(1), run%cx_params(2), error)
      if (allocated(error)) return
      call lpb_point_complex_weights(st, z%gr, z%gn, z%gk1, z%gk2, z%ga, error)
      if (allocated(error)) return
      call lpb_point_complex_tangent(st, inp%dr, inp%dn, inp%dk1, inp%dk2, inp%da, &
         & cmplx(inp%dv, 0.0_wp, wp), z%hr, z%hn, z%hk1, z%hk2, z%ha, z%hv, error)
      if (allocated(error)) return
      call keep_complex_derivatives(run, z)
      call lpb_point_complex_tangent(st, inp%zr, inp%zr, inp%zs, inp%zs, inp%zs, &
         & cmplx(inp%dv, 0.0_wp, wp), z%hr, z%hn, z%hk1, z%hk2, z%ha, z%hv, error)
      if (allocated(error)) return
      run%cx_response_tangent = stack_complex(z%hv, z%hvn)

      call lpb_point_solve(st%lpb_point_state, inp%v, e, w%q, w%gv, error)
      if (allocated(error)) return
      call keep_again(run, e, w)
      do k = 1, size(fd6_offsets)
         call lpb_point_solve(st%lpb_point_state, inp%v + fd6_offsets(k)*inp%h*inp%dv(:, 1), e, &
            & w%q, w%gv, error)
         if (allocated(error)) return
         run%stencil(:, k) = [w%gv, w%gvn]
      end do
      do k = 1, size(fd6_offsets)
         call lpb_point_complex_solve(st, &
            & cmplx(inp%v + fd6_offsets(k)*inp%h*inp%dv(:, 1), 0.0_wp, wp), ez, z%q, z%gv, error)
         if (allocated(error)) return
         run%cx_stencil(:, k) = [z%gv, z%gvn]
      end do
   end subroutine run_lpb_point

   !> SLIC, Gaussian kernel: real phases and stencil
   !>
   !> @param[in]     inp    Kernel inputs
   !> @param[in]     solver Linear solver
   !> @param[in]     eps    Permittivity of the complex phases, unused (no complex variant)
   !> @param[in,out] run    Kernel results
   !> @param[out]    error  Library error
   subroutine run_slic_gauss(inp, solver, eps, run, error)
      !> Kernel inputs
      type(kernel_inputs), intent(in) :: inp
      !> Linear solver
      class(pcm_solver), intent(in) :: solver
      !> Permittivity of the complex phases
      complex(wp), intent(in) :: eps
      !> Kernel results
      type(kernel_run), intent(inout) :: run
      !> Library error
      type(moist_error_type), allocatable, intent(out) :: error

      !> Kernel state
      type(slic_state), target :: st
      !> Output buffers
      type(real_work) :: w
      !> Energy
      real(wp) :: e
      !> Stencil index
      integer :: k

      associate (unused_eps => eps)
      end associate
      call new_real_work(inp%npt, w)
      allocate (run%params(5))
      call slic_assemble(st, inp%r, inp%n, inp%k1, inp%k2, inp%xi, inp%a, error, solver)
      if (allocated(error)) return
      call slic_medium(st, inp%phys%eps, inp%phys%slic_alpha, inp%phys%slic_beta, &
         & inp%phys%slic_gamma, inp%phys%slic_mu, error)
      if (allocated(error)) return
      call slic_solve(st, inp%v, inp%vn, e, w%q, w%gv, w%gvn, error)
      if (allocated(error)) return
      call keep_solve(run, e, w)
      call slic_response(st, inp%dv, inp%dvn, w%hv, w%hvn, error)
      if (allocated(error)) return
      run%response = stack(w%hv, w%hvn)
      call slic_param_weights(st, run%params(1), run%params(2), run%params(3), run%params(4), &
         & run%params(5), error)
      if (allocated(error)) return
      call slic_weights(st, w%gr, w%gn, w%gk1, w%gk2, w%gxi, w%ga, error)
      if (allocated(error)) return
      call keep_weights(run, w)
      call slic_tangent(st, inp%dr, inp%dn, inp%dk1, inp%dk2, inp%dxi, inp%da, inp%dv, inp%dvn, &
         & w%hr, w%hn, w%hk1, w%hk2, w%hxi, w%ha, w%hv, w%hvn, error)
      if (allocated(error)) return
      call keep_tangent(run, w)
      call slic_tangent(st, inp%zr, inp%zr, inp%zs, inp%zs, inp%zs, inp%zs, inp%dv, inp%dvn, w%hr, &
         & w%hn, w%hk1, w%hk2, w%hxi, w%ha, w%hv, w%hvn, error)
      if (allocated(error)) return
      run%response_tangent = stack(w%hv, w%hvn)

      call slic_solve(st, inp%v, inp%vn, e, w%q, w%gv, w%gvn, error)
      if (allocated(error)) return
      call keep_again(run, e, w)
      do k = 1, size(fd6_offsets)
         call slic_solve(st, inp%v + fd6_offsets(k)*inp%h*inp%dv(:, 1), &
            & inp%vn + fd6_offsets(k)*inp%h*inp%dvn(:, 1), e, w%q, w%gv, w%gvn, error)
         if (allocated(error)) return
         run%stencil(:, k) = [w%gv, w%gvn]
      end do
   end subroutine run_slic_gauss

   !> SLIC, point kernel: real phases and stencil
   !>
   !> @param[in]     inp    Kernel inputs
   !> @param[in]     solver Linear solver
   !> @param[in]     eps    Permittivity of the complex phases, unused (no complex variant)
   !> @param[in,out] run    Kernel results
   !> @param[out]    error  Library error
   subroutine run_slic_point(inp, solver, eps, run, error)
      !> Kernel inputs
      type(kernel_inputs), intent(in) :: inp
      !> Linear solver
      class(pcm_solver), intent(in) :: solver
      !> Permittivity of the complex phases
      complex(wp), intent(in) :: eps
      !> Kernel results
      type(kernel_run), intent(inout) :: run
      !> Library error
      type(moist_error_type), allocatable, intent(out) :: error

      !> Kernel state
      type(slic_point_state), target :: st
      !> Output buffers
      type(real_work) :: w
      !> Energy
      real(wp) :: e
      !> Stencil index
      integer :: k

      associate (unused_eps => eps)
      end associate
      call new_real_work(inp%npt, w)
      allocate (run%params(5))
      call slic_point_assemble(st, inp%r, inp%n, inp%k1, inp%k2, inp%a, error, solver)
      if (allocated(error)) return
      call slic_point_medium(st, inp%phys%eps, inp%phys%slic_alpha, inp%phys%slic_beta, &
         & inp%phys%slic_gamma, inp%phys%slic_mu, error)
      if (allocated(error)) return
      call slic_point_solve(st, inp%v, inp%vn, e, w%q, w%gv, w%gvn, error)
      if (allocated(error)) return
      call keep_solve(run, e, w)
      call slic_point_response(st, inp%dv, inp%dvn, w%hv, w%hvn, error)
      if (allocated(error)) return
      run%response = stack(w%hv, w%hvn)
      call slic_point_param_weights(st, run%params(1), run%params(2), run%params(3), &
         & run%params(4), run%params(5), error)
      if (allocated(error)) return
      call slic_point_weights(st, w%gr, w%gn, w%gk1, w%gk2, w%ga, error)
      if (allocated(error)) return
      call keep_weights(run, w)
      call slic_point_tangent(st, inp%dr, inp%dn, inp%dk1, inp%dk2, inp%da, inp%dv, inp%dvn, w%hr, &
         & w%hn, w%hk1, w%hk2, w%ha, w%hv, w%hvn, error)
      if (allocated(error)) return
      call keep_tangent(run, w)
      call slic_point_tangent(st, inp%zr, inp%zr, inp%zs, inp%zs, inp%zs, inp%dv, inp%dvn, w%hr, &
         & w%hn, w%hk1, w%hk2, w%ha, w%hv, w%hvn, error)
      if (allocated(error)) return
      run%response_tangent = stack(w%hv, w%hvn)

      call slic_point_solve(st, inp%v, inp%vn, e, w%q, w%gv, w%gvn, error)
      if (allocated(error)) return
      call keep_again(run, e, w)
      do k = 1, size(fd6_offsets)
         call slic_point_solve(st, inp%v + fd6_offsets(k)*inp%h*inp%dv(:, 1), &
            & inp%vn + fd6_offsets(k)*inp%h*inp%dvn(:, 1), e, w%q, w%gv, w%gvn, error)
         if (allocated(error)) return
         run%stencil(:, k) = [w%gv, w%gvn]
      end do
   end subroutine run_slic_point

   !> Nonlocal, Gaussian kernel: real phases and stencil
   !>
   !> @param[in]     inp    Kernel inputs
   !> @param[in]     solver Linear solver
   !> @param[in]     eps    Permittivity of the complex phases, unused (no complex variant)
   !> @param[in,out] run    Kernel results
   !> @param[out]    error  Library error
   subroutine run_nonlocal_gauss(inp, solver, eps, run, error)
      !> Kernel inputs
      type(kernel_inputs), intent(in) :: inp
      !> Linear solver
      class(pcm_solver), intent(in) :: solver
      !> Permittivity of the complex phases
      complex(wp), intent(in) :: eps
      !> Kernel results
      type(kernel_run), intent(inout) :: run
      !> Library error
      type(moist_error_type), allocatable, intent(out) :: error

      !> Kernel state
      type(nonlocal_pcm_state), target :: st
      !> Output buffers
      type(real_work) :: w
      !> Energy
      real(wp) :: e
      !> Stencil index
      integer :: k

      associate (unused_eps => eps)
      end associate
      call new_real_work(inp%npt, w)
      allocate (run%params(4))
      call nonlocal_pcm_assemble(st, inp%r, inp%n, inp%k1, inp%k2, inp%xi, inp%f, inp%a, error, &
         & solver)
      if (allocated(error)) return
      call nonlocal_pcm_medium(st, inp%phys%eps, inp%phys%eps_inf, inp%phys%lambda, &
         & inp%phys%eps_solute, error)
      if (allocated(error)) return
      call nonlocal_pcm_solve(st, inp%v, inp%vn, e, w%u, w%q, w%wext, w%gv, w%gvn, error)
      if (allocated(error)) return
      call keep_solve(run, e, w)
      call nonlocal_pcm_response(st, inp%dv, inp%dvn, w%hv, w%hvn, error)
      if (allocated(error)) return
      run%response = stack(w%hv, w%hvn)
      call nonlocal_pcm_param_weights(st, run%params(1), run%params(2), run%params(3), &
         & run%params(4), error)
      if (allocated(error)) return
      call nonlocal_pcm_weights(st, w%gr, w%gn, w%gk1, w%gk2, w%gxi, w%gf, w%ga, error)
      if (allocated(error)) return
      call keep_weights(run, w)
      call nonlocal_pcm_tangent(st, inp%dr, inp%dn, inp%dk1, inp%dk2, inp%dxi, inp%df, inp%da, &
         & inp%dv, inp%dvn, w%hr, w%hn, w%hk1, w%hk2, w%hxi, w%hf, w%ha, w%hv, w%hvn, error)
      if (allocated(error)) return
      call keep_tangent(run, w)
      call nonlocal_pcm_tangent(st, inp%zr, inp%zr, inp%zs, inp%zs, inp%zs, inp%zs, inp%zs, &
         & inp%dv, inp%dvn, w%hr, w%hn, w%hk1, w%hk2, w%hxi, w%hf, w%ha, w%hv, w%hvn, error)
      if (allocated(error)) return
      run%response_tangent = stack(w%hv, w%hvn)

      call nonlocal_pcm_solve(st, inp%v, inp%vn, e, w%u, w%q, w%wext, w%gv, w%gvn, error)
      if (allocated(error)) return
      call keep_again(run, e, w)
      do k = 1, size(fd6_offsets)
         call nonlocal_pcm_solve(st, inp%v + fd6_offsets(k)*inp%h*inp%dv(:, 1), &
            & inp%vn + fd6_offsets(k)*inp%h*inp%dvn(:, 1), e, w%u, w%q, w%wext, w%gv, w%gvn, &
            & error)
         if (allocated(error)) return
         run%stencil(:, k) = [w%gv, w%gvn]
      end do
   end subroutine run_nonlocal_gauss

   !> Nonlocal, point kernel: real phases and stencil
   !>
   !> @param[in]     inp    Kernel inputs
   !> @param[in]     solver Linear solver
   !> @param[in]     eps    Permittivity of the complex phases, unused (no complex variant)
   !> @param[in,out] run    Kernel results
   !> @param[out]    error  Library error
   subroutine run_nonlocal_point(inp, solver, eps, run, error)
      !> Kernel inputs
      type(kernel_inputs), intent(in) :: inp
      !> Linear solver
      class(pcm_solver), intent(in) :: solver
      !> Permittivity of the complex phases
      complex(wp), intent(in) :: eps
      !> Kernel results
      type(kernel_run), intent(inout) :: run
      !> Library error
      type(moist_error_type), allocatable, intent(out) :: error

      !> Kernel state
      type(nonlocal_pcm_point_state), target :: st
      !> Output buffers
      type(real_work) :: w
      !> Energy
      real(wp) :: e
      !> Stencil index
      integer :: k

      associate (unused_eps => eps)
      end associate
      call new_real_work(inp%npt, w)
      allocate (run%params(4))
      call nonlocal_pcm_point_assemble(st, inp%r, inp%n, inp%k1, inp%k2, inp%a, error, solver)
      if (allocated(error)) return
      call nonlocal_pcm_point_medium(st, inp%phys%eps, inp%phys%eps_inf, inp%phys%lambda, &
         & inp%phys%eps_solute, error)
      if (allocated(error)) return
      call nonlocal_pcm_point_solve(st, inp%v, inp%vn, e, w%u, w%q, w%wext, w%gv, w%gvn, error)
      if (allocated(error)) return
      call keep_solve(run, e, w)
      call nonlocal_pcm_point_response(st, inp%dv, inp%dvn, w%hv, w%hvn, error)
      if (allocated(error)) return
      run%response = stack(w%hv, w%hvn)
      call nonlocal_pcm_point_param_weights(st, run%params(1), run%params(2), run%params(3), &
         & run%params(4), error)
      if (allocated(error)) return
      call nonlocal_pcm_point_weights(st, w%gr, w%gn, w%gk1, w%gk2, w%ga, error)
      if (allocated(error)) return
      call keep_weights(run, w)
      call nonlocal_pcm_point_tangent(st, inp%dr, inp%dn, inp%dk1, inp%dk2, inp%da, inp%dv, &
         & inp%dvn, w%hr, w%hn, w%hk1, w%hk2, w%ha, w%hv, w%hvn, error)
      if (allocated(error)) return
      call keep_tangent(run, w)
      call nonlocal_pcm_point_tangent(st, inp%zr, inp%zr, inp%zs, inp%zs, inp%zs, inp%dv, inp%dvn, &
         & w%hr, w%hn, w%hk1, w%hk2, w%ha, w%hv, w%hvn, error)
      if (allocated(error)) return
      run%response_tangent = stack(w%hv, w%hvn)

      call nonlocal_pcm_point_solve(st, inp%v, inp%vn, e, w%u, w%q, w%wext, w%gv, w%gvn, error)
      if (allocated(error)) return
      call keep_again(run, e, w)
      do k = 1, size(fd6_offsets)
         call nonlocal_pcm_point_solve(st, inp%v + fd6_offsets(k)*inp%h*inp%dv(:, 1), &
            & inp%vn + fd6_offsets(k)*inp%h*inp%dvn(:, 1), e, w%u, w%q, w%wext, w%gv, w%gvn, &
            & error)
         if (allocated(error)) return
         run%stencil(:, k) = [w%gv, w%gvn]
      end do
   end subroutine run_nonlocal_point

   !> Dielectric factor (eps - 1)/eps of a conductor-like model
   !>
   !> @param[in] eps Permittivity
   pure function feps_of(eps) result(feps)
      !> Permittivity
      real(wp), intent(in) :: eps
      !> Dielectric factor
      real(wp) :: feps

      feps = (eps - 1.0_wp)/eps
   end function feps_of

   !> Dielectric factor (eps - 1)/eps of a conductor-like model, complex
   !>
   !> @param[in] eps Permittivity
   pure function cfeps_of(eps) result(feps)
      !> Permittivity
      complex(wp), intent(in) :: eps
      !> Dielectric factor
      complex(wp) :: feps

      feps = (eps - 1.0_wp)/eps
   end function cfeps_of

end module test_model_component_pcm_kernel_variants
