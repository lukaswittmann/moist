!> Adjoint-level unit tests for the Gaussian PCM interaction matrix
!>
!> [[moist_model_component_pcm_amat_adjoint]] contracts the derivative of the
!> matrix into O(ngrid) surface-variable weights and then folds those weights
!> against the cavity's nuclear derivative arrays
!>
!> Comparisons:
!>  * the surface-variable weights against 4-point central differences of the
!>    contracted energy `q1^T A q2` on a molecular surface, where a grid-point has
!>    hundreds of partners spread across both kernel branches
!>  * the nuclear gradient against a contraction of the dense derivative tensor
!>    [[assemble_pcm_amat_with_gradient]] builds entry by entry, which is a
!>    second, independent route to the same number
module test_model_component_pcm_amat_adjoint
   use mctc_env, only: wp
   use mctc_env_error, only: moist_error_type => error_type
   use mctc_io, only: structure_type
   use moist_cavity_iswig, only: cavity_type_iswig
   use moist_model_component_pcm_amat, only: assemble_pcm_amat, &
                                             assemble_pcm_amat_with_gradient, &
                                             pcm_amat_surface_weights, &
                                             pcm_amat_nuclear_gradient
   use test_helpers, only: get_test_structures, center_at_origin, fd4_scalar, &
                           fd4_offsets, get_test_cavity_iswig, rel_deviation
   use test_model_component_pcm_amat, only: count_branches, nmol, nleb_survey
   use testdrive, only: new_unittest, unittest_type, error_type, check, test_failed
   implicit none(type, external)
   private

   public :: collect_model_component_pcm_amat_adjoint

   !> Finite-difference steps for the surface-variable channels
   real(wp), parameter :: rel_step = 3.0e-3_wp
   real(wp), parameter :: xyz_step = 3.0e-3_wp

   !> 4-point central FD tolerances for the surface-weight comparison
   real(wp), parameter :: fd_atol = 1.0e-10_wp
   real(wp), parameter :: fd_rtol = 1.0e-9_wp

   !> Smallest switching factor used for finite-difference, applied as the
   !> cavity's own cutoff so only the exposed surface is ever built
   real(wp), parameter :: fd_min_f = 0.1_wp

   !> Smallest exposed surface the finite-difference tests accept
   integer, parameter :: fd_min_points = 50

   !> Number of  grid points finite-differenced per structure
   integer, parameter :: n_fd_points = 50

contains

   !> Collect the adjoint test suite
   subroutine collect_model_component_pcm_amat_adjoint(testsuite)
      !> Collected unit tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("surface_weights_vs_fd", test_surface_weights_vs_fd), &
                  new_unittest("weight_paths_agree", test_weight_paths_agree) &
                  ]

   end subroutine collect_model_component_pcm_amat_adjoint

   !* --------------------------------- Local helpers --------------------------------- *!

   !> Deterministic contraction vectors for the adjoint tests
   !>
   !> Two different, sign-alternating profiles so that q1 /= q2 and the
   !> symmetrized charge product q1_i q2_j + q1_j q2_i is not simply q1_i q1_j;
   !> an accidental symmetrization in the weights would otherwise pass
   !>
   !> @param[in]  ngrid  Number of  grid points
   !> @param[out] q1     Left contraction vector
   !> @param[out] q2     Right contraction vector
   subroutine make_charges(ngrid, q1, q2)
      !> Number of  grid points
      integer, intent(in) :: ngrid
      !> Left and right contraction vectors
      real(wp), allocatable, intent(out) :: q1(:), q2(:)

      !> Grid-point index
      integer :: i

      allocate (q1(ngrid), q2(ngrid))
      do i = 1, ngrid
         q1(i) = 0.3_wp*sin(0.7_wp*real(i, wp)) + 0.1_wp
         q2(i) = 0.2_wp*cos(0.31_wp*real(i, wp)) - 0.05_wp
      end do

   end subroutine make_charges

   !> Contracted energy q1^T A q2 on a given surface
   !>
   !> @param[in]  xi     Gaussian widths (ngrid)
   !> @param[in]  f      Gaussian switching factors (ngrid)
   !> @param[in]  xyz    Surface positions (3, ngrid)
   !> @param[in]  q1     Left contraction vector
   !> @param[in]  q2     Right contraction vector
   !> @param[out] energy Contracted energy
   !> @param[out] err    Error handling
   subroutine contracted_energy(xi, f, xyz, q1, q2, energy, err)
      !> Gaussian widths, switching factors and surface positions
      real(wp), intent(in) :: xi(:), f(:), xyz(:, :)
      !> Left and right contraction vectors
      real(wp), intent(in) :: q1(:), q2(:)
      !> Contracted energy
      real(wp), intent(out) :: energy
      !> Error handling
      type(moist_error_type), allocatable, intent(out) :: err

      !> Assembled matrix and its product with the right vector
      real(wp), allocatable :: amat(:, :), aq(:)
      !> Grid size
      integer :: ngrid

      ngrid = size(xi)
      allocate (amat(ngrid, ngrid), aq(ngrid))
      call assemble_pcm_amat(xi, f, xyz, amat, err)
      energy = 0.0_wp
      if (allocated(err)) return
      aq = matmul(amat, q2)
      energy = dot_product(q1, aq)

   end subroutine contracted_energy

   !* ------------------------------------- Tests ------------------------------------- *!

   !> Surface-variable adjoints against 4-point central differences
   subroutine test_surface_weights_vs_fd(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      !> Sampled structures
      type(structure_type), allocatable :: mols(:)
      !> Surface under test
      type(cavity_type_iswig) :: cavity
      !> Library error handling
      type(moist_error_type), allocatable :: err
      !> Working copies of the surface variables
      real(wp), allocatable :: xi(:), f(:), xyz(:, :)
      !> Contraction vectors
      real(wp), allocatable :: q1(:), q2(:)
      !> Analytic adjoints
      real(wp), allocatable :: w_xi(:), w_f(:), w_xyz(:, :)
      !> Grid-point, axis, stencil and stride indices; grid size
      integer :: ip, ig, iax, k, stride, ngrid
      !> Per-branch pair counts of the differenced surface
      integer :: nfar, nnear
      !> Stencil values, step size, finite difference and saved coordinate
      real(wp) :: vals(4), step, fd, saved
      !> Failure context
      character(len=128) :: context

      call get_test_structures(mols, nmol)
      call center_at_origin(mols(1))
      ! `cut_f` drops the buried  grid points at construction: their near-zero
      ! switching factor makes both the contracted energy and the relative FD
      ! step degenerate. The working copies below are what the stencils poke.
      call get_test_cavity_iswig(mols(1), cavity, err, nleb=nleb_survey, &
                                 cut_f=fd_min_f)
      if (allocated(err)) then
         call test_failed(error, "surface construction failed: "//err%message)
         return
      end if

      ngrid = cavity%ngrid
      xi = cavity%xi0
      f = cavity%f
      xyz = cavity%xyz
      call check(error, ngrid >= fd_min_points, &
                 more="exposed surface too small to difference meaningfully")
      if (allocated(error)) return

      call count_branches(xi, xyz, nfar, nnear)
      call check(error, nnear > 0 .and. nfar > 0, &
                 more="exposed surface does not span both kernel branches")
      if (allocated(error)) return

      call make_charges(ngrid, q1, q2)

      allocate (w_xi(ngrid), w_f(ngrid), w_xyz(3, ngrid))
      call pcm_amat_surface_weights(xi, f, xyz, q1, q2, w_xi, w_f, w_xyz, err)
      if (allocated(err)) then
         call test_failed(error, "surface-weight assembly failed: "//err%message)
         return
      end if

      stride = max(1, ngrid/n_fd_points)
      do ip = 1, n_fd_points
         ig = 1 + (ip - 1)*stride
         if (ig > ngrid) exit

         ! Gaussian width channel.
         step = rel_step*xi(ig)
         saved = xi(ig)
         do k = 1, 4
            xi(ig) = saved + fd4_offsets(k)*step
            call contracted_energy(xi, f, xyz, q1, q2, vals(k), err)
            if (allocated(err)) then
               call test_failed(error, "perturbed assembly failed: "//err%message)
               return
            end if
         end do
         xi(ig) = saved
         fd = fd4_scalar(vals(1), vals(2), vals(3), vals(4), step)
         write (context, "(a,i0)") "dE/dxi at grid-point ", ig
         call check(error, w_xi(ig), fd, thr=fd_atol + fd_rtol*abs(fd), &
                    more=trim(context))
         if (allocated(error)) return

         ! Switching-factor channel. It enters only through the diagonal, so a
         ! sign or prefactor slip there shows up here and nowhere else.
         step = rel_step*f(ig)
         saved = f(ig)
         do k = 1, 4
            f(ig) = saved + fd4_offsets(k)*step
            call contracted_energy(xi, f, xyz, q1, q2, vals(k), err)
            if (allocated(err)) then
               call test_failed(error, "perturbed assembly failed: "//err%message)
               return
            end if
         end do
         f(ig) = saved
         fd = fd4_scalar(vals(1), vals(2), vals(3), vals(4), step)
         write (context, "(a,i0)") "dE/df at grid-point ", ig
         call check(error, w_f(ig), fd, thr=fd_atol + fd_rtol*abs(fd), &
                    more=trim(context))
         if (allocated(error)) return

         ! Position channel.
         do iax = 1, 3
            step = xyz_step
            saved = xyz(iax, ig)
            do k = 1, 4
               xyz(iax, ig) = saved + fd4_offsets(k)*step
               call contracted_energy(xi, f, xyz, q1, q2, vals(k), err)
               if (allocated(err)) then
                  call test_failed(error, "perturbed assembly failed: "//err%message)
                  return
               end if
            end do
            xyz(iax, ig) = saved
            fd = fd4_scalar(vals(1), vals(2), vals(3), vals(4), step)
            write (context, "(a,i0,a,i0)") "dE/dxyz axis ", iax, " at grid-point ", ig
            call check(error, w_xyz(iax, ig), fd, thr=fd_atol + fd_rtol*abs(fd), &
                       more=trim(context))
            if (allocated(error)) return
         end do
      end do

   end subroutine test_surface_weights_vs_fd

   !> The dense derivative tensor and the O(ngrid) adjoint agree
   subroutine test_weight_paths_agree(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      !> Sampled structures
      type(structure_type), allocatable :: mols(:)
      !> Surface under test
      type(cavity_type_iswig) :: cavity
      !> Library error handling
      type(moist_error_type), allocatable :: err
      !> Plain and dense-path matrices
      real(wp), allocatable :: amat(:, :), amat_dense(:, :)
      !> Dense nuclear derivative tensor
      real(wp), allocatable :: amat1_rA(:, :, :, :)
      !> Contraction vectors
      real(wp), allocatable :: q1(:), q2(:)
      !> Analytic adjoints
      real(wp), allocatable :: w_xi(:), w_f(:), w_xyz(:, :)
      !> Gradients from the adjoint and the dense route
      real(wp), allocatable :: grad(:, :), grad_ref(:, :)
      !> Matrix, atom, axis indices and extents
      integer :: i, j, iatom, iaxis, ngrid, nsph
      !> Largest deviations of the matrix and of the gradient
      real(wp) :: dev_amat, dev_grad, scale

      !> The dense tensor is (3, nsph, ngrid, ngrid), so this test only runs on
      !> the coarsest iSwiG grid; nleb 14 keeps it in the tens of megabytes.
      integer, parameter :: nleb_dense_tensor = 14

      call get_test_structures(mols, nmol)
      call center_at_origin(mols(1))
      call get_test_cavity_iswig(mols(1), cavity, err, nleb=nleb_dense_tensor)
      if (allocated(err)) then
         call test_failed(error, "surface construction failed: "//err%message)
         return
      end if
      call cavity%get_gradient(err)
      if (allocated(err)) then
         call test_failed(error, "forward gradient failed: "//err%message)
         return
      end if

      ngrid = cavity%ngrid
      nsph = cavity%nsph
      call make_charges(ngrid, q1, q2)

      allocate (amat(ngrid, ngrid), amat_dense(ngrid, ngrid))
      allocate (amat1_rA(3, nsph, ngrid, ngrid))
      call assemble_pcm_amat(cavity%xi0, cavity%f, cavity%xyz, amat, err)
      if (allocated(err)) then
         call test_failed(error, "assembly failed: "//err%message)
         return
      end if
      call assemble_pcm_amat_with_gradient(cavity%xi0, cavity%f, cavity%xyz, &
                                           cavity%xi1_rA, cavity%f1_rA, cavity%xyz1_rA, &
                                           amat_dense, amat1_rA, err)
      if (allocated(err)) then
         call test_failed(error, "dense gradient assembly failed: "//err%message)
         return
      end if

      do j = 1, ngrid
         do i = 1, ngrid
            call check(error, amat_dense(i, j), amat(i, j), &
                       thr_abs=1.0e-15_wp, thr_rel=1.0e-15_wp, &
                       more="dense-path matrix differs from the plain assembly")
            if (allocated(error)) return
         end do
      end do

      allocate (w_xi(ngrid), w_f(ngrid), w_xyz(3, ngrid))
      call pcm_amat_surface_weights(cavity%xi0, cavity%f, cavity%xyz, q1, q2, &
                                    w_xi, w_f, w_xyz, err)
      if (allocated(err)) then
         call test_failed(error, "surface-weight assembly failed: "//err%message)
         return
      end if

      allocate (grad(3, nsph), grad_ref(3, nsph))
      call pcm_amat_nuclear_gradient(cavity%xi1_rA, cavity%f1_rA, cavity%xyz1_rA, &
                                     w_xi, w_f, w_xyz, grad, err)
      if (allocated(err)) then
         call test_failed(error, "adjoint contraction failed: "//err%message)
         return
      end if

      grad_ref = 0.0_wp
      do j = 1, ngrid
         do i = 1, ngrid
            do iatom = 1, nsph
               do iaxis = 1, 3
                  grad_ref(iaxis, iatom) = grad_ref(iaxis, iatom) &
                                           + q1(i)*q2(j)*amat1_rA(iaxis, iatom, i, j)
               end do
            end do
         end do
      end do

      scale = maxval(abs(grad_ref))
      call check(error, scale > 0.0_wp, &
                 more="reference gradient is identically zero, the test is vacuous")
      if (allocated(error)) return
      dev_grad = maxval(abs(grad - grad_ref))/scale
      call check(error, dev_grad, 0.0_wp, thr=1.0e-13_wp, &
                 more="adjoint gradient differs from the dense-tensor contraction")

   end subroutine test_weight_paths_agree

end module test_model_component_pcm_amat_adjoint
