!> Second-order unit tests for the Gaussian PCM interaction matrix
!>
!> [[moist_model_component_pcm_amat_hessian]] moves the matrix along a batch of
!> surface tangents and differentiates the surface-variable weights of
!> [[pcm_amat_surface_weights]] along the same tangents. Both are compared to
!> 4-point central differences of the first-order routines on a molecular
!> surface where a grid point has partners across both kernel branches; the
!> charge response, the bilinear weight of the two charge vectors, is compared
!> to that weight directly.
module test_model_component_pcm_amat_hessian
   use mctc_env, only: wp
   use mctc_env_error, only: moist_error_type => error_type
   use mctc_io, only: structure_type
   use moist_cavity_iswig, only: cavity_type_iswig
   use moist_model_component_pcm_amat, only: assemble_pcm_amat, pcm_amat_surface_weights, &
                                             pcm_amat_tangent_apply, &
                                             pcm_amat_surface_weights_response
   use test_helpers, only: get_test_structures, center_at_origin, fd4_scalar, fd4_offsets, &
                           get_test_cavity_iswig
   use test_model_component_pcm_amat, only: count_branches, nmol, nleb_survey
   use testdrive, only: new_unittest, unittest_type, error_type, check, test_failed
   implicit none(type, external)
   private

   public :: collect_model_component_pcm_amat_hessian

   !> Number of directions in the tangent batch
   integer, parameter :: ndir = 2

   !> Step along the tangent; the tangents are O(0.1), so the surface moves by
   !> O(1e-4) bohr per unit step
   real(wp), parameter :: step = 1.0e-3_wp

   !> 4-point central FD tolerances
   real(wp), parameter :: fd_atol = 1.0e-9_wp
   real(wp), parameter :: fd_rtol = 1.0e-8_wp

   !> Round-off bound of the charge-response comparison, relative to the
   !> largest weight: the two routines sum the same pair terms in another order
   real(wp), parameter :: exact_rtol = 1.0e-12_wp

   !> Below this a reference carries no information
   real(wp), parameter :: vacuity_thr = 1.0e-4_wp

   !> Smallest switching factor kept, applied as the cavity's own cutoff
   real(wp), parameter :: fd_min_f = 0.1_wp

   !> Smallest exposed surface the tests accept
   integer, parameter :: fd_min_points = 50

contains

   !> Collect the second-order test suite
   subroutine collect_model_component_pcm_amat_hessian(testsuite)
      !> Collected unit tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("tangent_apply_vs_fd", test_tangent_apply_vs_fd), &
                  new_unittest("weights_response_surface_vs_fd", &
                               test_weights_response_surface_vs_fd), &
                  new_unittest("weights_response_charge_exact", &
                               test_weights_response_charge_exact), &
                  new_unittest("rejects_invalid_shapes", test_rejects_invalid_shapes) &
                  ]

   end subroutine collect_model_component_pcm_amat_hessian

   !* --------------------------------- Local helpers --------------------------------- *!

   !> A molecular surface spanning both kernel branches, with charges, a
   !> charge response and a surface tangent batch, all deterministic
   !>
   !> @param[out] xi     Gaussian widths (ngrid)
   !> @param[out] f      Switching factors (ngrid)
   !> @param[out] xyz    Surface positions (3, ngrid)
   !> @param[out] q      Charges (ngrid)
   !> @param[out] dq     Charge response (ngrid, ndir)
   !> @param[out] d_xi   Width tangents (ngrid, ndir)
   !> @param[out] d_f    Switching-factor tangents (ngrid, ndir)
   !> @param[out] d_xyz  Position tangents (3, ngrid, ndir)
   !> @param[out] error  Test failure
   subroutine build_fixture(xi, f, xyz, q, dq, d_xi, d_f, d_xyz, error)
      !> Surface variables
      real(wp), allocatable, intent(out) :: xi(:), f(:), xyz(:, :)
      !> Charges and their response
      real(wp), allocatable, intent(out) :: q(:), dq(:, :)
      !> Surface tangent batch
      real(wp), allocatable, intent(out) :: d_xi(:, :), d_f(:, :), d_xyz(:, :, :)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      type(structure_type), allocatable :: mols(:)
      type(cavity_type_iswig) :: cavity
      type(moist_error_type), allocatable :: err
      integer :: ngrid, nfar, nnear, i, v, iaxis

      call get_test_structures(mols, nmol)
      call center_at_origin(mols(1))
      call get_test_cavity_iswig(mols(1), cavity, err, nleb=nleb_survey, cut_f=fd_min_f)
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

      allocate (q(ngrid), dq(ngrid, ndir), d_xi(ngrid, ndir), d_f(ngrid, ndir), &
                d_xyz(3, ngrid, ndir))
      do i = 1, ngrid
         q(i) = 0.3_wp*sin(0.7_wp*real(i, wp)) + 0.1_wp
         do v = 1, ndir
            dq(i, v) = 0.2_wp*cos(0.31_wp*real(i, wp) + real(v, wp)) - 0.05_wp
            ! Relative width and switching tangents keep both positive under
            ! every stencil point
            d_xi(i, v) = 0.1_wp*xi(i)*cos(0.5_wp*real(i, wp) + real(v, wp))
            d_f(i, v) = 0.1_wp*f(i)*sin(0.9_wp*real(i, wp) - real(v, wp))
            do iaxis = 1, 3
               d_xyz(iaxis, i, v) = 0.3_wp*cos(0.3_wp*real(i, wp) + 0.8_wp*real(iaxis, wp) &
                                               + real(v, wp))
            end do
         end do
      end do

   end subroutine build_fixture

   !> The surface displaced along one direction of the tangent batch
   !>
   !> @param[in]  xi, f, xyz        Base surface
   !> @param[in]  d_xi, d_f, d_xyz  Tangent batch
   !> @param[in]  v                 Direction
   !> @param[in]  s                 Signed step
   !> @param[out] xi_t, f_t, xyz_t  Displaced surface
   subroutine displace(xi, f, xyz, d_xi, d_f, d_xyz, v, s, xi_t, f_t, xyz_t)
      real(wp), intent(in) :: xi(:), f(:), xyz(:, :)
      real(wp), intent(in) :: d_xi(:, :), d_f(:, :), d_xyz(:, :, :)
      integer, intent(in) :: v
      real(wp), intent(in) :: s
      real(wp), allocatable, intent(out) :: xi_t(:), f_t(:), xyz_t(:, :)

      xi_t = xi + s*d_xi(:, v)
      f_t = f + s*d_f(:, v)
      xyz_t = xyz + s*d_xyz(:, :, v)

   end subroutine displace

   !* ------------------------------------- Tests ------------------------------------- *!

   !> The moved matrix applied to the charges against differenced products
   subroutine test_tangent_apply_vs_fd(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      type(moist_error_type), allocatable :: err
      real(wp), allocatable :: xi(:), f(:), xyz(:, :), q(:), dq(:, :)
      real(wp), allocatable :: d_xi(:, :), d_f(:, :), d_xyz(:, :, :)
      real(wp), allocatable :: xi_t(:), f_t(:), xyz_t(:, :), amat(:, :)
      real(wp), allocatable :: daq(:, :), aq(:, :), fd(:)
      integer :: ngrid, v, k
      character(len=64) :: context

      call build_fixture(xi, f, xyz, q, dq, d_xi, d_f, d_xyz, error)
      if (allocated(error)) return
      ngrid = size(xi)

      allocate (daq(ngrid, ndir))
      call pcm_amat_tangent_apply(xi, f, xyz, q, d_xi, d_f, d_xyz, daq, err)
      if (allocated(err)) then
         call test_failed(error, "tangent apply failed: "//err%message)
         return
      end if

      allocate (amat(ngrid, ngrid), aq(4, ngrid), fd(ngrid))
      do v = 1, ndir
         do k = 1, 4
            call displace(xi, f, xyz, d_xi, d_f, d_xyz, v, fd4_offsets(k)*step, xi_t, f_t, xyz_t)
            call assemble_pcm_amat(xi_t, f_t, xyz_t, amat, err)
            if (allocated(err)) then
               call test_failed(error, "perturbed assembly failed: "//err%message)
               return
            end if
            aq(k, :) = matmul(amat, q)
         end do
         fd = fd4_scalar(aq(1, :), aq(2, :), aq(3, :), aq(4, :), step)
         call check(error, maxval(abs(fd)) > vacuity_thr, &
                    more="differenced product is vacuous")
         if (allocated(error)) return
         write (context, "(a,i0)") "moved matrix applied to q, direction ", v
         call check(error, maxval(abs(daq(:, v) - fd)), 0.0_wp, &
                    thr=fd_atol + fd_rtol*maxval(abs(fd)), more=trim(context))
         if (allocated(error)) return
      end do

   end subroutine test_tangent_apply_vs_fd

   !> The surface part of the weight response against differenced weights
   subroutine test_weights_response_surface_vs_fd(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      type(moist_error_type), allocatable :: err
      real(wp), allocatable :: xi(:), f(:), xyz(:, :), q(:), dq(:, :)
      real(wp), allocatable :: d_xi(:, :), d_f(:, :), d_xyz(:, :, :)
      real(wp), allocatable :: xi_t(:), f_t(:), xyz_t(:, :)
      real(wp), allocatable :: dw_xi(:, :), dw_f(:, :), dw_xyz(:, :, :)
      real(wp), allocatable :: w_xi(:, :), w_f(:, :), w_xyz(:, :, :)
      real(wp), allocatable :: fd_xi(:), fd_f(:), fd_xyz(:, :)
      integer :: ngrid, v, k
      character(len=64) :: context

      call build_fixture(xi, f, xyz, q, dq, d_xi, d_f, d_xyz, error)
      if (allocated(error)) return
      ngrid = size(xi)

      ! Charges held fixed: the response is the surface part alone
      dq = 0.0_wp
      allocate (dw_xi(ngrid, ndir), dw_f(ngrid, ndir), dw_xyz(3, ngrid, ndir))
      call pcm_amat_surface_weights_response(xi, f, xyz, q, dq, d_xi, d_f, d_xyz, &
                                             dw_xi, dw_f, dw_xyz, err)
      if (allocated(err)) then
         call test_failed(error, "weight response failed: "//err%message)
         return
      end if

      allocate (w_xi(4, ngrid), w_f(4, ngrid), w_xyz(4, 3, ngrid))
      allocate (fd_xi(ngrid), fd_f(ngrid), fd_xyz(3, ngrid))
      do v = 1, ndir
         do k = 1, 4
            call displace(xi, f, xyz, d_xi, d_f, d_xyz, v, fd4_offsets(k)*step, xi_t, f_t, xyz_t)
            call pcm_amat_surface_weights(xi_t, f_t, xyz_t, q, q, w_xi(k, :), w_f(k, :), &
                                          w_xyz(k, :, :), err)
            if (allocated(err)) then
               call test_failed(error, "perturbed weights failed: "//err%message)
               return
            end if
         end do
         fd_xi = fd4_scalar(w_xi(1, :), w_xi(2, :), w_xi(3, :), w_xi(4, :), step)
         fd_f = fd4_scalar(w_f(1, :), w_f(2, :), w_f(3, :), w_f(4, :), step)
         fd_xyz = fd4_scalar(w_xyz(1, :, :), w_xyz(2, :, :), w_xyz(3, :, :), w_xyz(4, :, :), step)
         call check(error, min(maxval(abs(fd_xi)), maxval(abs(fd_f)), maxval(abs(fd_xyz))) &
                    > vacuity_thr, more="differenced weights are vacuous")
         if (allocated(error)) return

         write (context, "(a,i0)") "width-weight response, direction ", v
         call check(error, maxval(abs(dw_xi(:, v) - fd_xi)), 0.0_wp, &
                    thr=fd_atol + fd_rtol*maxval(abs(fd_xi)), more=trim(context))
         if (allocated(error)) return
         write (context, "(a,i0)") "switching-weight response, direction ", v
         call check(error, maxval(abs(dw_f(:, v) - fd_f)), 0.0_wp, &
                    thr=fd_atol + fd_rtol*maxval(abs(fd_f)), more=trim(context))
         if (allocated(error)) return
         write (context, "(a,i0)") "position-weight response, direction ", v
         call check(error, maxval(abs(dw_xyz(:, :, v) - fd_xyz)), 0.0_wp, &
                    thr=fd_atol + fd_rtol*maxval(abs(fd_xyz)), more=trim(context))
         if (allocated(error)) return
      end do

   end subroutine test_weights_response_surface_vs_fd

   !> The charge part of the weight response is twice the bilinear weight
   subroutine test_weights_response_charge_exact(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      type(moist_error_type), allocatable :: err
      real(wp), allocatable :: xi(:), f(:), xyz(:, :), q(:), dq(:, :)
      real(wp), allocatable :: d_xi(:, :), d_f(:, :), d_xyz(:, :, :)
      real(wp), allocatable :: dw_xi(:, :), dw_f(:, :), dw_xyz(:, :, :)
      real(wp), allocatable :: w_xi(:), w_f(:), w_xyz(:, :)
      real(wp) :: scale
      integer :: ngrid, v
      character(len=64) :: context

      call build_fixture(xi, f, xyz, q, dq, d_xi, d_f, d_xyz, error)
      if (allocated(error)) return
      ngrid = size(xi)

      ! Surface held fixed: the response is the charge part alone
      d_xi = 0.0_wp
      d_f = 0.0_wp
      d_xyz = 0.0_wp
      allocate (dw_xi(ngrid, ndir), dw_f(ngrid, ndir), dw_xyz(3, ngrid, ndir))
      call pcm_amat_surface_weights_response(xi, f, xyz, q, dq, d_xi, d_f, d_xyz, &
                                             dw_xi, dw_f, dw_xyz, err)
      if (allocated(err)) then
         call test_failed(error, "weight response failed: "//err%message)
         return
      end if

      allocate (w_xi(ngrid), w_f(ngrid), w_xyz(3, ngrid))
      do v = 1, ndir
         call pcm_amat_surface_weights(xi, f, xyz, dq(:, v), q, w_xi, w_f, w_xyz, err)
         if (allocated(err)) then
            call test_failed(error, "bilinear weights failed: "//err%message)
            return
         end if
         scale = max(maxval(abs(w_xi)), maxval(abs(w_f)), maxval(abs(w_xyz)))
         call check(error, scale > vacuity_thr, more="bilinear weights are vacuous")
         if (allocated(error)) return
         write (context, "(a,i0)") "charge response, direction ", v
         call check(error, max(maxval(abs(dw_xi(:, v) - 2.0_wp*w_xi)), &
                               maxval(abs(dw_f(:, v) - 2.0_wp*w_f)), &
                               maxval(abs(dw_xyz(:, :, v) - 2.0_wp*w_xyz))), 0.0_wp, &
                    thr=exact_rtol*scale, more=trim(context))
         if (allocated(error)) return
      end do

   end subroutine test_weights_response_charge_exact

   !> Malformed batches are refused and leave a defined output behind
   subroutine test_rejects_invalid_shapes(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      type(moist_error_type), allocatable :: err
      real(wp) :: xi(3), f(3), xyz(3, 3), q(3)
      real(wp) :: d_xi(3, 2), d_f(3, 2), d_xyz(3, 3, 2), dq(3, 2), dq_bad(2, 2)
      real(wp) :: daq(3, 2), daq_bad(3, 1)
      real(wp) :: dw_xi(3, 2), dw_f(3, 2), dw_xyz(3, 3, 2), dw_xyz_bad(3, 3, 1)

      xi = [0.5_wp, 0.6_wp, 0.7_wp]
      f = [0.9_wp, 0.8_wp, 0.7_wp]
      xyz = reshape([0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 1.5_wp, 0.0_wp], [3, 3])
      q = [0.1_wp, -0.2_wp, 0.3_wp]
      d_xi = 0.1_wp
      d_f = 0.1_wp
      d_xyz = 0.1_wp
      dq = 0.1_wp
      dq_bad = 0.1_wp
      daq = 1.0_wp
      daq_bad = 1.0_wp
      dw_xyz_bad = 1.0_wp

      call pcm_amat_tangent_apply(xi, f, xyz, q, d_xi, d_f, d_xyz, daq_bad, err)
      call check(error, allocated(err), more="output shape mismatch was accepted")
      if (allocated(error)) return
      call check(error, maxval(abs(daq_bad)), 0.0_wp, thr=0.0_wp, &
                 more="rejected tangent apply left its output undefined")
      if (allocated(error)) return
      deallocate (err)

      call pcm_amat_surface_weights_response(xi, f, xyz, q, dq_bad, d_xi, d_f, d_xyz, &
                                             dw_xi, dw_f, dw_xyz, err)
      call check(error, allocated(err), more="charge-response shape mismatch was accepted")
      if (allocated(error)) return
      deallocate (err)

      call pcm_amat_surface_weights_response(xi, f, xyz, q, dq, d_xi, d_f, d_xyz, &
                                             dw_xi, dw_f, dw_xyz_bad, err)
      call check(error, allocated(err), more="output shape mismatch was accepted")
      if (allocated(error)) return
      call check(error, maxval(abs(dw_xyz_bad)), 0.0_wp, thr=0.0_wp, &
                 more="rejected weight response left its output undefined")
      if (allocated(error)) return
      deallocate (err)

      call pcm_amat_tangent_apply(xi, f, xyz, q, d_xi, d_f, d_xyz, daq, err)
      call check(error, .not. allocated(err), more="a valid batch was refused")

   end subroutine test_rejects_invalid_shapes

end module test_model_component_pcm_amat_hessian
