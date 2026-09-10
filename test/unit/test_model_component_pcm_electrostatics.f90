!> Unit tests for PCM electrostatic nuclear-gradient contractions
module test_model_component_pcm_electrostatics
   use mctc_env, only: wp
   use mctc_env_error, only: moist_error_type => error_type
   use mctc_io, only: structure_type
   use moist_model_component_pcm_electrostatics, only: pcm_electrostatic_nuclear_gradient, &
                                                        pcm_electrostatic_surface_weights, &
                                                        pcm_electrostatic_potential_tangent, &
                                                        pcm_electrostatic_surface_weights_response
   use test_helpers, only: get_test_structures, get_test_points, center_at_origin, &
                           fd4_scalar, fd4_offsets
   use testdrive, only: new_unittest, unittest_type, error_type, check, test_failed
   implicit none (type, external)
   private

   public :: collect_model_component_pcm_electrostatics

   integer, parameter :: nmol = 5
   real(wp), parameter :: fd_atol = 1.0e-10_wp
   real(wp), parameter :: fd_rtol = 1.0e-9_wp

contains

   !> Collect the PCM electrostatics test suite
   subroutine collect_model_component_pcm_electrostatics(testsuite)
      !> Collected unit tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("nuclear_gradient_vs_fd", test_nuclear_gradient_vs_fd), &
                  new_unittest("channels_isolated", test_channels_isolated), &
                  new_unittest("rejects_invalid_shapes", test_rejects_invalid_shapes), &
                  new_unittest("potential_tangent_vs_fd", test_potential_tangent_vs_fd), &
                  new_unittest("surface_weights_response_vs_fd", &
                               test_surface_weights_response_vs_fd) &
                  ]
   end subroutine collect_model_component_pcm_electrostatics

   !> A synthetic surface over a sampled structure, with surface charges, a
   !> charge response, surface tangents and nuclear directions, deterministic
   !>
   !> @param[out] xyz    Surface positions (3, ngrid)
   !> @param[out] sphxyz Nuclear positions (3, nsph)
   !> @param[out] za     Nuclear charges (nsph)
   !> @param[out] q      Surface charges (ngrid)
   !> @param[out] dq     Surface-charge response (ngrid, ndir)
   !> @param[out] d_xyz  Surface-position tangents (3, ngrid, ndir)
   !> @param[out] dirs   Nuclear directions (3, nsph, ndir)
   subroutine second_order_fixture(xyz, sphxyz, za, q, dq, d_xyz, dirs)
      real(wp), allocatable, intent(out) :: xyz(:, :), sphxyz(:, :), za(:), q(:)
      real(wp), allocatable, intent(out) :: dq(:, :), d_xyz(:, :, :), dirs(:, :, :)

      integer, parameter :: n_surface = 12, ndir = 2
      type(structure_type), allocatable :: mols(:)
      integer :: i, iatom, iaxis, v, ngrid, nsph

      call get_test_structures(mols, nmol)
      call center_at_origin(mols(1))
      nsph = mols(1)%nat
      allocate (sphxyz, source=mols(1)%xyz)
      call get_test_points(mols(1), xyz, n_surface)
      ngrid = size(xyz, 2)

      allocate (za(nsph), dirs(3, nsph, ndir))
      do iatom = 1, nsph
         za(iatom) = real(mols(1)%num(mols(1)%id(iatom)), wp)
         do v = 1, ndir
            do iaxis = 1, 3
               dirs(iaxis, iatom, v) = 0.5_wp + 0.4_wp*sin(0.37_wp*iaxis + 0.61_wp*iatom + 1.1_wp*v)
            end do
         end do
      end do
      allocate (q(ngrid), dq(ngrid, ndir), d_xyz(3, ngrid, ndir))
      do i = 1, ngrid
         q(i) = 0.1_wp*sin(0.9_wp*real(i, wp)) - 0.02_wp
         do v = 1, ndir
            dq(i, v) = 0.05_wp*cos(0.4_wp*real(i, wp) + real(v, wp))
            do iaxis = 1, 3
               d_xyz(iaxis, i, v) = 0.3_wp*cos(0.3_wp*real(i, wp) + 0.8_wp*real(iaxis, wp) + real(v, wp))
            end do
         end do
      end do

   end subroutine second_order_fixture

   !> Point-charge potential at the surface, `phi_i = sum_k za_k / |r_i - R_k|`
   pure function point_charge_potential(xyz, sphxyz, za) result(phi)
      real(wp), intent(in) :: xyz(:, :), sphxyz(:, :), za(:)
      real(wp) :: phi(size(xyz, 2))
      integer :: i, k

      phi = 0.0_wp
      do i = 1, size(xyz, 2)
         do k = 1, size(za)
            phi(i) = phi(i) + za(k)/norm2(xyz(:, i) - sphxyz(:, k))
         end do
      end do
   end function point_charge_potential

   !> The moved point-charge potential against differenced potentials
   subroutine test_potential_tangent_vs_fd(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      type(moist_error_type), allocatable :: err
      real(wp), allocatable :: xyz(:, :), sphxyz(:, :), za(:), q(:), dq(:, :)
      real(wp), allocatable :: d_xyz(:, :, :), dirs(:, :, :)
      real(wp), allocatable :: dphi(:, :), vals(:, :), fd(:)
      real(wp) :: s
      integer :: ngrid, ndir, v, k
      real(wp), parameter :: step = 1.0e-3_wp
      character(len=64) :: context

      call second_order_fixture(xyz, sphxyz, za, q, dq, d_xyz, dirs)
      ngrid = size(xyz, 2)
      ndir = size(dirs, 3)

      allocate (dphi(ngrid, ndir), vals(4, ngrid), fd(ngrid))
      call pcm_electrostatic_potential_tangent(xyz, sphxyz, za, d_xyz, dirs, dphi, err)
      if (allocated(err)) then
         call test_failed(error, "potential tangent failed: "//err%message)
         return
      end if

      do v = 1, ndir
         do k = 1, 4
            s = fd4_offsets(k)*step
            vals(k, :) = point_charge_potential(xyz + s*d_xyz(:, :, v), sphxyz + s*dirs(:, :, v), za)
         end do
         fd = fd4_scalar(vals(1, :), vals(2, :), vals(3, :), vals(4, :), step)
         call check(error, maxval(abs(fd)) > 1.0e-4_wp, more="differenced potential is vacuous")
         if (allocated(error)) return
         write (context, "(a,i0)") "potential tangent, direction ", v
         call check(error, maxval(abs(dphi(:, v) - fd)), 0.0_wp, &
                    thr=fd_atol + fd_rtol*maxval(abs(fd)), more=trim(context))
         if (allocated(error)) return
      end do

   end subroutine test_potential_tangent_vs_fd

   !> The adjoint and direct-gradient responses against differenced weights
   subroutine test_surface_weights_response_vs_fd(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      type(moist_error_type), allocatable :: err
      real(wp), allocatable :: xyz(:, :), sphxyz(:, :), za(:), q(:), dq(:, :)
      real(wp), allocatable :: d_xyz(:, :, :), dirs(:, :, :), qefield(:, :)
      real(wp), allocatable :: dw_xyz(:, :, :), dgrad(:, :, :)
      real(wp), allocatable :: w_xyz(:, :, :), grad(:, :, :), fd_w(:, :), fd_g(:, :)
      real(wp) :: s
      integer :: ngrid, nsph, ndir, v, k
      real(wp), parameter :: step = 1.0e-3_wp
      character(len=64) :: context

      call second_order_fixture(xyz, sphxyz, za, q, dq, d_xyz, dirs)
      ngrid = size(xyz, 2)
      nsph = size(za)
      ndir = size(dirs, 3)
      ! The electronic field is host data and is held fixed by both routes
      allocate (qefield(3, ngrid), source=0.0_wp)

      allocate (dw_xyz(3, ngrid, ndir), dgrad(3, nsph, ndir))
      call pcm_electrostatic_surface_weights_response(xyz, sphxyz, q, dq, za, d_xyz, dirs, &
                                                      dw_xyz, dgrad, err)
      if (allocated(err)) then
         call test_failed(error, "adjoint response failed: "//err%message)
         return
      end if

      allocate (w_xyz(4, 3, ngrid), grad(4, 3, nsph), fd_w(3, ngrid), fd_g(3, nsph))
      do v = 1, ndir
         do k = 1, 4
            s = fd4_offsets(k)*step
            call pcm_electrostatic_surface_weights(xyz + s*d_xyz(:, :, v), &
                                                   sphxyz + s*dirs(:, :, v), q + s*dq(:, v), &
                                                   qefield, za, w_xyz(k, :, :), grad(k, :, :), err)
            if (allocated(err)) then
               call test_failed(error, "perturbed adjoint failed: "//err%message)
               return
            end if
         end do
         fd_w = fd4_scalar(w_xyz(1, :, :), w_xyz(2, :, :), w_xyz(3, :, :), w_xyz(4, :, :), step)
         fd_g = fd4_scalar(grad(1, :, :), grad(2, :, :), grad(3, :, :), grad(4, :, :), step)
         call check(error, min(maxval(abs(fd_w)), maxval(abs(fd_g))) > 1.0e-4_wp, &
                    more="differenced adjoint is vacuous")
         if (allocated(error)) return
         write (context, "(a,i0)") "surface-position adjoint response, direction ", v
         call check(error, maxval(abs(dw_xyz(:, :, v) - fd_w)), 0.0_wp, &
                    thr=fd_atol + fd_rtol*maxval(abs(fd_w)), more=trim(context))
         if (allocated(error)) return
         write (context, "(a,i0)") "direct-gradient response, direction ", v
         call check(error, maxval(abs(dgrad(:, :, v) - fd_g)), 0.0_wp, &
                    thr=fd_atol + fd_rtol*maxval(abs(fd_g)), more=trim(context))
         if (allocated(error)) return
      end do

   end subroutine test_surface_weights_response_vs_fd

   !> Check the nuclear/electronic field contraction by finite differences
   subroutine test_nuclear_gradient_vs_fd(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      !> Sampled structures
      type(structure_type), allocatable :: mols(:)
      !> Library error handling
      type(moist_error_type), allocatable :: err
      !> Reference nuclear positions and perturbed copy
      real(wp), allocatable :: sphxyz(:, :), sphxyz_trial(:, :)
      !> Reference surface positions
      real(wp), allocatable :: xyz(:, :)
      !> Surface-position derivatives
      real(wp), allocatable :: xyz1_rA(:, :, :, :)
      !> Surface charges, field weights, and nuclear charges
      real(wp), allocatable :: surface_q(:), qefield(:, :), za(:)
      !> Analytic nuclear gradient
      real(wp), allocatable :: grad(:, :)
      !> Atom, axis, surface, stencil, and extent indices
      integer :: iatom, iaxis, i, k, ngrid, nsph
      !> Stencil values, finite difference, and saved coordinate
      real(wp) :: vals(4), fd, saved
         !> Finite-difference step
      real(wp), parameter :: step = 1.0e-3_wp
      !> Failure context
      character(len=128) :: context
      !> Number of synthetic surface points
      integer, parameter :: n_surface = 12

      call get_test_structures(mols, nmol)
      call center_at_origin(mols(1))
      nsph = mols(1)%nat
      allocate (sphxyz, source=mols(1)%xyz)
      allocate (sphxyz_trial, source=sphxyz)
      call get_test_points(mols(1), xyz, n_surface)
      ngrid = size(xyz, 2)

      allocate (za(nsph))
      do iatom = 1, nsph
         za(iatom) = real(mols(1)%num(mols(1)%id(iatom)), wp)
      end do

      allocate (surface_q(ngrid), qefield(3, ngrid))
      allocate (xyz1_rA(3, 3, nsph, ngrid))
      do i = 1, ngrid
         surface_q(i) = 0.1_wp*sin(0.9_wp*real(i, wp)) - 0.02_wp
         qefield(1, i) = 0.05_wp*cos(0.4_wp*real(i, wp))
         qefield(2, i) = 0.03_wp*sin(1.1_wp*real(i, wp))
         qefield(3, i) = -0.04_wp*cos(0.7_wp*real(i, wp))
         do iatom = 1, nsph
            do iaxis = 1, 3
               xyz1_rA(1, iaxis, iatom, i) = 0.1_wp*sin(0.3_wp*real(i + iatom + iaxis, wp))
               xyz1_rA(2, iaxis, iatom, i) = 0.1_wp*cos(0.5_wp*real(i + 2*iatom - iaxis, wp))
               xyz1_rA(3, iaxis, iatom, i) = 0.1_wp*sin(0.7_wp*real(i - iatom + 3*iaxis, wp))
            end do
         end do
      end do

      allocate (grad(3, nsph))
      call pcm_electrostatic_nuclear_gradient(xyz, sphxyz, xyz1_rA, surface_q, &
                                              qefield, za, grad, err)
      if (allocated(err)) then
         call test_failed(error, "nuclear/electronic contraction failed: "//err%message)
         return
      end if

      do iatom = 1, min(nsph, 3)
         do iaxis = 1, 3
            saved = sphxyz(iaxis, iatom)
            do k = 1, 4
               sphxyz_trial = sphxyz
               sphxyz_trial(iaxis, iatom) = saved + fd4_offsets(k)*step
               vals(k) = nuc_elec_energy(xyz, sphxyz, sphxyz_trial, xyz1_rA, &
                                         surface_q, qefield, za)
            end do
            fd = fd4_scalar(vals(1), vals(2), vals(3), vals(4), step)
            write (context, "(a,i0,a,i0)") "nuc/elec gradient atom ", iatom, &
               " axis ", iaxis
            call check(error, grad(iaxis, iatom), fd, &
                       thr=fd_atol + fd_rtol*abs(fd), more=trim(context))
            if (allocated(error)) return
         end do
      end do
   end subroutine test_nuclear_gradient_vs_fd

   !> Each contraction channel in isolation, against its closed form
   !>
   !> The routine folds two physically distinct terms into one accumulator: the
   !> direct force between the surface charges and the nuclei, and the chain-rule
   !> term carrying the surface's response to nuclear motion. The combined
   !> finite-difference test above would still pass if the two were swapped or
   !> if one absorbed a sign that the other cancelled, so each is switched on
   !> alone here and compared with its closed form to machine precision.
   !> @param[out] error Test failure
   subroutine test_channels_isolated(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      !> Sampled structures
      type(structure_type), allocatable :: mols(:)
      !> Library error handling
      type(moist_error_type), allocatable :: err
      !> Nuclear and surface positions
      real(wp), allocatable :: sphxyz(:, :), xyz(:, :)
      !> Surface-position derivatives
      real(wp), allocatable :: xyz1_rA(:, :, :, :)
      !> Surface charges, field weights, and nuclear charges
      real(wp), allocatable :: surface_q(:), qefield(:, :), za(:)
      !> Analytic and closed-form gradients
      real(wp), allocatable :: grad(:, :), grad_ref(:, :)
      !> Surface, atom, and axis indices; extents
      integer :: i, iatom, iaxis, ngrid, nsph
      !> Displacement and its cubed length
      real(wp) :: rvec(3), r3
      !> Number of synthetic surface points
      integer, parameter :: n_surface = 12
      !> Tolerance for a sum of a few dozen terms in double precision
      real(wp), parameter :: exact_thr = 1.0e-13_wp

      call get_test_structures(mols, nmol)
      call center_at_origin(mols(1))
      nsph = mols(1)%nat
      allocate (sphxyz, source=mols(1)%xyz)
      call get_test_points(mols(1), xyz, n_surface)
      ngrid = size(xyz, 2)

      allocate (za(nsph))
      do iatom = 1, nsph
         za(iatom) = real(mols(1)%num(mols(1)%id(iatom)), wp)
      end do
      allocate (surface_q(ngrid), qefield(3, ngrid))
      allocate (xyz1_rA(3, 3, nsph, ngrid))
      allocate (grad(3, nsph), grad_ref(3, nsph))
      do i = 1, ngrid
         surface_q(i) = 0.1_wp*sin(0.9_wp*real(i, wp)) - 0.02_wp
         qefield(1, i) = 0.05_wp*cos(0.4_wp*real(i, wp))
         qefield(2, i) = 0.03_wp*sin(1.1_wp*real(i, wp))
         qefield(3, i) = -0.04_wp*cos(0.7_wp*real(i, wp))
      end do

      ! Direct channel: a rigid surface (no response) and no electronic field
      ! leave the plain Coulomb force of the surface charges on each nucleus.
      xyz1_rA = 0.0_wp
      call pcm_electrostatic_nuclear_gradient(xyz, sphxyz, xyz1_rA, surface_q, &
                                              0.0_wp*qefield, za, grad, err)
      if (allocated(err)) then
         call test_failed(error, "direct-channel contraction failed: "//err%message)
         return
      end if

      grad_ref = 0.0_wp
      do iatom = 1, nsph
         do i = 1, ngrid
            rvec = xyz(:, i) - sphxyz(:, iatom)
            r3 = norm2(rvec)**3
            grad_ref(:, iatom) = grad_ref(:, iatom) + surface_q(i)*za(iatom)*rvec/r3
         end do
      end do
      call check(error, maxval(abs(grad_ref)) > 0.0_wp, &
                 more="direct channel is identically zero, the test is vacuous")
      if (allocated(error)) return
      call check(error, maxval(abs(grad - grad_ref)), 0.0_wp, &
                 thr=exact_thr*maxval(abs(grad_ref)), &
                 more="direct nuclear force deviates from its closed form")
      if (allocated(error)) return

      ! Chain-rule channel: uncharged nuclei and an uncharged surface leave only
      ! the electronic field contracted with the surface response.
      do i = 1, ngrid
         do iatom = 1, nsph
            do iaxis = 1, 3
               xyz1_rA(1, iaxis, iatom, i) = 0.1_wp*sin(0.3_wp*real(i + iatom + iaxis, wp))
               xyz1_rA(2, iaxis, iatom, i) = 0.1_wp*cos(0.5_wp*real(i + 2*iatom - iaxis, wp))
               xyz1_rA(3, iaxis, iatom, i) = 0.1_wp*sin(0.7_wp*real(i - iatom + 3*iaxis, wp))
            end do
         end do
      end do
      call pcm_electrostatic_nuclear_gradient(xyz, sphxyz, xyz1_rA, &
                                              0.0_wp*surface_q, qefield, &
                                              0.0_wp*za, grad, err)
      if (allocated(err)) then
         call test_failed(error, "chain-rule contraction failed: "//err%message)
         return
      end if

      grad_ref = 0.0_wp
      do i = 1, ngrid
         do iatom = 1, nsph
            do iaxis = 1, 3
               grad_ref(iaxis, iatom) = grad_ref(iaxis, iatom) &
                                        + dot_product(xyz1_rA(:, iaxis, iatom, i), &
                                                      qefield(:, i))
            end do
         end do
      end do
      call check(error, maxval(abs(grad_ref)) > 0.0_wp, &
                 more="chain-rule channel is identically zero, the test is vacuous")
      if (allocated(error)) return
      call check(error, maxval(abs(grad - grad_ref)), 0.0_wp, &
                 thr=exact_thr*maxval(abs(grad_ref)), &
                 more="surface chain rule deviates from its closed form")

   end subroutine test_channels_isolated

   !> Inconsistent array shapes are rejected and the gradient is left defined.
   !> @param[out] error Test failure
   subroutine test_rejects_invalid_shapes(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      !> Library error handling
      type(moist_error_type), allocatable :: err
      !> Synthetic two-nucleus, three-grid point problem
      real(wp) :: xyz(3, 3), sphxyz(3, 2), xyz1_rA(3, 3, 2, 3)
      real(wp) :: surface_q(3), qefield(3, 3), za(2)
      real(wp) :: grad(3, 2), grad_small(3, 1)

      xyz(:, 1) = [1.0_wp, 0.0_wp, 0.0_wp]
      xyz(:, 2) = [0.0_wp, 1.5_wp, 0.3_wp]
      xyz(:, 3) = [-1.2_wp, 0.4_wp, -0.8_wp]
      sphxyz(:, 1) = [0.0_wp, 0.0_wp, 0.0_wp]
      sphxyz(:, 2) = [0.0_wp, 0.0_wp, 2.1_wp]
      xyz1_rA = 0.05_wp
      surface_q = [0.2_wp, -0.1_wp, 0.05_wp]
      qefield = 0.01_wp
      za = [8.0_wp, 1.0_wp]

      ! A gradient sized for the wrong number of nuclei.
      grad_small = 1.0_wp
      call pcm_electrostatic_nuclear_gradient(xyz, sphxyz, xyz1_rA, surface_q, &
                                              qefield, za, grad_small, err)
      call check(error, allocated(err), more="a mis-shaped gradient was accepted")
      if (allocated(error)) return
      call check(error, maxval(abs(grad_small)), 0.0_wp, thr=0.0_wp, &
                 more="rejected contraction left the gradient untouched")
      if (allocated(error)) return
      deallocate (err)

      ! A surface response that does not match the surface it belongs to.
      grad = 1.0_wp
      call pcm_electrostatic_nuclear_gradient(xyz, sphxyz, xyz1_rA(:, :, :, 1:2), &
                                              surface_q, qefield, za, grad, err)
      call check(error, allocated(err), &
                 more="a mis-shaped surface response was accepted")
      if (allocated(error)) return
      call check(error, maxval(abs(grad)), 0.0_wp, thr=0.0_wp, &
                 more="rejected contraction left the gradient untouched")

   end subroutine test_rejects_invalid_shapes

   !> Evaluate the nuclear plus external-field reference energy.
   !> @param[in] xyz Reference grid point positions
   !> @param[in] sphxyz0 Reference nuclear positions
   !> @param[in] sphxyz Displaced nuclear positions
   !> @param[in] xyz1_rA Surface-position response
   !> @param[in] surface_q Surface charges
   !> @param[in] qefield Charge-weighted electronic field
   !> @param[in] za Nuclear charges
   !> @return energy Reference energy
   function nuc_elec_energy(xyz, sphxyz0, sphxyz, xyz1_rA, surface_q, qefield, za) &
      result(energy)
      !> Reference surface positions, reference nuclei, and displaced nuclei
      real(wp), intent(in) :: xyz(:, :), sphxyz0(:, :), sphxyz(:, :)
      !> Surface-position response
      real(wp), intent(in) :: xyz1_rA(:, :, :, :)
      !> Surface charges, electronic field weights, and nuclear charges
      real(wp), intent(in) :: surface_q(:), qefield(:, :), za(:)
      !> Nuclear plus external-field reference energy
      real(wp) :: energy

      !> Surface, atom, axis, and extent indices
      integer :: i, iatom, katom, iaxis, ngrid, nsph
      !> Displaced surface position and response contribution
      real(wp) :: ri(3), shift(3)

      ngrid = size(surface_q)
      nsph = size(za)
      energy = 0.0_wp

      do i = 1, ngrid
         ri = xyz(:, i)
         do iatom = 1, nsph
            do iaxis = 1, 3
               shift = xyz1_rA(:, iaxis, iatom, i)
               ri = ri + shift*(sphxyz(iaxis, iatom) - sphxyz0(iaxis, iatom))
            end do
         end do
         do katom = 1, nsph
            energy = energy + surface_q(i)*za(katom)/norm2(ri - sphxyz(:, katom))
         end do
         energy = energy + dot_product(qefield(:, i), ri)
      end do
   end function nuc_elec_energy

end module test_model_component_pcm_electrostatics
