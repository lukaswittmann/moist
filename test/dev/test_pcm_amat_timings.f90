!> Timing example isolating the Gaussian PCM A-matrix kernels on polyalanine.
!!
!! The three O(ngrid^2)/O(ngrid^3) kernels of [[moist_model_component_pcm_amat]]
!! are timed separately so an optimization of the pair kernel can be judged
!! against the dense Cholesky solve that follows it:
!!
!!  * `assemble_pcm_amat`         - once per cavity update (once per SCF cycle in
!!                                  the density-dependent isodensity workflow)
!!  * `pcm_amat_surface_weights`  - once per nuclear gradient
!!  * `solve_pcm_cholesky`        - once per charge solve
!!
!! The table also reports the fraction of surface pairs that fall into the
!! saturated regime `u = xi_ij * r_ij >= u_far`, where `erf(u) == 1` to well
!! below machine precision and the kernel collapses to plain `1/r`. That
!! fraction is the head-room available to a near/far split of the pair loop.
module test_pcm_amat_timings
   use, intrinsic :: iso_fortran_env, only: int64
   use mctc_env, only: wp
   use mctc_env_error, only: moist_error_type => error_type
   use mctc_io, only: structure_type
   use mstore, only: get_structure
   use testdrive, only: new_unittest, unittest_type, error_type, check, test_failed
   use moist_cavity, only: cavity_type_iswig, new_cavity_iswig
   use moist_context, only: moist_context_type, new_context
   use moist_radii, only: default_cpcm_radii
   use moist_model_component_pcm_amat, only: assemble_pcm_amat, pcm_amat_surface_weights
   use moist_model_component_pcm_solvers, only: solve_pcm_cholesky
   use moist_math_boys, only: dboysfun1
   use mctc_io_constants, only: pi
!$ use omp_lib, only: omp_get_max_threads
   implicit none
   private

   public :: collect_pcm_amat_timings

   !> Atomic angular grids sampled by the scaling example
   integer, parameter :: nlebs(2) = [110, 194]
   !> Polyalanine chain lengths sampled by the scaling example
   integer, parameter :: nalas(4) = [4, 12, 20, 28]
   !> Water bulk dielectric constant used for the right-hand side scaling
   real(wp), parameter :: eps_bulk = 78.4_wp
   !> Saturation threshold on u = xi_ij * r_ij; erf(9) = 1 - 4e-37
   real(wp), parameter :: u_far = 9.0_wp
   !> Off-diagonal Gaussian Coulomb prefactor of the reference kernel
   real(wp), parameter :: pref = 2.0_wp/sqrt(pi)
   !> Diagonal Gaussian self-interaction prefactor of the reference kernel
   real(wp), parameter :: pref_diag = sqrt(2.0_wp/pi)

contains

   !> Collect the polyalanine A-matrix timing example.
   !> @param[out] testsuite  Collected test-drive timing test
   subroutine collect_pcm_amat_timings(testsuite)
      !> Collected test-drive timing test
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [new_unittest("polyalanine pcm amat kernel timings", &
                                test_polyalanine_amat_timings)]

   end subroutine collect_pcm_amat_timings

   !> Time the three PCM A-matrix kernels over a polyalanine size series.
   !> @param[out] error  Test failure, if a cavity build or a kernel fails
   subroutine test_polyalanine_amat_timings(error)
      !> Test failure, if a cavity build or a kernel fails
      type(error_type), allocatable, intent(out) :: error

      !> Molecular structure of the current chain
      type(structure_type) :: mol
      !> iSwiG cavity carrying the Gaussian surface (xi, f, xyz)
      type(cavity_type_iswig), allocatable :: cav
      !> Shared run context; must outlive the cavity
      type(moist_context_type), target :: ctx
      !> Cavity/kernel error channel
      type(moist_error_type), allocatable :: err
      !> Dense interaction matrix and its reference counterpart
      real(wp), allocatable :: amat(:, :), amat_ref(:, :)
      !> Right-hand side, surface charges, and adjoint weights
      real(wp), allocatable :: rhs(:), q(:), w_xi(:), w_f(:), w_xyz(:, :)
      !> Reference adjoint weights
      real(wp), allocatable :: w_xi_ref(:), w_f_ref(:), w_xyz_ref(:, :)
      !> Wall-clock timings of the production and reference kernels
      real(wp) :: t_cavity, t_amat, t_weights, t_solve, t_amat_ref, t_weights_ref
      !> Relative deviations from the reference kernels
      real(wp) :: dev_amat, dev_w
      !> Fraction of pairs in the saturated (far-field) regime
      real(wp) :: far_frac
      !> Clock samples and tick rate
      integer(int64) :: c0, c1, rate
      !> Grid, structure, and grid-point loop indices
      integer :: ileb, ial, i, ngrid
      !> Structure name buffer
      character(len=12) :: name

      call new_context(ctx, verbosity=0, debug=.false.)
      call system_clock(count_rate=rate)

      write (*, '(a)') ""
      write (*, '(a)') "  Polyalanine Gaussian PCM A-matrix kernel timings"
      write (*, '(a)') "  mstore POLYALANINE; iSwiG cavity; default CPCM radii"
      write (*, '(a)') "  warm wall-clock timings; *_ref is the pre-optimization Boys kernel"
!$    write (*, '(a,i0)') "  OpenMP threads: ", omp_get_max_threads()
      write (*, '(a,f4.1,a)') "  far-field fraction counts pairs with u = xi_ij*r_ij >= ", u_far, &
         & " (erf(u) == 1 to below eps)"
      write (*, '(2x,a10,3a8,5a11,3a9,2a11)') "structure", "nleb", "nat", "ngrid", &
         & "cavity[s]", "amat[s]", "amat_ref", "weight[s]", "wgt_ref", &
         & "amat_x", "wgt_x", "far_frac", "dev_amat", "dev_wgt"
      write (*, '(2x,a10,3a8,5a11,3a9,2a11)') repeat("-", 10), repeat("-", 8), &
         & repeat("-", 8), repeat("-", 8), repeat("-", 11), repeat("-", 11), &
         & repeat("-", 11), repeat("-", 11), repeat("-", 11), repeat("-", 9), &
         & repeat("-", 9), repeat("-", 9), repeat("-", 11), repeat("-", 11)

      do ileb = 1, size(nlebs)
         do ial = 1, size(nalas)
            write (name, '(a,i2.2)') "polyala_", nalas(ial)
            call get_structure(mol, "POLYALANINE", trim(name))

            if (allocated(cav)) deallocate (cav)
            allocate (cav)
            call system_clock(c0)
            call new_cavity_iswig(cav, ctx, nleb=nlebs(ileb), &
                                  radius_model=default_cpcm_radii(), error=err)
            if (allocated(err)) then
               call test_failed(error, trim(name)//": cavity setup failed: "//err%message)
               return
            end if
            call cav%update(mol, error=err)
            if (allocated(err)) then
               call test_failed(error, trim(name)//": cavity update failed: "//err%message)
               return
            end if
            call system_clock(c1)
            t_cavity = elapsed_seconds(c0, c1, rate)

            ngrid = cav%ngrid
            if (allocated(amat)) then
               deallocate (amat, amat_ref, rhs, q, w_xi, w_f, w_xyz)
               deallocate (w_xi_ref, w_f_ref, w_xyz_ref)
            end if
            allocate (amat(ngrid, ngrid), amat_ref(ngrid, ngrid), rhs(ngrid), q(ngrid))
            allocate (w_xi(ngrid), w_f(ngrid), w_xyz(3, ngrid))
            allocate (w_xi_ref(ngrid), w_f_ref(ngrid), w_xyz_ref(3, ngrid))

            ! Warm timings: the first call also pays the first-touch page
            ! faults of the freshly allocated matrix, which the production path
            ! (a cavity matrix reused across SCF cycles) does not.
            call assemble_pcm_amat(cav%xi0, cav%f, cav%xyz, amat, err)
            call system_clock(c0)
            call assemble_pcm_amat(cav%xi0, cav%f, cav%xyz, amat, err)
            call system_clock(c1)
            if (allocated(err)) then
               call test_failed(error, trim(name)//": amat assembly failed: "//err%message)
               return
            end if
            t_amat = elapsed_seconds(c0, c1, rate)

            ! Deterministic, sign-alternating charges: no pair is screened away
            do i = 1, ngrid
               q(i) = (-1.0_wp)**i*real(i, wp)/real(ngrid + 1, wp)
            end do

            call pcm_amat_surface_weights(cav%xi0, cav%f, cav%xyz, q, q, &
                                          w_xi, w_f, w_xyz, err)
            call system_clock(c0)
            call pcm_amat_surface_weights(cav%xi0, cav%f, cav%xyz, q, q, &
                                          w_xi, w_f, w_xyz, err)
            call system_clock(c1)
            if (allocated(err)) then
               call test_failed(error, trim(name)//": surface weights failed: "//err%message)
               return
            end if
            t_weights = elapsed_seconds(c0, c1, rate)

            ! Reference implementation, timed warm the same way, and compared
            ! entry by entry so the speedup is never bought with accuracy. This
            ! runs before the solve, which consumes `q` as its solution vector.
            call reference_assemble(cav%xi0, cav%f, cav%xyz, amat_ref)
            call system_clock(c0)
            call reference_assemble(cav%xi0, cav%f, cav%xyz, amat_ref)
            call system_clock(c1)
            t_amat_ref = elapsed_seconds(c0, c1, rate)

            call reference_weights(cav%xi0, cav%f, cav%xyz, q, w_xi_ref, w_f_ref, w_xyz_ref)
            call system_clock(c0)
            call reference_weights(cav%xi0, cav%f, cav%xyz, q, w_xi_ref, w_f_ref, w_xyz_ref)
            call system_clock(c1)
            t_weights_ref = elapsed_seconds(c0, c1, rate)

            ! Elementwise relative deviation, accumulated column by column so no
            ! ngrid**2 temporary is materialized. A norm-wise measure would be
            ! meaningless here: the diagonal spans ten orders of magnitude
            ! because A_ii = sqrt(2/pi)*xi_i/f_i blows up for a barely-exposed
            ! tessera, which would hide any off-diagonal disagreement.
            dev_amat = 0.0_wp
            do i = 1, ngrid
               dev_amat = max(dev_amat, maxval(abs(amat(:, i) - amat_ref(:, i)) &
                                               /(1.0_wp + abs(amat_ref(:, i)))))
            end do
            dev_w = max(maxval(abs(w_xi - w_xi_ref)), maxval(abs(w_xyz - w_xyz_ref))) &
                    /max(maxval(abs(w_xi_ref)), maxval(abs(w_xyz_ref)))

            rhs = -(eps_bulk - 1.0_wp)/eps_bulk
            call system_clock(c0)
            call solve_pcm_cholesky(amat, rhs, q, err)
            call system_clock(c1)
            if (allocated(err)) then
               call test_failed(error, trim(name)//": Cholesky solve failed: "//err%message)
               return
            end if
            t_solve = elapsed_seconds(c0, c1, rate)

            far_frac = far_field_fraction(cav%xi0, cav%xyz)

            write (*, '(2x,a10,3i8,5f11.4,2f9.2,f9.4,2es11.2)') trim(name), nlebs(ileb), &
               & mol%nat, ngrid, t_cavity, t_amat, t_amat_ref, t_weights, &
               & t_weights_ref, t_amat_ref/max(t_amat, tiny(1.0_wp)), &
               & t_weights_ref/max(t_weights, tiny(1.0_wp)), far_frac, dev_amat, dev_w

            call check(error, ngrid > 0, message=trim(name)//": empty cavity")
            if (allocated(error)) return
            call check(error, all(ieee_is_finite_vec(w_xi)), &
                       message=trim(name)//": non-finite adjoint weights")
            if (allocated(error)) return
            call check(error, dev_amat, 0.0_wp, thr=1.0e-14_wp, &
                       message=trim(name)//": matrix deviates from the reference kernel")
            if (allocated(error)) return
            call check(error, dev_w, 0.0_wp, thr=1.0e-12_wp, &
                       message=trim(name)//": weights deviate from the reference kernel")
            if (allocated(error)) return
         end do
      end do

   end subroutine test_polyalanine_amat_timings

   !> Reference A-matrix assembly: one Boys call per pair, no near/far split.
   !>
   !> This is the straightforward implementation the production kernel replaced.
   !> It is kept here, rather than only in the git history, so the speedup and
   !> the agreement are both measured in the same run on the same cavity.
   !>
   !> @param[in]  xi    Gaussian widths (ngrid)
   !> @param[in]  f     Gaussian switching factors (ngrid)
   !> @param[in]  xyz   Surface positions (3, ngrid)
   !> @param[out] amat  Interaction matrix (ngrid, ngrid)
   subroutine reference_assemble(xi, f, xyz, amat)
      !> Gaussian widths (ngrid)
      real(wp), intent(in) :: xi(:)
      !> Gaussian switching factors (ngrid)
      real(wp), intent(in) :: f(:)
      !> Surface positions (3, ngrid)
      real(wp), intent(in) :: xyz(:, :)
      !> Interaction matrix (ngrid, ngrid)
      real(wp), intent(out) :: amat(:, :)

      !> Surface-point indices and total number of points
      integer :: i, j, ngrid
      !> Pair width, squared-width sum, squared separation, Boys data
      real(wp) :: xi_ij, xi2_sum, r_ij2, t, boys(0:1)

      ngrid = size(xi)
      amat = 0.0_wp
      !$omp parallel do default(none) shared(xi, f, xyz, amat, ngrid) &
      !$omp private(i, j, xi2_sum, xi_ij, r_ij2, t, boys) schedule(static)
      do i = 1, ngrid
         amat(i, i) = pref_diag*xi(i)/f(i)
         do j = 1, i - 1
            xi2_sum = xi(i)*xi(i) + xi(j)*xi(j)
            xi_ij = xi(i)*xi(j)/sqrt(xi2_sum)
            r_ij2 = sum((xyz(:, i) - xyz(:, j))**2)
            t = xi_ij*xi_ij*r_ij2
            call dboysfun1(t, boys)
            amat(i, j) = pref*xi_ij*boys(0)
            amat(j, i) = amat(i, j)
         end do
      end do
      !$omp end parallel do

   end subroutine reference_assemble

   !> Reference adjoint weights: triangular sweep with a scatter to the partner.
   !>
   !> @param[in]  xi     Gaussian widths (ngrid)
   !> @param[in]  f      Gaussian switching factors (ngrid)
   !> @param[in]  xyz    Surface positions (3, ngrid)
   !> @param[in]  q      Contraction vector, used on both sides (ngrid)
   !> @param[out] w_xi   Width weights (ngrid)
   !> @param[out] w_f    Switching-factor weights (ngrid)
   !> @param[out] w_xyz  Position weights (3, ngrid)
   subroutine reference_weights(xi, f, xyz, q, w_xi, w_f, w_xyz)
      !> Gaussian widths (ngrid)
      real(wp), intent(in) :: xi(:)
      !> Gaussian switching factors (ngrid)
      real(wp), intent(in) :: f(:)
      !> Surface positions (3, ngrid)
      real(wp), intent(in) :: xyz(:, :)
      !> Contraction vector, used on both sides (ngrid)
      real(wp), intent(in) :: q(:)
      !> Width and switching-factor weights (ngrid)
      real(wp), intent(out) :: w_xi(:), w_f(:)
      !> Position weights (3, ngrid)
      real(wp), intent(out) :: w_xyz(:, :)

      !> Surface-point indices and total number of points
      integer :: i, j, ngrid
      !> Charge products, pair width, and reusable width denominators
      real(wp) :: qdiag, qsym, xi_ij, xi2_sum, inv_sum, inv_sqrt
      !> Pair-width derivatives, displacement, Boys data, and contracted terms
      real(wp) :: dxi_i, dxi_j, rvec(3), t, boys(0:1), term_xi, term_r

      ngrid = size(xi)
      w_xi = 0.0_wp
      w_f = 0.0_wp
      w_xyz = 0.0_wp
      do i = 1, ngrid
         qdiag = q(i)*q(i)
         w_xi(i) = w_xi(i) + qdiag*pref_diag/f(i)
         w_f(i) = w_f(i) - qdiag*pref_diag*xi(i)/(f(i)*f(i))

         do j = 1, i - 1
            qsym = 2.0_wp*q(i)*q(j)
            xi2_sum = xi(i)*xi(i) + xi(j)*xi(j)
            inv_sum = 1.0_wp/xi2_sum
            inv_sqrt = sqrt(inv_sum)
            xi_ij = xi(i)*xi(j)*inv_sqrt
            dxi_i = xi(j)*inv_sqrt - xi_ij*xi(i)*inv_sum
            dxi_j = xi(i)*inv_sqrt - xi_ij*xi(j)*inv_sum
            rvec = xyz(:, i) - xyz(:, j)
            t = xi_ij*xi_ij*sum(rvec*rvec)
            call dboysfun1(t, boys)
            term_xi = qsym*pref*exp(-t)
            term_r = -2.0_wp*qsym*pref*xi_ij**3*boys(1)
            w_xi(i) = w_xi(i) + term_xi*dxi_i
            w_xi(j) = w_xi(j) + term_xi*dxi_j
            w_xyz(:, i) = w_xyz(:, i) + term_r*rvec
            w_xyz(:, j) = w_xyz(:, j) - term_r*rvec
         end do
      end do

   end subroutine reference_weights

   !> Fraction of surface pairs whose Gaussian argument saturates erf.
   !>
   !> Counts the off-diagonal pairs with u = xi_ij * r_ij >= u_far, i.e. those
   !> for which the kernel is exactly 1/r_ij in double precision.
   !>
   !> @param[in]  xi   Gaussian widths (ngrid)
   !> @param[in]  xyz  Surface positions (3, ngrid)
   !> @return     far  Fraction of off-diagonal pairs in the saturated regime
   function far_field_fraction(xi, xyz) result(far)
      !> Gaussian widths (ngrid)
      real(wp), intent(in) :: xi(:)
      !> Surface positions (3, ngrid)
      real(wp), intent(in) :: xyz(:, :)
      !> Fraction of off-diagonal pairs in the saturated regime
      real(wp) :: far

      !> Surface-point indices and total number of points
      integer :: i, j, ngrid
      !> Saturated-pair counter
      integer(int64) :: nfar
      !> Pair width and squared separation
      real(wp) :: xi_ij, r2

      ngrid = size(xi)
      nfar = 0_int64
      !$omp parallel do default(none) shared(xi, xyz, ngrid) &
      !$omp private(i, j, xi_ij, r2) reduction(+:nfar) schedule(dynamic, 16)
      do i = 1, ngrid
         do j = 1, i - 1
            xi_ij = xi(i)*xi(j)/sqrt(xi(i)*xi(i) + xi(j)*xi(j))
            r2 = sum((xyz(:, i) - xyz(:, j))**2)
            if (xi_ij*xi_ij*r2 >= u_far*u_far) nfar = nfar + 1_int64
         end do
      end do
      !$omp end parallel do

      far = real(nfar, wp)/max(1.0_wp, 0.5_wp*real(ngrid, wp)*real(ngrid - 1, wp))

   end function far_field_fraction

   !> Elementwise finite check usable inside `all(...)`.
   !> @param[in]  x    Values to test
   !> @return     ok   True where the value is finite
   pure function ieee_is_finite_vec(x) result(ok)
      !> Values to test
      real(wp), intent(in) :: x(:)
      !> True where the value is finite
      logical :: ok(size(x))

      ok = (x == x) .and. (abs(x) < huge(1.0_wp))

   end function ieee_is_finite_vec

   !> Convert a system_clock interval to seconds.
   !> @param[in]  c0       Start tick
   !> @param[in]  c1       End tick
   !> @param[in]  rate     Ticks per second
   !> @return     seconds  Elapsed wall-clock time
   pure function elapsed_seconds(c0, c1, rate) result(seconds)
      !> Start tick
      integer(int64), intent(in) :: c0
      !> End tick
      integer(int64), intent(in) :: c1
      !> Ticks per second
      integer(int64), intent(in) :: rate
      !> Elapsed wall-clock time
      real(wp) :: seconds

      seconds = real(c1 - c0, wp)/real(rate, wp)

   end function elapsed_seconds

end module test_pcm_amat_timings
