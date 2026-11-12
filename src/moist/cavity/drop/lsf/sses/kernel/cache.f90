!> Per-thread geometry and converged moment cache for SSES
module moist_cavity_drop_lsf_sses_cache
   use mctc_env_accuracy, only: wp
   use moist_cavity_drop_lsf_sses_algebra, only: sses_jet_type, jet_constant
   use moist_cavity_drop_lsf_sses_native, only: moment_state_type, moment_init, &
                                                moment_eval_native, moment_third_native
   use moist_cavity_drop_lsf_sses_kernel, only: sses_jets
   implicit none(type, external)
   private
   public :: sses_cache_type, sses_cache_build, sses_cached_evaluate, sses_cache_matches, sses_cache_moment
   !> Value-owned arrays and native state, safe under intrinsic thread-clone assignment
   type :: sses_cache_type
      !> Successful finite evaluation
      logical :: valid = .false.
      !> Probe radius in Bohr
      real(wp) :: probe = -1.0_wp
      !> Smoothing length in Bohr
      real(wp) :: epsilon = -1.0_wp
      !> Boltzmann temperature in Bohr^2
      real(wp) :: tau = -1.0_wp
      !> Atom centers in Bohr
      real(wp), allocatable :: centers(:, :)
      !> Bare radii in Bohr
      real(wp), allocatable :: radii(:)
      !> Native moment state
      type(moment_state_type) :: state
   end type sses_cache_type
contains
   !> Compare geometry and physical parameters with the immutable snapshot
   !>
   !> @param[in] cache Per-thread moment cache
   !> @param[in] centers Atom centers in Bohr
   !> @param[in] radii Bare radii in Bohr
   !> @param[in] probe Probe radius in Bohr
   !> @param[in] epsilon Smoothing length in Bohr
   !> @param[in] tau Boltzmann temperature in Bohr^2
   function sses_cache_matches(cache, centers, radii, probe, epsilon, tau) result(matches)
      !> Per-thread moment cache
      type(sses_cache_type), intent(in) :: cache
      !> Atom centers in Bohr
      real(wp), intent(in) :: centers(:, :)
      !> Bare radii in Bohr
      real(wp), intent(in) :: radii(:)
      !> Probe radius in Bohr
      real(wp), intent(in) :: probe
      !> Smoothing length in Bohr
      real(wp), intent(in) :: epsilon
      !> Boltzmann temperature in Bohr^2
      real(wp), intent(in) :: tau
      !> Exact geometry snapshot match
      logical :: matches
      matches = .false.
      if (.not. cache%valid) return
      if (size(radii) /= size(cache%radii)) return
      if (any(shape(centers) /= shape(cache%centers))) return
      matches = all(centers == cache%centers) .and. all(radii == cache%radii) &
                .and. probe == cache%probe .and. epsilon == cache%epsilon .and. tau == cache%tau
   end function sses_cache_matches
   !> Initialize a geometry snapshot and value-owned moment state
   !>
   !> @param[out] cache Per-thread moment cache
   !> @param[in] centers Atom centers in Bohr
   !> @param[in] radii Bare radii in Bohr
   !> @param[in] probe Probe radius in Bohr
   !> @param[in] epsilon Smoothing length in Bohr
   !> @param[in] tau Boltzmann temperature in Bohr^2
   !> @param[out] valid Successful finite evaluation
   subroutine sses_cache_build(cache, centers, radii, probe, epsilon, tau, valid)
      !> Per-thread moment cache
      type(sses_cache_type), intent(out) :: cache
      !> Atom centers in Bohr
      real(wp), intent(in) :: centers(:, :)
      !> Bare radii in Bohr
      real(wp), intent(in) :: radii(:)
      !> Probe radius in Bohr
      real(wp), intent(in) :: probe
      !> Smoothing length in Bohr
      real(wp), intent(in) :: epsilon
      !> Boltzmann temperature in Bohr^2
      real(wp), intent(in) :: tau
      !> Successful finite evaluation
      logical, intent(out) :: valid
      !> Native status code
      integer :: status
      valid = .false.
      cache%centers = centers
      cache%radii = radii
      cache%probe = probe
      cache%epsilon = epsilon
      cache%tau = tau
      status = moment_init(cache%state, size(radii), centers, radii, probe, epsilon, tau)
      valid = status == 0
      cache%valid = valid
   end subroutine sses_cache_build
   !> Evaluate only the requested spatial order, reusing converged point state
   !>
   !> @param[in,out] cache Per-thread moment cache
   !> @param[in] point Evaluation point in Bohr
   !> @param[in] order Highest total derivative order
   !> @param[out] result Normalized Taylor coefficients
   !> @param[out] valid Successful finite evaluation
   subroutine sses_cached_evaluate(cache, point, order, result, valid)
      !> Per-thread moment cache
      type(sses_cache_type), intent(inout) :: cache
      !> Evaluation point in Bohr
      real(wp), intent(in) :: point(3)
      !> Highest total derivative order
      integer, intent(in) :: order
      !> Normalized Taylor coefficients
      type(sses_jet_type), intent(out) :: result
      !> Successful finite evaluation
      logical, intent(out) :: valid
      real(wp) :: f, g(3), h(3, 3), third(3, 3, 3), stats(3)
      real(wp) :: ps(3, 5), cs(3, size(cache%radii), 5), rs(size(cache%radii), 5)
      integer :: status, i
      valid = .false.
      result = jet_constant(0.0_wp, order)
      if (.not. cache%valid .or. order < 0 .or. order > 4) return
      if (size(cache%radii) > 1) then
         status = moment_eval_native(cache%state, point, min(order, 2), f, g, h, stats)
         if (status /= 0) return
      end if
      if (order == 4 .or. size(cache%radii) == 1) then
         ps = 0.0_wp
         cs = 0.0_wp
         rs = 0.0_wp
         do i = 1, 3
            ps(i, i) = 1.0_wp
         end do
         call sses_jets(point, cache%centers, cache%radii, cache%probe, cache%epsilon, cache%tau, &
                        ps, cs, rs, order, cache%state%y, result, valid)
         return
      end if
      result%c(0) = f
      if (order >= 1) result%c(1:3) = g
      if (order >= 2) then
         result%c(6) = h(1, 1)/2.0_wp
         result%c(7) = h(1, 2)
         result%c(8) = h(1, 3)
         result%c(11) = h(2, 2)/2.0_wp
         result%c(12) = h(2, 3)
         result%c(15) = h(3, 3)/2.0_wp
      end if
      if (order >= 3) then
         status = moment_third_native(cache%state, point, third, stats)
         if (status /= 0) return
         result%c(21) = third(1, 1, 1)/6.0_wp
         result%c(22) = third(1, 1, 2)/2.0_wp
         result%c(23) = third(1, 1, 3)/2.0_wp
         result%c(26) = third(1, 2, 2)/2.0_wp
         result%c(27) = third(1, 2, 3)
         result%c(30) = third(1, 3, 3)/2.0_wp
         result%c(36) = third(2, 2, 2)/6.0_wp
         result%c(37) = third(2, 2, 3)/2.0_wp
         result%c(40) = third(2, 3, 3)/2.0_wp
         result%c(46) = third(3, 3, 3)/6.0_wp
      end if
      result%constant = order == 0
      valid = .true.
   end subroutine sses_cached_evaluate
   !> Converged moment variable of the cached point, for warm Taylor solves
   !>
   !> @param[in] cache Per-thread moment cache
   !> @param[in] point Evaluation point in Bohr
   !> @param[out] y Moment optimizer variable
   function sses_cache_moment(cache, point, y) result(available)
      !> Per-thread moment cache
      type(sses_cache_type), intent(in) :: cache
      !> Evaluation point in Bohr
      real(wp), intent(in) :: point(3)
      !> Moment optimizer variable
      real(wp), intent(out) :: y(3)
      !> Cached moment belongs to point
      logical :: available
      y = cache%state%y
      available = cache%valid .and. cache%state%valid .and. all(point == cache%state%x)
   end function sses_cache_moment
end module moist_cavity_drop_lsf_sses_cache
