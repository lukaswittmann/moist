!> Coarse-space deflation (DEF1) for Newton-Krylov linear solves
!> The driver supplies basis profiles and weights; the evaluator fills
!> E = Z^T Omega J Z. With Q = Z E^{-1} Z^T Omega, GMRES solves the
!> projected system (I-JQ)Jx = (I-JQ)b and returns Qb + (I-QJ)x
module moist_model_moz_solver_deflation
   use mctc_env, only: wp
   use moist_context, only: resolve_num_threads
   use moist_math_lapack_getrf, only: wrap_getrf
   use moist_math_lapack_getrs, only: wrap_getrs
   implicit none(type, external)
   private

   public :: moist_model_moz_solver_deflation_type

   !> Fixed coarse space and factorized coarse matrix for DEF1
   type :: moist_model_moz_solver_deflation_type
      !> Whether the solver may use the coarse space
      logical :: enabled = .false.
      !> Residual RMS below which the coarse space is used
      real(wp) :: rmse_activate = 1.0e-3_wp
      !> Basis profiles per site pair (0 = no basis)
      integer :: nfun = 0
      !> Grid points and site pairs the basis was built for (0 = no basis)
      integer :: npts = 0
      integer :: npair = 0
      !> Leading grid points on which the profiles are nonzero
      integer :: ncore = 0
      !> Coarse dimension nfun*npair; coarse index of (profile j, pair p) is
      !> (p - 1)*nfun + j
      integer :: kz = 0
      !> Basis profiles g_j, (npts, nfun); zero beyond ncore
      real(wp), allocatable :: profile(:, :)
      !> Weighted profiles Omega g_j, (npts, nfun); zero beyond ncore
      real(wp), allocatable :: wprofile(:, :)
      !> Coarse matrix E = Z^T Omega J Z (kz, kz), filled by the evaluator;
      !> holds its LU factors once `factorize` succeeds
      real(wp), allocatable :: coarse(:, :)
      !> LU pivots of the coarse matrix
      integer, allocatable :: ipiv(:)
      !> Whether E is factorized at the current Newton point
      logical :: active = .false.
   contains
      procedure :: wanted => deflation_wanted
      procedure :: has_basis => deflation_has_basis
      procedure :: init_basis => deflation_init_basis
      procedure :: factorize => deflation_factorize
      procedure :: coordinates => deflation_coordinates
      procedure :: expand => deflation_expand
      procedure :: deactivate => deflation_deactivate
      procedure :: release => deflation_release
   end type moist_model_moz_solver_deflation_type

contains

   !> Whether this Newton step should be deflated
   !>
   !> @param[in] self Deflation settings
   !> @param[in] rmse Residual RMS at the current Newton point
   !> @return want True when enabled and below the activation threshold
   pure logical function deflation_wanted(self, rmse) result(want)
      !> Deflation settings
      class(moist_model_moz_solver_deflation_type), intent(in) :: self
      !> Residual RMS at the current Newton point
      real(wp), intent(in) :: rmse
      want = self%enabled .and. rmse < self%rmse_activate
   end function deflation_wanted

   !> Whether the basis matches a (npts, npair) iterate
   !>
   !> @param[in] self  Deflation state
   !> @param[in] npts  Grid points of the iterate
   !> @param[in] npair Site pairs of the iterate
   !> @return match True when `init_basis` was called for this shape
   pure logical function deflation_has_basis(self, npts, npair) result(match)
      !> Deflation state
      class(moist_model_moz_solver_deflation_type), intent(in) :: self
      !> Grid points of the iterate
      integer, intent(in) :: npts
      !> Site pairs of the iterate
      integer, intent(in) :: npair
      match = self%npts == npts .and. self%npair == npair .and. allocated(self%coarse)
   end function deflation_has_basis

   !> Take the basis profiles of the coarse space and size the coarse matrix
   !>
   !> - Rows past the last nonzero profile row are skipped by the coarse algebra
   !>
   !> @param[inout] self     Deflation state; receives the basis
   !> @param[in]    profile  basis profiles g_j (npts, nfun)
   !> @param[in]    wprofile weighted profiles Omega g_j (npts, nfun)
   !> @param[in]    npair    site pairs of the iterate
   subroutine deflation_init_basis(self, profile, wprofile, npair)
      !> Deflation state
      class(moist_model_moz_solver_deflation_type), intent(inout) :: self
      !> Basis profiles and weighted profiles
      real(wp), intent(in) :: profile(:, :), wprofile(:, :)
      !> Site pairs of the iterate
      integer, intent(in) :: npair

      !> Grid counter
      integer :: i

      call self%release()
      self%npts = size(profile, 1)
      self%npair = npair
      self%nfun = size(profile, 2)
      self%kz = self%nfun*npair
      self%profile = profile
      self%wprofile = wprofile
      allocate (self%coarse(self%kz, self%kz), self%ipiv(self%kz))
      self%ncore = 0
      do i = self%npts, 1, -1
         if (any(profile(i, :) /= 0.0_wp)) then
            self%ncore = i
            exit
         end if
      end do
      self%active = .false.
   end subroutine deflation_init_basis

   !> Factorize the coarse matrix filled by the evaluator and activate the
   !> deflation for this Newton step; a singular or non-finite E leaves it
   !> inactive; the step then runs undeflated
   !>
   !> @param[inout] self Deflation state (coarse matrix in, LU factors out)
   subroutine deflation_factorize(self)
      !> Deflation state
      class(moist_model_moz_solver_deflation_type), intent(inout) :: self

      !> LAPACK status
      integer :: info

      self%active = .false.
      if (.not. allocated(self%coarse)) return
      if (.not. all(abs(self%coarse) <= huge(1.0_wp))) return
      call wrap_getrf(self%coarse, self%ipiv, info)
      self%active = info == 0 .and. all(abs(self%coarse) <= huge(1.0_wp))
   end subroutine deflation_factorize

   !> Coarse coordinates y = E^{-1} Z^T Omega w of a flat (npts*npair) vector
   !>
   !> @param[in]  self Deflation state with an active factorization
   !> @param[in]  w    Flat vector, pair-major columns of length npts
   !> @param[out] y    Coarse coordinates (kz, 1)
   !> @param[in]  nthreads OpenMP team size; absent or <= 0 takes omp_get_max_threads
   subroutine deflation_coordinates(self, w, y, nthreads)
      !> Deflation state
      class(moist_model_moz_solver_deflation_type), intent(in) :: self
      !> Flat vector
      real(wp), intent(in) :: w(:)
      !> Coarse coordinates
      real(wp), intent(out) :: y(:, :)
      !> OpenMP team size; absent or <= 0 takes omp_get_max_threads
      integer, intent(in), optional :: nthreads
      !> Effective team size
      integer :: nt

      !> Pair, profile and pair-offset indices, LAPACK status
      integer :: p, j, off, info

      nt = merge(resolve_num_threads(nthreads), 1, size(w) >= 65536)
      !$omp parallel do schedule(static) private(j, off) num_threads(nt)
      do p = 1, self%npair
         off = (p - 1)*self%npts
         do j = 1, self%nfun
            y((p - 1)*self%nfun + j, 1) = dot_product(self%wprofile(1:self%ncore, j), &
               & w(off + 1:off + self%ncore))
         end do
      end do
      !$omp end parallel do
      ! Caller supplies active factors and matching shapes; getrs only reports illegal arguments
      call wrap_getrs(self%coarse, y, self%ipiv, info)
   end subroutine deflation_coordinates

   !> Add a coarse-space vector: x := x + alpha Z y
   !>
   !> @param[in]    self  Deflation state
   !> @param[in]    y     Coarse coordinates (kz, 1)
   !> @param[in]    alpha Scale
   !> @param[inout] x     Flat vector, pair-major columns of length npts
   !> @param[in]    nthreads OpenMP team size; absent or <= 0 takes omp_get_max_threads
   subroutine deflation_expand(self, y, alpha, x, nthreads)
      !> Deflation state
      class(moist_model_moz_solver_deflation_type), intent(in) :: self
      !> Coarse coordinates
      real(wp), intent(in) :: y(:, :)
      !> Scale
      real(wp), intent(in) :: alpha
      !> Flat vector
      real(wp), intent(inout) :: x(:)
      !> OpenMP team size; absent or <= 0 takes omp_get_max_threads
      integer, intent(in), optional :: nthreads
      !> Effective team size
      integer :: nt

      !> Pair, profile and pair-offset indices
      integer :: p, j, off

      nt = merge(resolve_num_threads(nthreads), 1, size(x) >= 65536)
      !$omp parallel do schedule(static) private(j, off) num_threads(nt)
      do p = 1, self%npair
         off = (p - 1)*self%npts
         do j = 1, self%nfun
            x(off + 1:off + self%ncore) = x(off + 1:off + self%ncore) + &
               & alpha*y((p - 1)*self%nfun + j, 1)*self%profile(1:self%ncore, j)
         end do
      end do
      !$omp end parallel do
   end subroutine deflation_expand

   !> Switch the deflation off for the current Newton step
   !>
   !> @param[inout] self Deflation state
   subroutine deflation_deactivate(self)
      !> Deflation state
      class(moist_model_moz_solver_deflation_type), intent(inout) :: self
      self%active = .false.
   end subroutine deflation_deactivate

   !> Drop the basis and the coarse matrix; the settings are kept
   !>
   !> @param[inout] self Deflation state
   subroutine deflation_release(self)
      !> Deflation state
      class(moist_model_moz_solver_deflation_type), intent(inout) :: self
      if (allocated(self%profile)) deallocate (self%profile)
      if (allocated(self%wprofile)) deallocate (self%wprofile)
      if (allocated(self%coarse)) deallocate (self%coarse)
      if (allocated(self%ipiv)) deallocate (self%ipiv)
      self%npts = 0
      self%npair = 0
      self%nfun = 0
      self%ncore = 0
      self%kz = 0
      self%active = .false.
   end subroutine deflation_release

end module moist_model_moz_solver_deflation
