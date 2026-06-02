!> Provides an inverface to symmetric divide-and-conquer solver

!> Wrapper to symmetric divide-and-conquer solver for general eigenvalue problems
module moist_math_lapack_sygvd
   use mctc_env, only: sp, dp, error_type, fatal_error
   use moist_output_format, only: format_string
   use moist_math_lapack_type, only: solver_type, context_solver
   implicit none
   private

   public :: new_sygvd

   interface lapack_sygvd
      pure subroutine ssygvd(itype, jobz, uplo, n, a, lda, b, ldb, w, work, lwork, &
            & iwork, liwork, info)
         import :: sp
         real(sp), intent(inout) :: a(lda, *)
         real(sp), intent(inout) :: b(ldb, *)
         real(sp), intent(out) :: w(*)
         integer, intent(in) :: itype
         character(len=1), intent(in) :: jobz
         character(len=1), intent(in) :: uplo
         integer, intent(out) :: info
         integer, intent(in) :: n
         integer, intent(in) :: lda
         integer, intent(in) :: ldb
         real(sp), intent(inout) :: work(*)
         integer, intent(in) :: lwork
         integer, intent(inout) :: iwork(*)
         integer, intent(in) :: liwork
      end subroutine ssygvd
      pure subroutine dsygvd(itype, jobz, uplo, n, a, lda, b, ldb, w, work, lwork, &
            & iwork, liwork, info)
         import :: dp
         real(dp), intent(inout) :: a(lda, *)
         real(dp), intent(inout) :: b(ldb, *)
         real(dp), intent(out) :: w(*)
         integer, intent(in) :: itype
         character(len=1), intent(in) :: jobz
         character(len=1), intent(in) :: uplo
         integer, intent(out) :: info
         integer, intent(in) :: n
         integer, intent(in) :: lda
         integer, intent(in) :: ldb
         real(dp), intent(inout) :: work(*)
         integer, intent(in) :: lwork
         integer, intent(inout) :: iwork(*)
         integer, intent(in) :: liwork
      end subroutine dsygvd
   end interface lapack_sygvd

   !> Wrapper class for solving symmetric general eigenvalue problems
   type, public, extends(solver_type) :: sygvd_solver
      private
      integer :: n = 0
      integer, allocatable :: iwork(:)
      real(sp), allocatable :: swork(:)
      real(sp), allocatable :: sbmat(:, :)
      real(dp), allocatable :: dwork(:)
      real(dp), allocatable :: dbmat(:, :)
   contains
      procedure :: solve_sp
      procedure :: solve_dp
   end type sygvd_solver

contains

   subroutine new_sygvd(self, ndim)
      type(sygvd_solver), intent(out) :: self
      integer, intent(in) :: ndim
      self%n = ndim
   end subroutine new_sygvd

   subroutine solve_sp(self, hmat, smat, eval, error)
      class(sygvd_solver), intent(inout) :: self
      real(sp), contiguous, intent(inout) :: hmat(:, :)
      real(sp), contiguous, intent(in) :: smat(:, :)
      real(sp), contiguous, intent(inout) :: eval(:)
      type(error_type), allocatable, intent(out) :: error
      integer :: info, lswork, liwork

      if (self%n == 0) then
         self%n = size(hmat, 1)
      end if
      if (.not. allocated(self%swork)) then
         allocate (self%swork(1 + 6*self%n + 2*self%n**2))
      end if
      if (.not. allocated(self%iwork)) then
         allocate (self%iwork(3 + 5*self%n))
      end if
      self%sbmat = smat
      lswork = size(self%swork)
      liwork = size(self%iwork)

      call lapack_sygvd(1, 'v', 'u', self%n, hmat, self%n, self%sbmat, self%n, eval, &
         & self%swork, lswork, self%iwork, liwork, info)

      call handle_info(error, info)

   end subroutine solve_sp

   subroutine solve_dp(self, hmat, smat, eval, error)
      class(sygvd_solver), intent(inout) :: self
      real(dp), contiguous, intent(inout) :: hmat(:, :)
      real(dp), contiguous, intent(in) :: smat(:, :)
      real(dp), contiguous, intent(inout) :: eval(:)
      type(error_type), allocatable, intent(out) :: error
      integer :: info, ldwork, liwork

      if (self%n == 0) then
         self%n = size(hmat, 1)
      end if
      if (.not. allocated(self%dwork)) then
         allocate (self%dwork(1 + 6*self%n + 2*self%n**2))
      end if
      if (.not. allocated(self%iwork)) then
         allocate (self%iwork(3 + 5*self%n))
      end if
      self%dbmat = smat
      ldwork = size(self%dwork)
      liwork = size(self%iwork)

      call lapack_sygvd(1, 'v', 'u', self%n, hmat, self%n, self%dbmat, self%n, eval, &
         & self%dwork, ldwork, self%iwork, liwork, info)

      call handle_info(error, info)

   end subroutine solve_dp

   subroutine handle_info(error, info)
      type(error_type), allocatable, intent(out) :: error
      integer, intent(in) :: info

      if (info /= 0) then
         call fatal_error(error, "(sygvd) failed to solve eigenvalue problem.&
            & info="//format_string(info, '(i0)'))
      end if
   end subroutine handle_info

end module moist_math_lapack_sygvd
