!*****************************************************************************************
!>
!  Wrapper for LUSOL. Will eventually be moved into a separate repo.

module moist_math_linalg_lusol_ez

   use moist_math_linalg_lusol, only: lu1fac, lu6sol
   use mctc_env_accuracy, only: wp, ip => i8

   implicit none

   private

   type, public :: lusol_settings
      integer(ip) :: nout = 6
      integer(ip) :: lprint = 0
      integer(ip) :: maxcol = 5
      integer(ip) :: method = 0 ! TPP
      integer(ip) :: keepLU = 1

      real(wp) :: Ltol1 = 100.0_wp
      real(wp) :: Ltol2 = 10.0_wp
      real(wp) :: small = epsilon(1.0_wp)**0.8_wp
      real(wp) :: Utol1 = epsilon(1.0_wp)**0.67_wp
      real(wp) :: Utol2 = epsilon(1.0_wp)**0.67_wp
      real(wp) :: Uspace = 3.0_wp
      real(wp) :: dens1 = 0.3_wp
      real(wp) :: dens2 = 0.5_wp

      integer(ip) :: mode = 5 ! for [[lu6sol]] : `w` solves `A w = v`.
   end type lusol_settings

   public :: solve

contains
!*****************************************************************************************

!*****************************************************************************************
!>
!  Wrapper for [[lu1fac]] + [[lu6sol]] to solve a linear system `A*x = b`.

   subroutine solve(n_cols, n_rows, n_nonzero, irow, icol, mat, b, x, istat, settings)

      integer, intent(in) :: n_cols !! `n`: number of columns in A.
      integer, intent(in) :: n_rows !! `m`: number of rows in A.
      integer, intent(in) :: n_nonzero !! number of nonzero elements of A.
      integer, dimension(n_nonzero), intent(in) :: irow, icol !! sparsity pattern (size is `n_nonzero`)
      real(wp), dimension(n_nonzero), intent(in) :: mat !! matrix elements (size is `n_nonzero`)
      real(wp), dimension(n_rows), intent(in) :: b !! right hand side (size is `m`)
      real(wp), dimension(n_cols), intent(out) :: x !! solution !size is `n`
      integer, intent(out) :: istat !! status code
      type(lusol_settings), intent(in), optional :: settings !! settings (if not present, defaults are used)

      integer(ip) :: nelem, n, m
      integer(ip) :: lena
      real(wp), dimension(:), allocatable :: a
      integer(ip), dimension(:), allocatable :: indc
      integer(ip), dimension(:), allocatable :: indr
      real(wp), dimension(:), allocatable :: ww
      real(wp), dimension(:), allocatable :: w
      real(wp), dimension(:), allocatable :: v
      integer(ip) :: inform
      integer(ip) :: luparm(30)
      real(wp) :: parmlu(30)
      integer(ip), dimension(:), allocatable :: p, q, &
                                                lenc, lenr, &
                                                iploc, iqloc, &
                                                ipinv, iqinv, &
                                                locc, locr

      type(lusol_settings) :: options

      if (present(settings)) options = settings ! use user-supplied settings

      n = n_cols
      m = n_rows
      nelem = n_nonzero
      lena = 1 + max(2*nelem, 10*m, 10*n, 10000_ip)

      allocate (a(lena))
      allocate (indc(lena))
      allocate (indr(lena))
      associate (n => n_cols, m => n_rows)
         allocate (p(m), q(n), lenc(n), lenr(m), &
                   iploc(n), iqloc(m), ipinv(m), &
                   iqinv(n), locc(n), locr(m))
         allocate (ww(n))
         allocate (w(n)) ! x
         allocate (v(m)) ! b
      end associate

      a = 0; indc = 0; indr = 0
      a(1:nelem) = mat
      indc(1:nelem) = irow
      indr(1:nelem) = icol

      ! settings
      luparm = 0
      luparm(1) = options%nout
      luparm(2) = options%lprint
      luparm(3) = options%maxcol
      luparm(6) = options%method
      luparm(8) = options%keepLU
      parmlu = 0
      parmlu(1) = options%Ltol1
      parmlu(2) = options%Ltol2
      parmlu(3) = options%small
      parmlu(4) = options%Utol1
      parmlu(5) = options%Utol2
      parmlu(6) = options%Uspace
      parmlu(7) = options%dens1
      parmlu(8) = options%dens2

      call lu1fac(m, n, nelem, lena, luparm, parmlu, &
                  a, indc, indr, p, q, &
                  lenc, lenr, locc, locr, &
                  iploc, iqloc, ipinv, iqinv, ww, inform)

      !write(*,*) 'lu1fac inform = ', inform

      !  inform = 0 if the LU factors were obtained successfully.
      !         = 1 if U appears to be singular, as judged by lu6chk.
      !         = 3 if some index pair indc(l), indr(l) lies outside
      !             the matrix dimensions 1:m , 1:n.
      !         = 4 if some index pair indc(l), indr(l) duplicates
      !             another such pair.
      !         = 7 if the arrays a, indc, indr were not large enough.
      !             Their length "lena" should be increase to at least
      !             the value "minlen" given in luparm(13).
      !         = 8 if there was some other fatal error.  (Shouldn't happen!)
      !         = 9 if no diagonal pivot could be found with TSP or TDP.
      !             The matrix must not be sufficiently definite
      !             or quasi-definite.
      !         =10 if there was some other fatal error.

      v = b ! right hand side

      ! solve `A w = v`.
      call lu6sol(options%mode, m, n, v, w, &
                  lena, luparm, parmlu, &
                  a, indc, indr, p, q, &
                  lenc, lenr, locc, locr, &
                  inform)

      !write(*,*) 'lu6sol inform = ', inform

      ! On exit, inform = 0 except as follows.
      ! If mode = 3,4,5,6 and if U (and hence A) is singular, then
      ! inform = 1 if there is a nonzero residual in solving the system
      ! involving U.  parmlu(20) returns the norm of the residual.

      x = w ! solution
      istat = int(inform)

   end subroutine solve
!*****************************************************************************************

!*****************************************************************************************
end module moist_math_linalg_lusol_ez
!*****************************************************************************************
