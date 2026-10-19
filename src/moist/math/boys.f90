!> Fast algorithm for computing the Boys function:
!> Beylkin, G. & Sharma, S., J. Chem. Phys. 155, 174117 (2021)
module moist_math_boys
   use mctc_env, only: wp
   use mctc_io_constants, only: pi

   implicit none(type, external)
   private

   public :: boys0, boys01, boys012, boys_x_far, boys_x_taylor
   public :: dboysfun1, dboysfun12, zboysfun00

   !* -------------------------- Beylkin & Beylkin fast boys -------------------------- *!

   real(wp), parameter :: tol = 1.0E-03_wp
   real(wp), parameter :: sqrtpio2 = 0.886226925452758014_wp
   real(wp), parameter :: t(0:11) = [ &
                          0.20000000000000000E+01_wp, &
                          0.66666666666666663E+00_wp, &
                          0.40000000000000002E+00_wp, &
                          0.28571428571428570E+00_wp, &
                          0.22222222222222221E+00_wp, &
                          0.18181818181818182E+00_wp, &
                          0.15384615384615385E+00_wp, &
                          0.13333333333333333E+00_wp, &
                          0.11764705882352941E+00_wp, &
                          0.10526315789473684E+00_wp, &
                          0.95238095238095233E-01_wp, &
                          0.86956521739130432E-01_wp]
   !> Complex-conjugate pole locations for the rational expansion
   complex(wp), parameter :: zz(1:10) = [ &
                             (0.64304020652330500E+01_wp, 0.18243694739308491E+02_wp), &
                             (0.64304020652330500E+01_wp, -0.18243694739308491E+02_wp), &
                             (-0.12572081889410178E+01_wp, 0.14121366415342502E+02_wp), &
                             (-0.12572081889410178E+01_wp, -0.14121366415342502E+02_wp), &
                             (-0.54103079551670268E+01_wp, 0.10457909575828442E+02_wp), &
                             (-0.54103079551670268E+01_wp, -0.10457909575828442E+02_wp), &
                             (-0.78720025594983341E+01_wp, 0.69309284623985663E+01_wp), &
                             (-0.78720025594983341E+01_wp, -0.69309284623985663E+01_wp), &
                             (-0.92069621609035313E+01_wp, 0.34559308619699376E+01_wp), &
                             (-0.92069621609035313E+01_wp, -0.34559308619699376E+01_wp)]
   !> Residues corresponding to each pole pair
   complex(wp), parameter :: fact(1:10) = [ &
                             (0.13249210991966042E-02_wp, 0.91787356295447745E-03_wp), &
                             (0.13249210991966042E-02_wp, -0.91787356295447745E-03_wp), &
                             (0.55545905103006735E-01_wp, -0.35151540664451613E+01_wp), &
                             (0.55545905103006735E-01_wp, 0.35151540664451613E+01_wp), &
                             (-0.11456407675096416E+03_wp, 0.19213789620924834E+03_wp), &
                             (-0.11456407675096416E+03_wp, -0.19213789620924834E+03_wp), &
                             (0.20915556220686653E+04_wp, -0.15825742912360638E+04_wp), &
                             (0.20915556220686653E+04_wp, 0.15825742912360638E+04_wp), &
                             (-0.94779394228935325E+04_wp, 0.30814443710192086E+04_wp), &
                             (-0.94779394228935325E+04_wp, -0.30814443710192086E+04_wp)]
   !> Weights for the continued-fraction quadrature
   complex(wp), parameter :: ww(1:10) = [ &
                             (-0.83418049867878959E-08_wp, -0.70958810331788253E-08_wp), &
                             (-0.83418050437598581E-08_wp, 0.70958810084577824E-08_wp), &
                             (0.82436739552884774E-07_wp, -0.27704117936134414E-06_wp), &
                             (0.82436739547688584E-07_wp, 0.27704117938414886E-06_wp), &
                             (0.19838416382728666E-05_wp, 0.78321058613942770E-06_wp), &
                             (0.19838416382681279E-05_wp, -0.78321058613180811E-06_wp), &
                             (-0.47372729839268780E-05_wp, 0.58076919074212929E-05_wp), &
                             (-0.47372729839287016E-05_wp, -0.58076919074154416E-05_wp), &
                             (-0.68186014282131608E-05_wp, -0.13515261354290787E-04_wp), &
                             (-0.68186014282138385E-05_wp, 0.13515261354295612E-04_wp)]
   real(wp), parameter :: rzz(1:1) = [ &
                          -0.96321934290343840E+01_wp]
   real(wp), parameter :: rfact(1:1) = [ &
                          0.15247844519077540E+05_wp]
   real(wp), parameter :: rww(1:1) = [ &
                          0.18995875677635889E-04_wp]
   real(wp), parameter :: asymcoef(1:7) = [ &
                          -0.499999999999999799_wp, &
                          0.249999999999993161_wp, &
                          -0.374999999999766599_wp, &
                          0.937499999992027020_wp, &
                          -3.28124999972738868_wp, &
                          14.7656249906697030_wp, &
                          -81.2109371803307752_wp]
   real(wp), parameter :: taylcoef(0:10) = [ &
                          1.0_wp, &
                          -0.333333333333333333_wp, &
                          0.1_wp, &
                          -0.238095238095238095E-01_wp, &
                          0.462962962962962963E-02_wp, &
                          -0.757575757575757576E-03_wp, &
                          0.106837606837606838E-03_wp, &
                          -0.132275132275132275E-04_wp, &
                          1.458916900093370682E-06_wp, &
                          -1.450385222315046877E-07_wp, &
                          1.3122532963802805073E-08_wp]
   real(wp), parameter :: sqpio2 = 0.886226925452758014_wp
   real(wp), parameter :: pp(1:22) = [ &
                          0.001477878263796956477_wp, &
                          0.013317276413725817441_wp, &
                          0.037063591452052541530_wp, &
                          0.072752512422882761543_wp, &
                          0.120236941228785688896_wp, &
                          0.179574293958937717967_wp, &
                          0.253534046984087292596_wp, &
                          0.350388652780721927513_wp, &
                          0.482109575931276669313_wp, &
                          0.663028993158374107103_wp, &
                          0.911814736856590885929_wp, &
                          1.2539502287919293_wp, &
                          1.7244634233573395_wp, &
                          2.3715248262781863_wp, &
                          3.2613796996078355_wp, &
                          4.485130169059591_wp, &
                          6.168062135122484_wp, &
                          8.48247187231787_wp, &
                          11.665305486296793_wp, &
                          16.042417132288328_wp, &
                          22.06192951814709_wp, &
                          30.340112094708307_wp]
   real(wp), parameter :: ff(1:22) = [ &
                          0.0866431027201416556_wp, &
                          0.0857720608434394764_wp, &
                          0.0839350436829178814_wp, &
                          0.0809661970413229146_wp, &
                          0.0769089548492978618_wp, &
                          0.0731552078711821626_wp, &
                          0.0726950035163157228_wp, &
                          0.0752842556089304050_wp, &
                          0.0770943953645196145_wp, &
                          0.0754250625677530441_wp, &
                          0.0689686192650315305_wp, &
                          0.05744480422143023_wp, &
                          0.04208199434694545_wp, &
                          0.025838539448223282_wp, &
                          0.012445024157255563_wp, &
                          0.004292541592599837_wp, &
                          0.0009354342987735969_wp, &
                          0.10840885466502504E-03_wp, &
                          5.271867966761674E-06_wp, &
                          7.765974039750418E-08_wp, &
                          2.2138172422680093E-10_wp, &
                          6.594161760037707E-14_wp]

   !* ------------------------------ Code generated boys ------------------------------ *!

   !> Boys argument below which `boys01` uses its Taylor series
   real(wp), parameter :: boys_x_taylor = 0.5_wp
   !> Boys argument above which erf(sqrt(x)) == 1 and exp(-x) == 0 in double precision
   real(wp), parameter :: boys_x_far = 81.0_wp

contains

   !* -------------------------- Beylkin & Beylkin fast boys -------------------------- *!

   !> Computes real Boys functions F_n(x) for n=0 and n=1
   !>
   !> Uses a short Taylor expansion near x=0 to avoid cancellation in F1, and
   !> an erf-based closed form otherwise
   subroutine dboysfun1(x, vals)
      implicit none(type, external)

      real(wp), intent(in) :: x
      real(wp), intent(out) :: vals(0:1)
      real(wp), parameter :: x_small = 1.0e-6_wp
      real(wp) :: x2, x3, x4, x5, x6, sqrtx, invx, expx

      if (x < 0.0_wp) then
         ! CPCM kernels pass x >= 0; safe fallback for unexpected values
         vals(0) = 1.0_wp
         vals(1) = 1.0_wp/3.0_wp
         return
      end if

      if (x <= x_small) then
         x2 = x*x
         x3 = x2*x
         x4 = x3*x
         x5 = x4*x
         x6 = x5*x

         vals(0) = 1.0_wp &
                   - x/3.0_wp &
                   + x2/10.0_wp &
                   - x3/42.0_wp &
                   + x4/216.0_wp &
                   - x5/1320.0_wp &
                   + x6/9360.0_wp

         vals(1) = 1.0_wp/3.0_wp &
                   - x/5.0_wp &
                   + x2/14.0_wp &
                   - x3/54.0_wp &
                   + x4/264.0_wp &
                   - x5/1560.0_wp &
                   + x6/10920.0_wp
         return
      end if

      sqrtx = sqrt(x)
      vals(0) = sqrtpio2*erf(sqrtx)/sqrtx
      expx = exp(-x)
      invx = 1.0_wp/x
      vals(1) = 0.5_wp*(vals(0) - expx)*invx
   end subroutine dboysfun1

   !> Computes real Boys functions F_n(x) for n=0..12 via quadrature and
   !> upward recursion, switching to asymptotic forms for large |x|
   subroutine dboysfun12(x, vals)
      implicit none(type, external)

      real(wp), intent(in) :: x
      real(wp), intent(out) :: vals(0:12)
      real(wp) :: y, yy, rtmp
      real(wp) ::   p, q, tmp
      integer :: n, k
!
      ! Precompute the exponential factor used throughout the recurrence
      y = exp(-x)
!
      ! Use asymptotic erf form to avoid cancellation for large arguments
      if (abs(x) >= 0.45425955121971775E+01_wp) then
         yy = sqrt(x)
         vals(0) = sqrtpio2*erf(yy)/yy
         yy = y/2.0_wp
         do n = 1, 12
            vals(n) = ((n - 0.5_wp)*vals(n - 1) - yy)/x
         end do
         return
      end if
!
      rtmp = 0
      ! Sum over one member of each conjugate pair and double afterwards
      do k = 1, 10, 2
         rtmp = rtmp + ww(k)*(1.0_wp - fact(k)*y)/(x + zz(k))
      end do
!
      tmp = 0
      do k = 1, 1
         ! Guard the quadrature pole with a Taylor fallback
         if (abs(x + rzz(k)) >= tol) then
            tmp = tmp + rww(k)*(1.0_wp - rfact(k)*y)/(x + rzz(k))
         else
            q = x + rzz(k)
            p = 1.0_wp - q/2.0_wp + q**2/6.0_wp - q**3/24.0_wp + q**4/120.0_wp
            tmp = tmp + rww(k)*p
         end if
      end do
!
      ! Combine conjugate contributions with the real pole correction
      vals(12) = 2*rtmp + tmp
      yy = y/2.0_wp
      do n = 11, 0, -1
         ! Downward recursion for F_n leverages the relation with F_{n+1}
         vals(n) = (x*vals(n + 1) + yy)*t(n)
      end do
!
      return
   end subroutine dboysfun12

   !> Evaluates the complex Boys F_0(z) using asymptotic, Taylor, or
   !> quadrature expansions depending on |z| for numerical stability
   subroutine zboysfun00(z, val)
      implicit none(type, external)

      complex(wp), intent(in) :: z
      complex(wp), intent(out) :: val
      complex(wp) :: z1, ez, y
      integer :: k

      ! z may be complex, so keep the exponential explicit for clarity
      ez = exp(-z)

      ! Large |z| branch relies on the asymptotic form of F_0
      if (abs(z) >= 100.0_wp) then
         z1 = 1.0_wp/sqrt(z)
         y = 1.0_wp/z
         val = asymcoef(7)
         do k = 6, 1, -1
            val = val*y + asymcoef(k)
         end do

         val = ez*val*y + z1*sqpio2

         return
      end if

      ! Series about zero avoids precision loss for small |z|
      if (abs(z) <= 0.35_wp) then
         val = taylcoef(10)
         do k = 9, 0, -1
            val = val*z + taylcoef(k)
         end do
         return
      end if

      ! Intermediate |z| uses the Pade-quadrature rational approximation
      val = sqpio2/sqrt(z) - 0.5_wp*ez*sum(ff(1:22)/(z + pp(1:22)))

      return
   end subroutine zboysfun00

   !* ------------------------------ Code generated boys ------------------------------ *!

   !> Boys function F0(x) = sqrt(pi/(4x))*erf(sqrt(x)) for x >= 0
   !>
   !> Uses the degree-16 Taylor series below boys_x_taylor and the
   !> closed erf form above it.  The saturated branch drops erf entirely
   !>
   !> @param[in]  x   Boys argument (p*r)**2
   !> @param[out] f0  Boys function F0(x)
   pure subroutine boys0(x, f0)
      !> Boys argument
      real(wp), intent(in) :: x
      !> Boys function F0(x)
      real(wp), intent(out) :: f0

      !> Square root of the Boys argument
      real(wp) :: sqrtx

      if (x <= boys_x_taylor) then
         f0 = 1.4483264643598138e-15_wp
         f0 = f0*x - 2.466827010264457e-14_wp
         f0 = f0*x + 3.9554295164585257e-13_wp
         f0 = f0*x - 5.9477940136376354e-12_wp
         f0 = f0*x + 8.35070279514724e-11_wp
         f0 = f0*x - 1.0892221037148573e-09_wp
         f0 = f0*x + 1.3122532963802806e-08_wp
         f0 = f0*x - 1.4503852223150468e-07_wp
         f0 = f0*x + 1.4589169000933706e-06_wp
         f0 = f0*x - 1.3227513227513228e-05_wp
         f0 = f0*x + 0.00010683760683760684_wp
         f0 = f0*x - 0.0007575757575757576_wp
         f0 = f0*x + 0.004629629629629629_wp
         f0 = f0*x - 0.023809523809523808_wp
         f0 = f0*x + 0.1_wp
         f0 = f0*x - 0.3333333333333333_wp
         f0 = f0*x + 1.0_wp
      else
         sqrtx = sqrt(x)
         if (x >= boys_x_far) then
            f0 = 0.88622692545275801_wp/sqrtx
         else
            f0 = 0.88622692545275801_wp*erf(sqrtx)/sqrtx
         end if
      end if

   end subroutine boys0

   !> Boys functions F0(x), F1(x) and the Gaussian factor exp(-x)
   !>
   !> @param[in]  x   Boys argument (p*r)**2
   !> @param[out] f0  Boys function F0(x)
   !> @param[out] f1  Boys function F1(x)
   !> @param[out] ex  Gaussian factor exp(-x)
   pure subroutine boys01(x, f0, f1, ex)
      !> Boys argument
      real(wp), intent(in) :: x
      !> Boys functions F0(x) and F1(x)
      real(wp), intent(out) :: f0, f1
      !> Gaussian factor exp(-x)
      real(wp), intent(out) :: ex

      ! Square root of the Boys argument and the recursion denominator
      real(wp) :: sqrtx, inv2x

      if (x <= boys_x_taylor) then
         ex = exp(-x)
         f0 = 1.4483264643598138e-15_wp
         f0 = f0*x - 2.466827010264457e-14_wp
         f0 = f0*x + 3.9554295164585257e-13_wp
         f0 = f0*x - 5.9477940136376354e-12_wp
         f0 = f0*x + 8.35070279514724e-11_wp
         f0 = f0*x - 1.0892221037148573e-09_wp
         f0 = f0*x + 1.3122532963802806e-08_wp
         f0 = f0*x - 1.4503852223150468e-07_wp
         f0 = f0*x + 1.4589169000933706e-06_wp
         f0 = f0*x - 1.3227513227513228e-05_wp
         f0 = f0*x + 0.00010683760683760684_wp
         f0 = f0*x - 0.0007575757575757576_wp
         f0 = f0*x + 0.004629629629629629_wp
         f0 = f0*x - 0.023809523809523808_wp
         f0 = f0*x + 0.1_wp
         f0 = f0*x - 0.3333333333333333_wp
         f0 = f0*x + 1.0_wp
         f1 = 1.3655649521106815e-15_wp
         f1 = f1*x - 2.317322342975702e-14_wp
         f1 = f1*x + 3.7002405153966856e-13_wp
         f1 = f1*x - 5.537601323041936e-12_wp
         f1 = f1*x + 7.732132217728926e-11_wp
         f1 = f1*x - 1.0020843354176688e-09_wp
         f1 = f1*x + 1.1981443140863431e-08_wp
         f1 = f1*x - 1.3122532963802806e-07_wp
         f1 = f1*x + 1.3053467000835422e-06_wp
         f1 = f1*x - 1.1671335200746965e-05_wp
         f1 = f1*x + 9.259259259259259e-05_wp
         f1 = f1*x - 0.000641025641025641_wp
         f1 = f1*x + 0.003787878787878788_wp
         f1 = f1*x - 0.018518518518518517_wp
         f1 = f1*x + 0.07142857142857142_wp
         f1 = f1*x - 0.2_wp
         f1 = f1*x + 0.3333333333333333_wp
      else if (x >= boys_x_far) then
         ex = 0.0_wp
         sqrtx = sqrt(x)
         f0 = 0.88622692545275801_wp/sqrtx
         f1 = 0.5_wp*f0/x
      else
         sqrtx = sqrt(x)
         ex = exp(-x)
         inv2x = 0.5_wp/x
         f0 = 0.88622692545275801_wp*erf(sqrtx)/sqrtx
         f1 = (f0 - ex)*inv2x
      end if

   end subroutine boys01

   !> Boys functions F0(x), F1(x), F2(x) and the Gaussian factor exp(-x)
   !>
   !> @param[in]  x   Boys argument (p*r)**2
   !> @param[out] f0  Boys function F0(x)
   !> @param[out] f1  Boys function F1(x)
   !> @param[out] f2  Boys function F2(x)
   !> @param[out] ex  Gaussian factor exp(-x)
   pure subroutine boys012(x, f0, f1, f2, ex)
      !> Boys argument
      real(wp), intent(in) :: x
      !> Boys functions F0(x), F1(x) and F2(x)
      real(wp), intent(out) :: f0, f1, f2
      !> Gaussian factor exp(-x)
      real(wp), intent(out) :: ex

      !> Square root of the Boys argument and the recursion denominator
      real(wp) :: sqrtx, inv2x

      if (x <= boys_x_taylor) then
         ex = exp(-x)
         f0 = 1.4483264643598138e-15_wp
         f0 = f0*x - 2.466827010264457e-14_wp
         f0 = f0*x + 3.9554295164585257e-13_wp
         f0 = f0*x - 5.9477940136376354e-12_wp
         f0 = f0*x + 8.35070279514724e-11_wp
         f0 = f0*x - 1.0892221037148573e-09_wp
         f0 = f0*x + 1.3122532963802806e-08_wp
         f0 = f0*x - 1.4503852223150468e-07_wp
         f0 = f0*x + 1.4589169000933706e-06_wp
         f0 = f0*x - 1.3227513227513228e-05_wp
         f0 = f0*x + 0.00010683760683760684_wp
         f0 = f0*x - 0.0007575757575757576_wp
         f0 = f0*x + 0.004629629629629629_wp
         f0 = f0*x - 0.023809523809523808_wp
         f0 = f0*x + 0.1_wp
         f0 = f0*x - 0.3333333333333333_wp
         f0 = f0*x + 1.0_wp
         f1 = 1.3655649521106815e-15_wp
         f1 = f1*x - 2.317322342975702e-14_wp
         f1 = f1*x + 3.7002405153966856e-13_wp
         f1 = f1*x - 5.537601323041936e-12_wp
         f1 = f1*x + 7.732132217728926e-11_wp
         f1 = f1*x - 1.0020843354176688e-09_wp
         f1 = f1*x + 1.1981443140863431e-08_wp
         f1 = f1*x - 1.3122532963802806e-07_wp
         f1 = f1*x + 1.3053467000835422e-06_wp
         f1 = f1*x - 1.1671335200746965e-05_wp
         f1 = f1*x + 9.259259259259259e-05_wp
         f1 = f1*x - 0.000641025641025641_wp
         f1 = f1*x + 0.003787878787878788_wp
         f1 = f1*x - 0.018518518518518517_wp
         f1 = f1*x + 0.07142857142857142_wp
         f1 = f1*x - 0.2_wp
         f1 = f1*x + 0.3333333333333333_wp
         f2 = 1.291750630374969e-15_wp
         f2 = f2*x - 2.1849039233770904e-14_wp
         f2 = f2*x + 3.475983514463553e-13_wp
         f2 = f2*x - 5.18033672155536e-12_wp
         f2 = f2*x + 7.198881719954517e-11_wp
         f2 = f2*x - 9.278558661274711e-10_wp
         f2 = f2*x + 1.1022927689594356e-08_wp
         f2 = f2*x - 1.198144314086343e-07_wp
         f2 = f2*x + 1.1810279667422524e-06_wp
         f2 = f2*x - 1.0442773600668338e-05_wp
         f2 = f2*x + 8.169934640522875e-05_wp
         f2 = f2*x - 0.0005555555555555556_wp
         f2 = f2*x + 0.003205128205128205_wp
         f2 = f2*x - 0.015151515151515152_wp
         f2 = f2*x + 0.05555555555555555_wp
         f2 = f2*x - 0.14285714285714285_wp
         f2 = f2*x + 0.2_wp
      else if (x >= boys_x_far) then
         ex = 0.0_wp
         sqrtx = sqrt(x)
         inv2x = 0.5_wp/x
         f0 = 0.88622692545275801_wp/sqrtx
         f1 = f0*inv2x
         f2 = 3.0_wp*f1*inv2x
      else
         sqrtx = sqrt(x)
         ex = exp(-x)
         inv2x = 0.5_wp/x
         f0 = 0.88622692545275801_wp*erf(sqrtx)/sqrtx
         f1 = (f0 - ex)*inv2x
         f2 = (3.0_wp*f1 - ex)*inv2x
      end if

   end subroutine boys012

end module moist_math_boys
