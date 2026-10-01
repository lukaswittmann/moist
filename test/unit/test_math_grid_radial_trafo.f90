!> Test suite for the radial Fourier-Bessel transforms and their factory
!>
!>   - Factory: tag dispatch and every mismatch error
!>   - DST-IV: Gaussian against its closed-form transform, round trips
!>   - General quadrature: convergence in the node count for r, k <= 10,
!>     sinc(0) = 1, agreement with DST-IV on a uniform pair
!>   - Both: adjoint dot products, batched against scalar, empty batches,
!>     size errors, independent clones; the base default batched loops
!>     through a minimal test trafo; one trafo per thread is tested in
!>     `math_grid_3d_threaded`, outside test-drive's team
module test_math_grid_radial_trafo
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use mctc_env, only: wp
   use mctc_env_error, only: mctc_error => error_type
   use mctc_io_constants, only: pi
   use testdrive, only: new_unittest, unittest_type, error_type, check, test_failed
   use moist_math_grid_radial_rule, only: moist_math_grid_radial_rule_chebyshev2_type, &
      & new_chebyshev2_rule
   use moist_math_grid_radial_mapping, only: moist_math_grid_radial_mapping_becke_type, &
      & new_becke_mapping
   use moist_math_grid_radial_grid, only: moist_math_grid_radial_type, &
      & moist_math_grid_radial_recipe_type, new_radial_grid, new_uniform_radial_pair, &
      & transform_quadrature, transform_dst4
   use moist_math_grid_radial_trafo, only: moist_math_grid_radial_trafo_type, &
      & moist_math_grid_radial_trafo_dst4_type, moist_math_grid_radial_trafo_quadrature_type, &
      & new_quadrature_trafo, new_radial_trafo
   implicit none(type, external)
   private

   public :: collect_math_grid_radial_trafo

   !> Transform directions, selectors for `apply_one` and `apply_all`
   integer, parameter :: op_r2k = 1, op_k2r = 2, op_r2k_adj = 3, op_k2r_adj = 4

   !> Minimal trafo that only scales, so the base default batched loops are exercised
   type, extends(moist_math_grid_radial_trafo_type) :: scale_trafo_type
   contains
      !> Transform tag (none)
      procedure :: implementation => scale_implementation
      !> 2*f
      procedure :: fbt_r2k => scale_r2k
      !> 3*f
      procedure :: fbt_k2r => scale_k2r
      !> 5*f
      procedure :: fbt_r2k_adj => scale_r2k_adj
      !> 7*f
      procedure :: fbt_k2r_adj => scale_k2r_adj
   end type scale_trafo_type

contains

   !> Collect all math_grid_radial_trafo tests
   !>
   !> @param[out] testsuite  Collected unit tests
   subroutine collect_math_grid_radial_trafo(testsuite)
      !> Collected unit tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
         new_unittest("factory_selects_dst4", test_factory_dst4), &
         new_unittest("factory_selects_quadrature", test_factory_quadrature), &
         new_unittest("factory_unset_tag", test_factory_unset_tag), &
         new_unittest("factory_unknown_tag", test_factory_unknown_tag), &
         new_unittest("factory_mixed_tags", test_factory_mixed_tags), &
         new_unittest("factory_dst4_npts_mismatch", test_factory_dst4_npts), &
         new_unittest("factory_dst4_spacing_mismatch", test_factory_dst4_spacing), &
         new_unittest("dst4_gaussian_analytic", test_dst4_gaussian), &
         new_unittest("dst4_round_trip", test_dst4_round_trip), &
         new_unittest("quadrature_gaussian_convergence", test_quadrature_convergence), &
         new_unittest("quadrature_sinc_at_zero", test_quadrature_sinc_zero), &
         new_unittest("quadrature_matches_dst4_on_uniform_pair", test_quadrature_vs_dst4), &
         new_unittest("adjoint_dot_product", test_adjoint), &
         new_unittest("batched_matches_scalar", test_batched), &
         new_unittest("empty_batch", test_empty_batch), &
         new_unittest("size_errors", test_size_errors), &
         new_unittest("base_default_batched_loops", test_base_default), &
         new_unittest("clone_independent", test_clone) &
         ]
   end subroutine collect_math_grid_radial_trafo

   ! --------------------------------------------------------------------------
   ! Helpers
   ! --------------------------------------------------------------------------

   !> Deterministic, well-conditioned test signal
   !>
   !> @param[in] i  1-based sample index
   pure function sample(i) result(x)
      !> Sample index
      integer, intent(in) :: i
      !> Sample value
      real(wp) :: x

      x = sin(1.7_wp*real(i, wp)) + 0.3_wp*cos(0.31_wp*real(i*i, wp)) &
         & + 0.1_wp*real(modulo(i, 7) - 3, wp)
   end function sample

   !> Largest deviation max|got - ref| relative to max|ref|; huge on a non-finite entry
   !>
   !> @param[in] got  Values under test
   !> @param[in] ref  Reference values
   pure function rel_dev(got, ref) result(err)
      !> Values under test
      real(wp), intent(in) :: got(:)
      !> Reference values
      real(wp), intent(in) :: ref(:)
      !> Relative max-norm deviation
      real(wp) :: err

      real(wp) :: scal

      if (.not. all(ieee_is_finite(got)) .or. .not. all(ieee_is_finite(ref))) then
         err = huge(1.0_wp)
         return
      end if
      if (size(ref) == 0) then
         err = 0.0_wp
         return
      end if
      scal = maxval(abs(ref))
      if (scal <= 0.0_wp) scal = 1.0_wp
      err = maxval(abs(got - ref))/scal
   end function rel_dev

   !> `rel_dev` over all columns of a 2D array
   !>
   !> @param[in] got  Values under test
   !> @param[in] ref  Reference values
   pure function rel_dev2(got, ref) result(err)
      !> Values under test
      real(wp), intent(in) :: got(:, :)
      !> Reference values
      real(wp), intent(in) :: ref(:, :)
      !> Relative max-norm deviation
      real(wp) :: err

      err = rel_dev(reshape(got, [size(got)]), reshape(ref, [size(ref)]))
   end function rel_dev2

   !> Chebyshev-II + Becke pair with fixed scales, the dense-trafo test pair
   !>
   !> @param[out] rgrid  r-space grid
   !> @param[out] kgrid  k-space grid
   !> @param[in]  nr     Number of r nodes
   !> @param[in]  p_r    r-space scale (bohr)
   !> @param[in]  nk     Number of k nodes
   !> @param[in]  p_k    k-space scale (1/bohr)
   !> @param[out] merr   Construction error
   subroutine make_cheb_pair(rgrid, kgrid, nr, p_r, nk, p_k, merr)
      !> r-space grid
      type(moist_math_grid_radial_type), intent(out) :: rgrid
      !> k-space grid
      type(moist_math_grid_radial_type), intent(out) :: kgrid
      !> Number of r nodes
      integer, intent(in) :: nr
      !> r-space scale
      real(wp), intent(in) :: p_r
      !> Number of k nodes
      integer, intent(in) :: nk
      !> k-space scale
      real(wp), intent(in) :: p_k
      !> Construction error
      type(mctc_error), allocatable, intent(out) :: merr

      call make_cheb_grid(rgrid, nr, p_r, merr)
      if (allocated(merr)) return
      call make_cheb_grid(kgrid, nk, p_k, merr)
   end subroutine make_cheb_pair

   !> One Chebyshev-II + Becke grid with a fixed scale
   !>
   !> @param[out] grid  Radial grid
   !> @param[in]  n     Number of nodes
   !> @param[in]  p     Becke scale
   !> @param[out] merr  Construction error
   subroutine make_cheb_grid(grid, n, p, merr)
      !> Radial grid
      type(moist_math_grid_radial_type), intent(out) :: grid
      !> Number of nodes
      integer, intent(in) :: n
      !> Becke scale
      real(wp), intent(in) :: p
      !> Construction error
      type(mctc_error), allocatable, intent(out) :: merr

      type(moist_math_grid_radial_rule_chebyshev2_type) :: rule
      type(moist_math_grid_radial_mapping_becke_type) :: becke
      type(moist_math_grid_radial_recipe_type) :: recipe

      call new_chebyshev2_rule(rule)
      call new_becke_mapping(becke, merr, scale=p)
      if (allocated(merr)) return
      allocate (recipe%rule, source=rule)
      allocate (recipe%mapping, source=becke)
      recipe%npts = n
      call new_radial_grid(grid, recipe, 1, merr)
   end subroutine make_cheb_grid

   !> Build a trafo with the factory, turning a library error into a test failure
   !>
   !> @param[out] error  Test failure
   !> @param[in]  rgrid  r-space grid
   !> @param[in]  kgrid  k-space grid
   !> @param[out] trafo  New trafo
   subroutine make_trafo(error, rgrid, kgrid, trafo)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> r-space grid
      type(moist_math_grid_radial_type), intent(in) :: rgrid
      !> k-space grid
      type(moist_math_grid_radial_type), intent(in) :: kgrid
      !> New trafo
      class(moist_math_grid_radial_trafo_type), allocatable, intent(out) :: trafo

      type(mctc_error), allocatable :: merr

      call new_radial_trafo(trafo, rgrid, kgrid, merr)
      if (allocated(merr)) call test_failed(error, merr%message)
   end subroutine make_trafo

   !> Uniform pair plus its trafo, or a Chebyshev-II + Becke pair plus its trafo
   !>
   !> @param[out] error  Test failure
   !> @param[in]  kind   transform_dst4 or transform_quadrature
   !> @param[out] rgrid  r-space grid
   !> @param[out] kgrid  k-space grid
   !> @param[out] trafo  New trafo
   subroutine make_case(error, kind, rgrid, kgrid, trafo)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> transform_dst4 or transform_quadrature
      integer, intent(in) :: kind
      !> r-space grid
      type(moist_math_grid_radial_type), intent(out) :: rgrid
      !> k-space grid
      type(moist_math_grid_radial_type), intent(out) :: kgrid
      !> New trafo
      class(moist_math_grid_radial_trafo_type), allocatable, intent(out) :: trafo

      type(mctc_error), allocatable :: merr

      if (kind == transform_dst4) then
         call new_uniform_radial_pair(rgrid, kgrid, 96, 0.1_wp, merr)
      else
         call make_cheb_pair(rgrid, kgrid, 64, 1.2_wp, 80, 1.5_wp, merr)
      end if
      if (allocated(merr)) then
         call test_failed(error, merr%message)
         return
      end if
      call make_trafo(error, rgrid, kgrid, trafo)
   end subroutine make_case

   !> Apply one scalar transform direction
   !>
   !> @param[in,out] trafo  Trafo instance
   !> @param[in]     op     Direction selector
   !> @param[in]     f_in   Input field
   !> @param[out]    f_out  Output field
   !> @param[out]    merr   Library error
   subroutine apply_one(trafo, op, f_in, f_out, merr)
      !> Trafo instance
      class(moist_math_grid_radial_trafo_type), intent(inout) :: trafo
      !> Direction selector
      integer, intent(in) :: op
      !> Input field
      real(wp), intent(in) :: f_in(:)
      !> Output field
      real(wp), intent(out) :: f_out(:)
      !> Library error
      type(mctc_error), allocatable, intent(out) :: merr

      select case (op)
      case (op_r2k)
         call trafo%fbt_r2k(f_in, f_out, merr)
      case (op_k2r)
         call trafo%fbt_k2r(f_in, f_out, merr)
      case (op_r2k_adj)
         call trafo%fbt_r2k_adj(f_in, f_out, merr)
      case default
         call trafo%fbt_k2r_adj(f_in, f_out, merr)
      end select
   end subroutine apply_one

   !> Apply one batched transform direction
   !>
   !> @param[in,out] trafo  Trafo instance
   !> @param[in]     op     Direction selector
   !> @param[in]     f_in   Input fields, one per column
   !> @param[out]    f_out  Output fields, one per column
   !> @param[out]    merr   Library error
   subroutine apply_all(trafo, op, f_in, f_out, merr)
      !> Trafo instance
      class(moist_math_grid_radial_trafo_type), intent(inout) :: trafo
      !> Direction selector
      integer, intent(in) :: op
      !> Input fields
      real(wp), intent(in) :: f_in(:, :)
      !> Output fields
      real(wp), intent(out) :: f_out(:, :)
      !> Library error
      type(mctc_error), allocatable, intent(out) :: merr

      select case (op)
      case (op_r2k)
         call trafo%fbt_r2k_all(f_in, f_out, merr)
      case (op_k2r)
         call trafo%fbt_k2r_all(f_in, f_out, merr)
      case (op_r2k_adj)
         call trafo%fbt_r2k_adj_all(f_in, f_out, merr)
      case default
         call trafo%fbt_k2r_adj_all(f_in, f_out, merr)
      end select
   end subroutine apply_all

   !> Input and output lengths of a direction: r-space nr, k-space nk
   !>
   !> @param[in]  trafo  Trafo instance
   !> @param[in]  op     Direction selector
   !> @param[out] n_in   Input length
   !> @param[out] n_out  Output length
   pure subroutine op_sizes(trafo, op, n_in, n_out)
      !> Trafo instance
      class(moist_math_grid_radial_trafo_type), intent(in) :: trafo
      !> Direction selector
      integer, intent(in) :: op
      !> Input length
      integer, intent(out) :: n_in
      !> Output length
      integer, intent(out) :: n_out

      if (op == op_r2k .or. op == op_k2r_adj) then
         n_in = trafo%nr
         n_out = trafo%nk
      else
         n_in = trafo%nk
         n_out = trafo%nr
      end if
   end subroutine op_sizes

   !> Deterministic input fields of shape (n, nb)
   !>
   !> @param[in]  n       Column length
   !> @param[in]  nb      Number of columns
   !> @param[in]  offset  Signal offset
   !> @param[out] f       Fields
   pure subroutine fill_fields(n, nb, offset, f)
      !> Column length
      integer, intent(in) :: n
      !> Number of columns
      integer, intent(in) :: nb
      !> Signal offset
      integer, intent(in) :: offset
      !> Fields
      real(wp), allocatable, intent(out) :: f(:, :)

      integer :: i, j

      allocate (f(n, nb))
      do j = 1, nb
         do i = 1, n
            f(i, j) = sample(i + offset + 17*j)
         end do
      end do
   end subroutine fill_fields

   !> Run a factory call that must fail with a message containing `expected`
   !>
   !> @param[out] error     Test failure
   !> @param[in]  rgrid     r-space grid
   !> @param[in]  kgrid     k-space grid
   !> @param[in]  expected  Substring of the expected message
   !> @param[in]  what      Case label for the failure message
   subroutine expect_factory_error(error, rgrid, kgrid, expected, what)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> r-space grid
      type(moist_math_grid_radial_type), intent(in) :: rgrid
      !> k-space grid
      type(moist_math_grid_radial_type), intent(in) :: kgrid
      !> Substring of the expected message
      character(len=*), intent(in) :: expected
      !> Case label
      character(len=*), intent(in) :: what

      class(moist_math_grid_radial_trafo_type), allocatable :: trafo
      type(mctc_error), allocatable :: merr

      call new_radial_trafo(trafo, rgrid, kgrid, merr)
      if (.not. allocated(merr)) then
         call test_failed(error, what//": factory accepted the pair")
         return
      end if
      call check(error, index(merr%message, expected) > 0, &
         & what//": unexpected message '"//merr%message//"'")
      if (allocated(error)) return
      call check(error, .not. allocated(trafo), what//": trafo allocated despite the error")
   end subroutine expect_factory_error

   ! --------------------------------------------------------------------------
   ! Factory
   ! --------------------------------------------------------------------------

   !> A uniform pair selects DST-IV and exposes the recorded spacings
   !>
   !> @param[out] error  Test failure
   subroutine test_factory_dst4(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      type(moist_math_grid_radial_type) :: rgrid, kgrid
      class(moist_math_grid_radial_trafo_type), allocatable :: trafo
      type(mctc_error), allocatable :: merr

      call new_uniform_radial_pair(rgrid, kgrid, 64, 0.05_wp, merr)
      if (allocated(merr)) then
         call test_failed(error, merr%message)
         return
      end if
      call make_trafo(error, rgrid, kgrid, trafo)
      if (allocated(error)) return
      call check(error, trafo%implementation() == transform_dst4, &
         & "Uniform pair must select the DST-IV trafo")
      if (allocated(error)) return
      call check(error, trafo%nr == 64 .and. trafo%nk == 64, "DST-IV trafo sizes")
      if (allocated(error)) return
      select type (trafo)
      type is (moist_math_grid_radial_trafo_dst4_type)
         call check(error, trafo%dr == rgrid%spacing .and. trafo%dk == kgrid%spacing, &
            & "DST-IV trafo must expose the recorded dr and dk")
      class default
         call test_failed(error, "Uniform pair did not produce a DST-IV trafo type")
      end select
   end subroutine test_factory_dst4

   !> A Chebyshev-II + Becke pair selects the general quadrature; nr and nk may differ
   !>
   !> @param[out] error  Test failure
   subroutine test_factory_quadrature(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      type(moist_math_grid_radial_type) :: rgrid, kgrid
      class(moist_math_grid_radial_trafo_type), allocatable :: trafo
      type(mctc_error), allocatable :: merr

      call make_cheb_pair(rgrid, kgrid, 40, 1.0_wp, 50, 1.5_wp, merr)
      if (allocated(merr)) then
         call test_failed(error, merr%message)
         return
      end if
      call make_trafo(error, rgrid, kgrid, trafo)
      if (allocated(error)) return
      call check(error, trafo%implementation() == transform_quadrature, &
         & "Chebyshev-II + Becke pair must select the quadrature trafo")
      if (allocated(error)) return
      call check(error, trafo%nr == 40 .and. trafo%nk == 50, "Quadrature trafo sizes")
      if (allocated(error)) return
      select type (trafo)
      type is (moist_math_grid_radial_trafo_quadrature_type)
      class default
         call test_failed(error, "Quadrature pair did not produce a quadrature trafo type")
      end select
   end subroutine test_factory_quadrature

   !> A default-initialized grid (tag 0) is refused on either side
   !>
   !> @param[out] error  Test failure
   subroutine test_factory_unset_tag(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      type(moist_math_grid_radial_type) :: unset, rgrid, kgrid
      type(mctc_error), allocatable :: merr

      call expect_factory_error(error, unset, unset, "r grid has no transform tag", "both unset")
      if (allocated(error)) return
      call new_uniform_radial_pair(rgrid, kgrid, 16, 0.1_wp, merr)
      if (allocated(merr)) then
         call test_failed(error, merr%message)
         return
      end if
      call expect_factory_error(error, rgrid, unset, "k grid has no transform tag", "k unset")
      if (allocated(error)) return
      call expect_factory_error(error, unset, kgrid, "r grid has no transform tag", "r unset")
   end subroutine test_factory_unset_tag

   !> A tag that is neither transform_quadrature nor transform_dst4 is refused
   !>
   !> @param[out] error  Test failure
   subroutine test_factory_unknown_tag(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      type(moist_math_grid_radial_type) :: rgrid, kgrid
      type(mctc_error), allocatable :: merr

      call make_cheb_pair(rgrid, kgrid, 16, 1.0_wp, 16, 1.0_wp, merr)
      if (allocated(merr)) then
         call test_failed(error, merr%message)
         return
      end if
      kgrid%transform = 99
      call expect_factory_error(error, rgrid, kgrid, "k grid has an unknown transform tag", &
         & "unknown k tag")
   end subroutine test_factory_unknown_tag

   !> Exactly one DST-IV grid is refused, in either order
   !>
   !> @param[out] error  Test failure
   subroutine test_factory_mixed_tags(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      type(moist_math_grid_radial_type) :: ur, uk, cr, ck
      type(mctc_error), allocatable :: merr

      call new_uniform_radial_pair(ur, uk, 32, 0.1_wp, merr)
      if (.not. allocated(merr)) call make_cheb_pair(cr, ck, 32, 1.0_wp, 32, 1.5_wp, merr)
      if (allocated(merr)) then
         call test_failed(error, merr%message)
         return
      end if
      call expect_factory_error(error, ur, ck, "exactly one grid is tagged for DST-IV", &
         & "DST-IV r, quadrature k")
      if (allocated(error)) return
      call expect_factory_error(error, cr, uk, "exactly one grid is tagged for DST-IV", &
         & "quadrature r, DST-IV k")
   end subroutine test_factory_mixed_tags

   !> Two DST-IV grids from pairs of different size are refused
   !>
   !> @param[out] error  Test failure
   subroutine test_factory_dst4_npts(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      type(moist_math_grid_radial_type) :: r1, k1, r2, k2
      type(mctc_error), allocatable :: merr

      call new_uniform_radial_pair(r1, k1, 32, 0.1_wp, merr)
      if (.not. allocated(merr)) call new_uniform_radial_pair(r2, k2, 48, 0.1_wp, merr)
      if (allocated(merr)) then
         call test_failed(error, merr%message)
         return
      end if
      call expect_factory_error(error, r1, k2, "differ in npts", "npts 32 vs 48")
   end subroutine test_factory_dst4_npts

   !> Two DST-IV grids whose recorded spacings do not satisfy dk = pi/(npts*dr) are refused
   !>
   !> Covers a k grid from a pair of another dr and a recorded dk off by one ULP
   !>
   !> @param[out] error  Test failure
   subroutine test_factory_dst4_spacing(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      type(moist_math_grid_radial_type) :: r1, k1, r2, k2
      type(mctc_error), allocatable :: merr

      call new_uniform_radial_pair(r1, k1, 32, 0.1_wp, merr)
      if (.not. allocated(merr)) call new_uniform_radial_pair(r2, k2, 32, 0.2_wp, merr)
      if (allocated(merr)) then
         call test_failed(error, merr%message)
         return
      end if
      call expect_factory_error(error, r1, k2, "k spacing does not match", "dr 0.1 vs 0.2")
      if (allocated(error)) return
      k1%spacing = nearest(k1%spacing, 1.0_wp)
      call expect_factory_error(error, r1, k1, "k spacing does not match", "dk one ULP off")
   end subroutine test_factory_dst4_spacing

   ! --------------------------------------------------------------------------
   ! DST-IV
   ! --------------------------------------------------------------------------

   !> DST-IV transforms of Gaussians match the closed-form 3D Fourier transform
   !>
   !> f(r) = exp(-a r^2), F(k) = (pi/a)^(3/2) exp(-k^2/(4a)); scalar and
   !> batched, two grids that resolve the Gaussian on both sides; deviation
   !> relative to the reference maximum
   !>
   !> @param[out] error  Test failure
   subroutine test_dst4_gaussian(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: sizes(2) = [256, 512]
      real(wp), parameter :: spacings(2) = [0.05_wp, 0.03_wp]
      real(wp), parameter :: expo(2) = [0.7_wp, 0.91_wp]
      real(wp), parameter :: thr = 3.0e-15_wp
      real(wp), parameter :: pi_dst = acos(-1.0_wp)
      type(moist_math_grid_radial_type) :: rgrid, kgrid
      class(moist_math_grid_radial_trafo_type), allocatable :: trafo
      type(mctc_error), allocatable :: merr
      real(wp), allocatable :: f(:, :), fk(:, :), got(:, :)
      integer :: ic, is, n

      do is = 1, size(sizes)
         n = sizes(is)
         call new_uniform_radial_pair(rgrid, kgrid, n, spacings(is), merr)
         if (allocated(merr)) then
            call test_failed(error, merr%message)
            return
         end if
         call make_trafo(error, rgrid, kgrid, trafo)
         if (allocated(error)) return
         allocate (f(n, 2), fk(n, 2), got(n, 2))
         do ic = 1, 2
            f(:, ic) = exp(-expo(ic)*rgrid%r**2)
            fk(:, ic) = (pi_dst/expo(ic))**1.5_wp*exp(-kgrid%r**2/(4.0_wp*expo(ic)))
         end do

         call trafo%fbt_r2k(f(:, 1), got(:, 1), merr)
         if (.not. allocated(merr)) call trafo%fbt_k2r(fk(:, 1), got(:, 2), merr)
         if (allocated(merr)) then
            call test_failed(error, merr%message)
            return
         end if
         call check(error, rel_dev(got(:, 1), fk(:, 1)) <= thr, &
            & "DST-IV fbt_r2k of a Gaussian deviates from its analytic transform")
         if (allocated(error)) return
         call check(error, rel_dev(got(:, 2), f(:, 1)) <= thr, &
            & "DST-IV fbt_k2r of the analytic transform deviates from the Gaussian")
         if (allocated(error)) return

         call trafo%fbt_r2k_all(f, got, merr)
         if (allocated(merr)) then
            call test_failed(error, merr%message)
            return
         end if
         call check(error, rel_dev(got(:, 1), fk(:, 1)) <= thr .and. &
            & rel_dev(got(:, 2), fk(:, 2)) <= thr, &
            & "DST-IV fbt_r2k_all of Gaussians deviates from the analytic transforms")
         if (allocated(error)) return
         call trafo%fbt_k2r_all(fk, got, merr)
         if (allocated(merr)) then
            call test_failed(error, merr%message)
            return
         end if
         call check(error, rel_dev(got(:, 1), f(:, 1)) <= thr .and. &
            & rel_dev(got(:, 2), f(:, 2)) <= thr, &
            & "DST-IV fbt_k2r_all of analytic transforms deviates from the Gaussians")
         if (allocated(error)) return
         deallocate (f, fk, got)
      end do
   end subroutine test_dst4_gaussian

   !> DST-IV round trips are the identity: k2r(r2k(f)) = f and r2k_adj(k2r_adj(g)) = g
   !>
   !> The weight diagonals telescope through the involution to exactly 1
   !> for dk = pi/(npts*dr); sizes include one node and odd lengths
   !>
   !> @param[out] error  Test failure
   subroutine test_dst4_round_trip(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: sizes(4) = [1, 7, 64, 255]
      type(moist_math_grid_radial_type) :: rgrid, kgrid
      class(moist_math_grid_radial_trafo_type), allocatable :: trafo
      type(mctc_error), allocatable :: merr
      real(wp), allocatable :: f(:, :), g(:, :), back(:, :)
      integer :: is, n

      do is = 1, size(sizes)
         n = sizes(is)
         call new_uniform_radial_pair(rgrid, kgrid, n, 0.07_wp, merr)
         if (allocated(merr)) then
            call test_failed(error, merr%message)
            return
         end if
         call make_trafo(error, rgrid, kgrid, trafo)
         if (allocated(error)) return
         call fill_fields(n, 2, 3, f)
         allocate (g(n, 2), back(n, 2))

         call trafo%fbt_r2k(f(:, 1), g(:, 1), merr)
         if (.not. allocated(merr)) call trafo%fbt_k2r(g(:, 1), back(:, 1), merr)
         if (.not. allocated(merr)) call trafo%fbt_k2r_adj(f(:, 2), g(:, 2), merr)
         if (.not. allocated(merr)) call trafo%fbt_r2k_adj(g(:, 2), back(:, 2), merr)
         if (allocated(merr)) then
            call test_failed(error, merr%message)
            return
         end if
         call check(error, rel_dev(back(:, 1), f(:, 1)) <= 1.0e-13_wp, &
            & "DST-IV k2r(r2k(f)) is not the identity")
         if (allocated(error)) return
         call check(error, rel_dev(back(:, 2), f(:, 2)) <= 1.0e-13_wp, &
            & "DST-IV r2k_adj(k2r_adj(g)) is not the identity")
         if (allocated(error)) return

         call trafo%fbt_r2k_all(f, g, merr)
         if (.not. allocated(merr)) call trafo%fbt_k2r_all(g, back, merr)
         if (allocated(merr)) then
            call test_failed(error, merr%message)
            return
         end if
         call check(error, rel_dev2(back, f) <= 1.0e-13_wp, &
            & "DST-IV batched k2r(r2k(f)) is not the identity")
         if (allocated(error)) return
         deallocate (f, g, back)
      end do
   end subroutine test_dst4_round_trip

   ! --------------------------------------------------------------------------
   ! General quadrature
   ! --------------------------------------------------------------------------

   !> Window errors of both quadrature transforms of a Gaussian at n nodes
   !>
   !> Chebyshev-II + Becke pair, p_r = 1, p_k = 1.5, exponent 0.7; forward
   !> error relative to F(0), backward error absolute (f(0) = 1), over the
   !> nodes with r, k <= 10
   !>
   !> @param[out] error  Test failure
   !> @param[in]  n      Number of r and k nodes
   !> @param[out] err_k  Forward window error
   !> @param[out] err_r  Backward window error
   subroutine gaussian_window_errors(error, n, err_k, err_r)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Number of r and k nodes
      integer, intent(in) :: n
      !> Forward window error
      real(wp), intent(out) :: err_k
      !> Backward window error
      real(wp), intent(out) :: err_r

      real(wp), parameter :: expo = 0.7_wp
      real(wp), parameter :: window = 10.0_wp
      type(moist_math_grid_radial_type) :: rgrid, kgrid
      class(moist_math_grid_radial_trafo_type), allocatable :: trafo
      type(mctc_error), allocatable :: merr
      real(wp), allocatable :: f(:), fk(:), got_k(:), got_r(:)
      real(wp) :: peak

      err_k = huge(1.0_wp)
      err_r = huge(1.0_wp)
      call make_cheb_pair(rgrid, kgrid, n, 1.0_wp, n, 1.5_wp, merr)
      if (allocated(merr)) then
         call test_failed(error, merr%message)
         return
      end if
      call make_trafo(error, rgrid, kgrid, trafo)
      if (allocated(error)) return

      peak = (pi/expo)**1.5_wp
      f = exp(-expo*rgrid%r**2)
      fk = peak*exp(-kgrid%r**2/(4.0_wp*expo))
      allocate (got_k(n), got_r(n))
      call trafo%fbt_r2k(f, got_k, merr)
      if (.not. allocated(merr)) call trafo%fbt_k2r(fk, got_r, merr)
      if (allocated(merr)) then
         call test_failed(error, merr%message)
         return
      end if
      if (all(ieee_is_finite(got_k)) .and. all(ieee_is_finite(got_r)) .and. &
         & any(kgrid%r <= window) .and. any(rgrid%r <= window)) then
         err_k = maxval(abs(got_k - fk), mask=kgrid%r <= window)/peak
         err_r = maxval(abs(got_r - f), mask=rgrid%r <= window)
      end if
   end subroutine gaussian_window_errors

   !> Quadrature transforms of a Gaussian converge in the node count for r, k <= 10
   !>
   !> - Chebyshev-II + Becke, p_r = 1, p_k = 1.5, exponent 0.7
   !> - Measured window errors at 64, 80, 128, 200 nodes: forward (relative
   !>   to F(0)) 1.7e-7, 9.5e-10, 5.4e-14, 3.7e-15; backward (absolute)
   !>   1.7e-5, 9.2e-7, 9.0e-12, 3.9e-15; bounds carry about a factor 3
   !> - The error must also shrink as the node count grows
   !> - Over all nodes the error stalls near 1e-3 (outer k nodes unresolved)
   !>
   !> @param[out] error  Test failure
   subroutine test_quadrature_convergence(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: sizes(4) = [64, 80, 128, 200]
      real(wp), parameter :: bounds_k(4) = [5.0e-7_wp, 3.0e-9_wp, 2.0e-13_wp, 1.5e-14_wp]
      real(wp), parameter :: bounds_r(4) = [5.0e-5_wp, 3.0e-6_wp, 3.0e-11_wp, 1.5e-14_wp]
      real(wp) :: err_k(4), err_r(4)
      integer :: is

      do is = 1, size(sizes)
         call gaussian_window_errors(error, sizes(is), err_k(is), err_r(is))
         if (allocated(error)) return
         call check(error, err_k(is) <= bounds_k(is), &
            & "Quadrature fbt_r2k window error above its bound")
         if (allocated(error)) return
         call check(error, err_r(is) <= bounds_r(is), &
            & "Quadrature fbt_k2r window error above its bound")
         if (allocated(error)) return
      end do
      call check(error, all(err_k(2:) < err_k(:3)) .and. all(err_r(2:) < err_r(:3)), &
         & "Quadrature window error does not decrease with the node count")
   end subroutine test_quadrature_convergence

   !> A node at k = 0 or r = 0 uses sinc(0) = 1 and gives the plain moments
   !>
   !> F(0) = 4*pi sum_i f_i r_i^2 w_i and f(0) = 1/(2*pi^2) sum_j F_j k_j^2 w_j;
   !> hand-made grids, since no rule places a node at 0
   !>
   !> @param[out] error  Test failure
   subroutine test_quadrature_sinc_zero(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      type(moist_math_grid_radial_type) :: rgrid, kgrid
      type(moist_math_grid_radial_trafo_quadrature_type) :: trafo
      type(mctc_error), allocatable :: merr
      real(wp) :: f(4), fk(3), got_k(3), got_r(4), ref

      rgrid%npts = 4
      rgrid%npts_requested = 4
      rgrid%r = [0.0_wp, 0.5_wp, 1.25_wp, 2.0_wp]
      rgrid%w = [0.25_wp, 0.5_wp, 0.625_wp, 0.4_wp]
      rgrid%transform = transform_quadrature
      kgrid%npts = 3
      kgrid%npts_requested = 3
      kgrid%r = [0.0_wp, 0.75_wp, 1.5_wp]
      kgrid%w = [0.375_wp, 0.75_wp, 0.5_wp]
      kgrid%transform = transform_quadrature
      f = [1.0_wp, 0.8_wp, 0.3_wp, 0.1_wp]
      fk = [2.0_wp, 1.0_wp, 0.25_wp]

      call new_quadrature_trafo(trafo, rgrid, kgrid, merr)
      if (.not. allocated(merr)) call trafo%fbt_r2k(f, got_k, merr)
      if (.not. allocated(merr)) call trafo%fbt_k2r(fk, got_r, merr)
      if (allocated(merr)) then
         call test_failed(error, merr%message)
         return
      end if
      call check(error, all(ieee_is_finite(got_k)) .and. all(ieee_is_finite(got_r)), &
         & "Quadrature trafo with a node at zero must stay finite")
      if (allocated(error)) return

      ref = 4.0_wp*pi*sum(f*rgrid%r**2*rgrid%w)
      call check(error, abs(got_k(1) - ref) <= 16.0_wp*epsilon(1.0_wp)*abs(ref), &
         & "Quadrature F(k = 0) must be the plain volume moment")
      if (allocated(error)) return
      ref = sum(fk*kgrid%r**2*kgrid%w)/(2.0_wp*pi*pi)
      call check(error, abs(got_r(1) - ref) <= 16.0_wp*epsilon(1.0_wp)*abs(ref), &
         & "Quadrature f(r = 0) must be the plain k-space moment")
   end subroutine test_quadrature_sinc_zero

   !> The dense kernel on a uniform pair equals the DST-IV transform to round-off
   !>
   !> On the midpoint nodes both evaluate the same discrete operator, so this
   !> pins the quadrature prefactors against the independent DST-IV path;
   !> the quadrature trafo is built directly, since the factory never selects
   !> it for a DST-IV pair; deviation relative to the largest output element,
   !> measured at most 9.7e-15 (the dense path reduces sin arguments up to
   !> about 300 in floating point, the DST-IV phase is exact)
   !>
   !> @param[out] error  Test failure
   subroutine test_quadrature_vs_dst4(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: n = 96
      type(moist_math_grid_radial_type) :: rgrid, kgrid
      class(moist_math_grid_radial_trafo_type), allocatable :: dst4
      type(moist_math_grid_radial_trafo_quadrature_type) :: quad
      type(mctc_error), allocatable :: merr
      real(wp), allocatable :: f(:, :)
      real(wp) :: got(n), ref(n)
      integer :: op

      call new_uniform_radial_pair(rgrid, kgrid, n, 0.1_wp, merr)
      if (.not. allocated(merr)) call new_quadrature_trafo(quad, rgrid, kgrid, merr)
      if (allocated(merr)) then
         call test_failed(error, merr%message)
         return
      end if
      call make_trafo(error, rgrid, kgrid, dst4)
      if (allocated(error)) return
      call fill_fields(n, 1, 5, f)

      do op = op_r2k, op_k2r_adj
         call apply_one(dst4, op, f(:, 1), ref, merr)
         if (.not. allocated(merr)) call apply_one(quad, op, f(:, 1), got, merr)
         if (allocated(merr)) then
            call test_failed(error, merr%message)
            return
         end if
         call check(error, rel_dev(got, ref) <= 1.0e-13_wp, &
            & "Quadrature trafo on a uniform pair deviates from DST-IV")
         if (allocated(error)) return
      end do
   end subroutine test_quadrature_vs_dst4

   ! --------------------------------------------------------------------------
   ! Contracts shared by both implementations
   ! --------------------------------------------------------------------------

   !> Both adjoints are the transposes: <T a, b> = <a, T^T b>
   !>
   !> The residual is scaled by ||T a|| ||b||, the natural size of the
   !> bilinear form; DST-IV on 96 nodes, quadrature on 64 r and 80 k nodes
   !>
   !> @param[out] error  Test failure
   subroutine test_adjoint(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: kinds(2) = [transform_dst4, transform_quadrature]
      integer, parameter :: fwd_ops(2) = [op_r2k, op_k2r]
      integer, parameter :: adj_ops(2) = [op_r2k_adj, op_k2r_adj]
      type(moist_math_grid_radial_type) :: rgrid, kgrid
      class(moist_math_grid_radial_trafo_type), allocatable :: trafo
      type(mctc_error), allocatable :: merr
      real(wp), allocatable :: a(:, :), b(:, :), ta(:), tb(:)
      real(wp) :: lhs, rhs, scal
      integer :: ik, io, n_in, n_out

      do ik = 1, size(kinds)
         call make_case(error, kinds(ik), rgrid, kgrid, trafo)
         if (allocated(error)) return
         do io = 1, size(fwd_ops)
            call op_sizes(trafo, fwd_ops(io), n_in, n_out)
            call fill_fields(n_in, 1, 0, a)
            call fill_fields(n_out, 1, 41, b)
            allocate (ta(n_out), tb(n_in))
            call apply_one(trafo, fwd_ops(io), a(:, 1), ta, merr)
            if (.not. allocated(merr)) call apply_one(trafo, adj_ops(io), b(:, 1), tb, merr)
            if (allocated(merr)) then
               call test_failed(error, merr%message)
               return
            end if
            lhs = dot_product(ta, b(:, 1))
            rhs = dot_product(a(:, 1), tb)
            scal = max(norm2(ta)*norm2(b(:, 1)), 1.0e-30_wp)
            call check(error, abs(lhs - rhs)/scal <= 1.0e-14_wp, &
               & "Adjoint dot-product identity failed")
            if (allocated(error)) return
            deallocate (ta, tb)
         end do
      end do
   end subroutine test_adjoint

   !> Batched transforms reproduce the scalar loop for all four directions
   !>
   !> @param[out] error  Test failure
   subroutine test_batched(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: kinds(2) = [transform_dst4, transform_quadrature]
      integer, parameter :: nb = 4
      type(moist_math_grid_radial_type) :: rgrid, kgrid
      class(moist_math_grid_radial_trafo_type), allocatable :: trafo
      type(mctc_error), allocatable :: merr
      real(wp), allocatable :: f(:, :), loop(:, :), batch(:, :)
      integer :: ik, op, j, n_in, n_out

      do ik = 1, size(kinds)
         call make_case(error, kinds(ik), rgrid, kgrid, trafo)
         if (allocated(error)) return
         do op = op_r2k, op_k2r_adj
            call op_sizes(trafo, op, n_in, n_out)
            call fill_fields(n_in, nb, 7*op, f)
            allocate (loop(n_out, nb), batch(n_out, nb))
            do j = 1, nb
               call apply_one(trafo, op, f(:, j), loop(:, j), merr)
               if (allocated(merr)) then
                  call test_failed(error, merr%message)
                  return
               end if
            end do
            call apply_all(trafo, op, f, batch, merr)
            if (allocated(merr)) then
               call test_failed(error, merr%message)
               return
            end if
            call check(error, rel_dev2(batch, loop) <= 1.0e-14_wp, &
               & "Batched transform disagrees with the scalar loop")
            if (allocated(error)) return
            deallocate (loop, batch)
         end do
      end do
   end subroutine test_batched

   !> A zero-width batch on a fresh trafo returns without an error
   !>
   !> @param[out] error  Test failure
   subroutine test_empty_batch(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: kinds(2) = [transform_dst4, transform_quadrature]
      type(moist_math_grid_radial_type) :: rgrid, kgrid
      class(moist_math_grid_radial_trafo_type), allocatable :: trafo
      type(mctc_error), allocatable :: merr
      real(wp), allocatable :: f(:, :), g(:, :)
      integer :: ik, op, n_in, n_out

      do ik = 1, size(kinds)
         call make_case(error, kinds(ik), rgrid, kgrid, trafo)
         if (allocated(error)) return
         do op = op_r2k, op_k2r_adj
            call op_sizes(trafo, op, n_in, n_out)
            allocate (f(n_in, 0), g(n_out, 0))
            call apply_all(trafo, op, f, g, merr)
            call check(error, .not. allocated(merr), "Empty batch reported an error")
            if (allocated(error)) return
            deallocate (f, g)
         end do
      end do
   end subroutine test_empty_batch

   !> Wrong input or output lengths and mismatched batch widths are refused
   !>
   !> @param[out] error  Test failure
   subroutine test_size_errors(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: kinds(2) = [transform_dst4, transform_quadrature]
      type(moist_math_grid_radial_type) :: rgrid, kgrid
      class(moist_math_grid_radial_trafo_type), allocatable :: trafo
      type(mctc_error), allocatable :: merr
      real(wp), allocatable :: f(:), g(:), f2(:, :), g2(:, :)
      integer :: ik, op, n_in, n_out

      do ik = 1, size(kinds)
         call make_case(error, kinds(ik), rgrid, kgrid, trafo)
         if (allocated(error)) return
         do op = op_r2k, op_k2r_adj
            call op_sizes(trafo, op, n_in, n_out)
            allocate (f(n_in + 1), g(n_out))
            call apply_one(trafo, op, f, g, merr)
            call check(error, allocated(merr), "Scalar transform accepted a long input")
            if (allocated(error)) return
            deallocate (f, g)
            allocate (f(n_in), g(n_out - 1))
            f(:) = 0.0_wp
            call apply_one(trafo, op, f, g, merr)
            call check(error, allocated(merr), "Scalar transform accepted a short output")
            if (allocated(error)) return
            deallocate (f, g)
            allocate (f2(n_in, 2), g2(n_out, 3))
            f2(:, :) = 0.0_wp
            call apply_all(trafo, op, f2, g2, merr)
            call check(error, allocated(merr), "Batched transform accepted mismatched widths")
            if (allocated(error)) return
            deallocate (f2, g2)
            allocate (f2(n_in - 1, 2), g2(n_out, 2))
            f2(:, :) = 0.0_wp
            call apply_all(trafo, op, f2, g2, merr)
            call check(error, allocated(merr), "Batched transform accepted a short input column")
            if (allocated(error)) return
            deallocate (f2, g2)
         end do
      end do
   end subroutine test_size_errors

   !> The base default batched loops apply the scalar form per column and check shapes
   !>
   !> @param[out] error  Test failure
   subroutine test_base_default(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      real(wp), parameter :: factors(4) = [2.0_wp, 3.0_wp, 5.0_wp, 7.0_wp]
      type(scale_trafo_type) :: trafo
      type(mctc_error), allocatable :: merr
      real(wp), allocatable :: f(:, :), g(:, :)
      integer :: op

      trafo%nr = 5
      trafo%nk = 5
      call check(error, trafo%implementation() == 0, "Scale trafo tag")
      if (allocated(error)) return
      call fill_fields(5, 3, 0, f)
      allocate (g(5, 3))
      do op = op_r2k, op_k2r_adj
         call apply_all(trafo, op, f, g, merr)
         if (allocated(merr)) then
            call test_failed(error, merr%message)
            return
         end if
         call check(error, all(g == factors(op)*f), &
            & "Default batched loop does not apply the scalar form per column")
         if (allocated(error)) return
      end do
      deallocate (g)
      allocate (g(5, 2))
      call apply_all(trafo, op_r2k, f, g, merr)
      call check(error, allocated(merr), "Default batched loop accepted mismatched widths")
      if (allocated(error)) return
      deallocate (g)
      allocate (g(4, 3))
      call apply_all(trafo, op_k2r, f, g, merr)
      call check(error, allocated(merr), "Default batched loop accepted a short output column")
   end subroutine test_base_default

   !> A clone is independent of its template, and a trafo of its grids
   !>
   !> - Template results recorded, grids overwritten afterwards, template
   !>   cloned with allocate(source=) and deallocated
   !> - The clone reproduces the template bit for bit, scalar and batched,
   !>   first with a batch width other than the one the template cached last
   !>
   !> @param[out] error  Test failure
   subroutine test_clone(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      integer, parameter :: kinds(2) = [transform_dst4, transform_quadrature]
      type(moist_math_grid_radial_type) :: rgrid, kgrid
      class(moist_math_grid_radial_trafo_type), allocatable :: template, clone
      type(mctc_error), allocatable :: merr
      real(wp), allocatable :: f3(:, :), f5(:, :), ref1(:), ref3(:, :), ref5(:, :)
      real(wp), allocatable :: got1(:), got3(:, :), got5(:, :)
      integer :: ik

      do ik = 1, size(kinds)
         call make_case(error, kinds(ik), rgrid, kgrid, template)
         if (allocated(error)) return
         call fill_fields(template%nr, 3, 1, f3)
         call fill_fields(template%nr, 5, 2, f5)
         allocate (ref1(template%nk), ref3(template%nk, 3), ref5(template%nk, 5))
         allocate (got1(template%nk), got3(template%nk, 3), got5(template%nk, 5))

         call template%fbt_r2k(f3(:, 1), ref1, merr)
         if (.not. allocated(merr)) call template%fbt_r2k_all(f5, ref5, merr)
         if (.not. allocated(merr)) call template%fbt_r2k_all(f3, ref3, merr)
         if (allocated(merr)) then
            call test_failed(error, merr%message)
            return
         end if

         ! The trafo copied what it needs; the grids may change or vanish
         rgrid%r(:) = -1.0_wp
         kgrid%r(:) = -1.0_wp
         deallocate (rgrid%r, rgrid%w, kgrid%r, kgrid%w)

         allocate (clone, source=template)
         deallocate (template)

         call clone%fbt_r2k_all(f5, got5, merr)
         if (.not. allocated(merr)) call clone%fbt_r2k(f3(:, 1), got1, merr)
         if (.not. allocated(merr)) call clone%fbt_r2k_all(f3, got3, merr)
         if (allocated(merr)) then
            call test_failed(error, merr%message)
            return
         end if
         call check(error, all(got1 == ref1) .and. all(got3 == ref3) .and. all(got5 == ref5), &
            & "Clone does not reproduce its template")
         if (allocated(error)) return
         deallocate (clone, f3, f5, ref1, ref3, ref5, got1, got3, got5)
      end do
   end subroutine test_clone

   ! --------------------------------------------------------------------------
   ! Minimal scale trafo for the base defaults
   ! --------------------------------------------------------------------------

   !> No implementation tag
   !>
   !> @param[in] self  Trafo instance
   pure function scale_implementation(self) result(tag)
      !> Trafo instance
      class(scale_trafo_type), intent(in) :: self
      !> 0
      integer :: tag

      tag = 0*self%nr
   end function scale_implementation

   !> f_out = 2*f_in
   !>
   !> @param[in,out] self   Trafo instance
   !> @param[in]     f_in   Input field
   !> @param[out]    f_out  Output field
   !> @param[out]    error  Never set
   subroutine scale_r2k(self, f_in, f_out, error)
      !> Trafo instance
      class(scale_trafo_type), intent(inout) :: self
      !> Input field
      real(wp), intent(in) :: f_in(:)
      !> Output field
      real(wp), intent(out) :: f_out(:)
      !> Error handling
      type(mctc_error), allocatable, intent(out) :: error

      f_out = 2.0_wp*f_in
   end subroutine scale_r2k

   !> f_out = 3*f_in
   !>
   !> @param[in,out] self   Trafo instance
   !> @param[in]     f_in   Input field
   !> @param[out]    f_out  Output field
   !> @param[out]    error  Never set
   subroutine scale_k2r(self, f_in, f_out, error)
      !> Trafo instance
      class(scale_trafo_type), intent(inout) :: self
      !> Input field
      real(wp), intent(in) :: f_in(:)
      !> Output field
      real(wp), intent(out) :: f_out(:)
      !> Error handling
      type(mctc_error), allocatable, intent(out) :: error

      f_out = 3.0_wp*f_in
   end subroutine scale_k2r

   !> f_out = 5*f_in
   !>
   !> @param[in,out] self   Trafo instance
   !> @param[in]     f_in   Input field
   !> @param[out]    f_out  Output field
   !> @param[out]    error  Never set
   subroutine scale_r2k_adj(self, f_in, f_out, error)
      !> Trafo instance
      class(scale_trafo_type), intent(inout) :: self
      !> Input field
      real(wp), intent(in) :: f_in(:)
      !> Output field
      real(wp), intent(out) :: f_out(:)
      !> Error handling
      type(mctc_error), allocatable, intent(out) :: error

      f_out = 5.0_wp*f_in
   end subroutine scale_r2k_adj

   !> f_out = 7*f_in
   !>
   !> @param[in,out] self   Trafo instance
   !> @param[in]     f_in   Input field
   !> @param[out]    f_out  Output field
   !> @param[out]    error  Never set
   subroutine scale_k2r_adj(self, f_in, f_out, error)
      !> Trafo instance
      class(scale_trafo_type), intent(inout) :: self
      !> Input field
      real(wp), intent(in) :: f_in(:)
      !> Output field
      real(wp), intent(out) :: f_out(:)
      !> Error handling
      type(mctc_error), allocatable, intent(out) :: error

      f_out = 7.0_wp*f_in
   end subroutine scale_k2r_adj

end module test_math_grid_radial_trafo
