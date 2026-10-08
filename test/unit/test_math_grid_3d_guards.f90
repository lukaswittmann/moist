!> Guards of the 3D volume grids whose absence lets a run continue with wrong numbers
!>
!> - Reciprocal period: a molecular update that outgrows the fixed k-grid period fails,
!>   instead of aliasing the field into a period the k-grid no longer resolves
!> - Cartesian box fit: a solute plus margin wider than the box fails, instead of
!>   wrapping density across the periodic boundary
!> - Transform blocks: mis-shaped real- or reciprocal-space blocks fail before the
!>   FFT or FINUFFT backend reads or writes past their extents
!> - Molecular transform state: an unprepared engine, or one prepared on an outdated
!>   grid geometry, fails instead of transforming with missing or stale plans
!> - Molecular plan lifecycle: re-creation and copies come out unprepared, so two
!>   engines never share or reuse FINUFFT plans
!> - Volume adjoints: mis-sized channels fail instead of being summed out of bounds
!> - Width adjoints: point grids and molecular Hessians reject Gaussian-width adjoints
!>   they cannot contract, instead of dropping that contribution
!> - Partition curvature: Becke stiffness one at a nucleus and coincident elliptic
!>   centers fail instead of returning an undefined Hessian
!> - Power-Voronoi radii: a radii list of the wrong length fails instead of being read
!>   out of bounds
module test_math_grid_3d_guards
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use mctc_env, only: wp
   use mctc_env_error, only: mctc_error => error_type
   use mctc_io, only: structure_type, new
   use mstore, only: get_structure
   use testdrive, only: new_unittest, unittest_type, error_type, check, test_failed
   use test_helpers, only: get_cartesian_gaussian_grid, center_at_origin, get_uniform_recipe, &
      & get_qc_handymod_recipe, check_moist_error
   use moist_math_grid_3d_base, only: moist_math_grid_3d_trafo_type
   use moist_math_grid_3d_adjoint, only: volume_adjoint_type
   use moist_math_grid_3d_cartesian, only: moist_math_grid_3d_cartesian_type, new_cartesian_point_grid
   use moist_math_grid_atomic_recipe, only: moist_math_grid_atomic_recipe_type, &
      & moist_math_grid_atomic_recipe_override_type
   use moist_math_grid_3d_kernel_becke, only: becke_partition_type
   use moist_math_grid_3d_kernel_pvoronoi, only: pvoronoi_partition_type
   use moist_math_grid_3d_molecular, only: moist_math_grid_3d_molecular_type, &
      & new_molecular_point_grid, new_molecular_gaussian_grid, molecular_grid_set_kgrid, &
      & moist_math_grid_3d_molecular_trafo_type, new_molecular_grid_trafo
   implicit none(type, external)
   private

   public :: collect_math_grid_3d_guards

   !> Reciprocal resolution (bohr) of the NUFFT probe grid, k_max = pi/dr
   real(wp), parameter :: probe_dr = 0.70_wp
   !> Outer radial clamp of the probe grid (bohr)
   real(wp), parameter :: probe_rmax = 16.0_wp
   !> Radial shells per atom of the probe grid
   integer, parameter :: probe_nrad = 50
   !> Raw Lebedev point count per shell of the probe grid
   integer, parameter :: probe_nang = 110

contains

   !> Collect all math_grid_3d_guards tests
   !>
   !> @param[out] testsuite  collected unit tests
   subroutine collect_math_grid_3d_guards(testsuite)
      !> Collected unit tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("molecular_motion_period", test_molecular_motion_period, should_fail=.true.), &
                  new_unittest("cartesian_bad_box_fit", test_cartesian_fit_fails, should_fail=.true.), &
                  new_unittest("cartesian_trafo_block_guards", test_cartesian_block_guards), &
                  new_unittest("molecular_trafo_block_guards", test_molecular_block_guards), &
                  new_unittest("molecular_trafo_bad_unprepared_forward", test_trafo_unprepared_r2k_fails, &
                     & should_fail=.true.), &
                  new_unittest("molecular_trafo_bad_unprepared_backward", test_trafo_unprepared_k2r_fails, &
                     & should_fail=.true.), &
                  new_unittest("molecular_trafo_lifecycle", test_molecular_trafo_lifecycle), &
                  new_unittest("nufft_trafo_stale_after_update", test_trafo_stale_after_update, &
                     & should_fail=.true.), &
                  new_unittest("nufft_trafo_stale_after_reconstruction", test_trafo_stale_after_reconstruction, &
                     & should_fail=.true.), &
                  new_unittest("volume_adjoint_sizes", test_volume_adjoint_sizes), &
                  new_unittest("point_grid_rejects_width_adjoint", test_point_grid_width_adjoint_fails, &
                     & should_fail=.true.), &
                  new_unittest("hessian_rejects_width_adjoint", test_hessian_width_adjoint_fails, &
                     & should_fail=.true.), &
                  new_unittest("hessian_rejects_becke1_at_nucleus", test_hessian_becke1_nucleus_fails, &
                     & should_fail=.true.), &
                  new_unittest("hessian_rejects_coincident_centers", test_hessian_coincident_fails, &
                     & should_fail=.true.), &
                  new_unittest("power_radii_size", test_power_radii_size_fails, should_fail=.true.) &
                  ]
   end subroutine collect_math_grid_3d_guards

   !> Fail an expected-failure test only on the targeted library error
   !>
   !> - Used by `should_fail=.true.` tests, where test-drive inverts the
   !>   verdict: a raised test failure passes, a clean return fails
   !> - Fails the test only when `err` names `expected`; no error or a
   !>   different one returns cleanly, so the test is reported as failed
   !> - Callers must therefore return cleanly (never `test_failed`) when their
   !>   own setup fails, or a broken fixture would pass as the expected error
   !>
   !> @param[out] error     test failure, set only on the expected error
   !> @param[in]  err       library error, possibly unallocated
   !> @param[in]  expected  substring the expected error message must contain
   subroutine expect_error(error, err, expected)
      !> Test failure, set only on the expected error
      type(error_type), allocatable, intent(out) :: error
      !> Library error, possibly unallocated
      type(mctc_error), allocatable, intent(in) :: err
      !> Substring the expected error message must contain
      character(len=*), intent(in) :: expected

      if (.not. allocated(err)) return
      if (index(err%message, expected) > 0) call test_failed(error, err%message)
   end subroutine expect_error

   !> Fail unless `merr` is set and names `expected`
   !>
   !> For an error in the middle of a test that continues afterwards; a test
   !> whose only purpose is one error uses `should_fail` instead
   !>
   !> @param[out] error     test failure
   !> @param[in]  merr      library error, possibly unallocated
   !> @param[in]  expected  substring the message must contain
   !> @param[in]  label     case label for the failure message
   subroutine require_message(error, merr, expected, label)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Library error, possibly unallocated
      type(mctc_error), allocatable, intent(in) :: merr
      !> Substring the message must contain
      character(len=*), intent(in) :: expected
      !> Case label
      character(len=*), intent(in) :: label

      if (.not. allocated(merr)) then
         call test_failed(error, label//": expected an error containing '"//expected//"'")
      else if (index(merr%message, expected) == 0) then
         call test_failed(error, label//": message '"//merr%message//"' lacks '"//expected//"'")
      end if
   end subroutine require_message

   !> Forward transform through the abstract engine, failing the test on error
   !>
   !> @param[out]    error  test error
   !> @param[in,out] trafo  transform engine, prepared where the backend needs it
   !> @param[in,out] f_r    real-space block
   !> @param[out]    f_k    reciprocal-space block
   subroutine forward(error, trafo, f_r, f_k)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      !> Transform engine (prepared where the backend needs it)
      class(moist_math_grid_3d_trafo_type), intent(inout) :: trafo
      !> Real-space block
      real(wp), intent(inout), contiguous, target :: f_r(:, :)
      !> Reciprocal-space block
      complex(wp), intent(out), contiguous, target :: f_k(:, :)

      type(mctc_error), allocatable :: merr

      call trafo%fft_r2k(f_r, f_k, merr)
      if (allocated(merr)) call test_failed(error, merr%message)
   end subroutine forward

   !> Backward transform through the abstract engine, failing the test on error
   !>
   !> @param[out]    error  test error
   !> @param[in,out] trafo  transform engine
   !> @param[in,out] f_k    reciprocal-space block, destroyed
   !> @param[out]    f_r    real-space block
   subroutine backward(error, trafo, f_k, f_r)
      !> Test error
      type(error_type), allocatable, intent(out) :: error
      !> Transform engine
      class(moist_math_grid_3d_trafo_type), intent(inout) :: trafo
      !> Reciprocal-space block (destroyed)
      complex(wp), intent(inout), contiguous, target :: f_k(:, :)
      !> Real-space block
      real(wp), intent(out), contiguous, target :: f_r(:, :)

      type(mctc_error), allocatable :: merr

      call trafo%fft_k2r(f_k, f_r, merr)
      if (allocated(merr)) call test_failed(error, merr%message)
   end subroutine backward

   !> Check that a complex value has finite real and imaginary parts
   !>
   !> @param[in] z  value to test
   pure function complex_is_finite(z) result(ok)
      !> Value to test
      complex(wp), intent(in) :: z
      !> True if both the real and imaginary parts are finite
      logical :: ok

      ok = ieee_is_finite(real(z, wp)) .and. ieee_is_finite(aimag(z))
   end function complex_is_finite

   !> Three heteronuclear sites away from membership crossings
   !>
   !> @param[out] mol  test molecule
   subroutine make_molecule(mol)
      !> Test molecule
      type(structure_type), intent(out) :: mol

      call new(mol, [8, 1, 6], reshape([-0.6_wp, -0.2_wp, 0.1_wp, &
         & 0.7_wp, 0.3_wp, -0.25_wp, 0.1_wp, 1.2_wp, 0.4_wp], [3, 3]))
   end subroutine make_molecule

   !> A stretch beyond the fixed reciprocal period fails the molecular update
   subroutine test_molecular_motion_period(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(mctc_error), allocatable :: err
      type(moist_math_grid_3d_molecular_type) :: domain
      type(moist_math_grid_atomic_recipe_type) :: recipe
      type(structure_type) :: mol
      !> Outer HandyMod radius (bohr)
      real(wp), parameter :: rmax = 4.0_wp

      call get_qc_handymod_recipe(recipe, err, nrad=12, rmax=rmax)
      if (.not. allocated(err)) call new_molecular_point_grid(domain, err, recipe=recipe, dr=1.0_wp, &
         & kbuffer=1.0_wp)
      if (allocated(err)) return
      call new(mol, [1, 1], reshape([-0.7_wp, 0.0_wp, 0.0_wp, 0.7_wp, 0.0_wp, 0.0_wp], [3, 2]))
      call domain%update(mol, err)
      if (allocated(err)) return
      mol%xyz(1, 2) = mol%xyz(1, 2) + 8.0_wp
      call domain%update(mol, err)
      call expect_error(error, err, "period exceeded")
   end subroutine test_molecular_motion_period

   !> Solute wider than the box minus its margin fails the update and names the needed length
   subroutine test_cartesian_fit_fails(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(moist_math_grid_3d_cartesian_type) :: grid
      type(mctc_error), allocatable :: err
      type(structure_type) :: mol

      ! Centroid 0.2, reach 3.2: needs 2*(3.2 + 2) bohr along x, the box has 10
      call new(mol, [1, 1], reshape([-3.0_wp, 0.0_wp, 0.0_wp, 3.4_wp, 0.0_wp, 0.0_wp], [3, 2]))
      call new_cartesian_point_grid(grid, err, nx=20, ny=20, nz=20, dr=0.5_wp, margin=2.0_wp)
      if (allocated(err)) return
      call grid%update(mol, err)
      call expect_error(error, err, "do not fit the box along x; nx*dr must be at least 10.400 bohr")
   end subroutine test_cartesian_fit_fails

   !> Build an H2 molecular grid with its k-grid and an unprepared engine bound to it
   !>
   !> @param[out] mgrid   molecular grid, must outlive mtrafo
   !> @param[out] mtrafo  unprepared engine bound to mgrid
   !> @param[out] err     library error from the construction
   subroutine molecular_trafo_fixture(mgrid, mtrafo, err)
      !> Molecular grid, must outlive mtrafo
      type(moist_math_grid_3d_molecular_type), intent(inout), target :: mgrid
      !> Unprepared engine bound to mgrid
      type(moist_math_grid_3d_molecular_trafo_type), intent(inout) :: mtrafo
      !> Library error from the construction
      type(mctc_error), allocatable, intent(out) :: err

      type(structure_type) :: mol
      type(moist_math_grid_atomic_recipe_type) :: recipe
      type(moist_math_grid_atomic_recipe_override_type), allocatable :: overrides(:)

      call get_structure(mol, "MB16-43", "H2")
      call center_at_origin(mol)
      call get_uniform_recipe(recipe, overrides, 12, 26, err, rmax=4.0_wp)
      if (allocated(err)) return
      call new_molecular_point_grid(mgrid, err, recipe=recipe, overrides=overrides, reciprocal=.false.)
      if (allocated(err)) return
      call mgrid%update(mol, err)
      if (allocated(err)) return
      call molecular_grid_set_kgrid(mgrid, 0.8_wp, err)
      if (allocated(err)) return
      call new_molecular_grid_trafo(mtrafo, mgrid)
   end subroutine molecular_trafo_fixture

   !> Cartesian engine rejects every mis-shaped block
   subroutine test_cartesian_block_guards(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(moist_math_grid_3d_cartesian_type), target :: grid
      class(moist_math_grid_3d_trafo_type), allocatable :: trafo
      type(mctc_error), allocatable :: merr

      call get_cartesian_gaussian_grid(grid, 12, 11, 10, 0.4_wp, error=merr)
      if (.not. allocated(merr)) call grid%new_trafo(trafo, merr)
      call check_moist_error(error, merr)
      if (allocated(error)) return
      call check_block_guards(error, trafo, grid%ngrid, grid%npts_k)
      call trafo%destroy()
   end subroutine test_cartesian_block_guards

   !> Prepared molecular engine rejects every mis-shaped block
   subroutine test_molecular_block_guards(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(moist_math_grid_3d_molecular_type), target :: mgrid
      type(moist_math_grid_3d_molecular_trafo_type) :: mtrafo
      type(mctc_error), allocatable :: merr

      call molecular_trafo_fixture(mgrid, mtrafo, merr)
      if (.not. allocated(merr)) call mtrafo%prepare(1, merr)
      call check_moist_error(error, merr)
      if (allocated(error)) return
      call check_block_guards(error, mtrafo, mgrid%ngrid, mgrid%npts_k)
      call mtrafo%destroy()
   end subroutine test_molecular_block_guards

   !> Run the mis-shaped block table through both directions of `trafo`
   !>
   !> - Unequal input and output widths, a wrong leading extent on either
   !>   block, and an empty batch are each rejected by name
   !> - A well-shaped single-column block then transforms cleanly
   !>
   !> @param[out]    error   test failure
   !> @param[in,out] trafo   engine bound to a grid, prepared for one column
   !> @param[in]     ngrid   real-space point count of the grid
   !> @param[in]     npts_k  reciprocal-space point count of the grid
   subroutine check_block_guards(error, trafo, ngrid, npts_k)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Engine bound to a grid, prepared for one column
      class(moist_math_grid_3d_trafo_type), intent(inout) :: trafo
      !> Real-space point count
      integer, intent(in) :: ngrid
      !> Reciprocal-space point count
      integer, intent(in) :: npts_k

      integer, parameter :: ncase = 5
      character(len=*), parameter :: width = "differ in batch width", extent = "field block size", &
         & empty = "at least one column"
      !> Real-space rows and columns, reciprocal-space rows and columns per case
      integer :: shapes(4, ncase)
      character(len=len(width)) :: expected(ncase)
      type(mctc_error), allocatable :: merr
      real(wp), allocatable :: f(:, :)
      complex(wp), allocatable :: fk(:, :)
      character(len=16) :: label
      integer :: ic

      shapes(:, 1) = [ngrid, 2, npts_k, 1]
      shapes(:, 2) = [ngrid, 1, npts_k, 2]
      shapes(:, 3) = [ngrid + 1, 1, npts_k, 1]
      shapes(:, 4) = [ngrid, 1, npts_k - 1, 1]
      shapes(:, 5) = [ngrid, 0, npts_k, 0]
      expected = [character(len=len(width)) :: width, width, extent, extent, empty]
      do ic = 1, ncase
         allocate (f(shapes(1, ic), shapes(2, ic)), fk(shapes(3, ic), shapes(4, ic)))
         write (label, "(a, i0)") "case ", ic
         f = 1.0_wp
         call trafo%fft_r2k(f, fk, merr)
         call require_message(error, merr, trim(expected(ic)), "forward "//trim(label))
         if (allocated(error)) return
         fk = (0.0_wp, 0.0_wp)
         call trafo%fft_k2r(fk, f, merr)
         call require_message(error, merr, trim(expected(ic)), "backward "//trim(label))
         if (allocated(error)) return
         deallocate (f, fk)
      end do

      allocate (f(ngrid, 1), fk(npts_k, 1))
      f = 1.0_wp
      call forward(error, trafo, f, fk)
      if (.not. allocated(error)) call backward(error, trafo, fk, f)
   end subroutine check_block_guards

   !> Forward transform on an unprepared molecular engine fails
   subroutine test_trafo_unprepared_r2k_fails(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(moist_math_grid_3d_molecular_type), target :: mgrid
      type(moist_math_grid_3d_molecular_trafo_type) :: mtrafo
      type(mctc_error), allocatable :: err
      real(wp), allocatable :: f(:, :)
      complex(wp), allocatable :: fk(:, :)

      call molecular_trafo_fixture(mgrid, mtrafo, err)
      if (allocated(err)) return
      allocate (f(mgrid%ngrid, 1), fk(mgrid%npts_k, 1))
      f = 1.0_wp
      call mtrafo%fft_r2k(f, fk, err)
      call expect_error(error, err, "call prepare(nv)")
      call mtrafo%destroy()
   end subroutine test_trafo_unprepared_r2k_fails

   !> Backward transform on an unprepared molecular engine fails
   subroutine test_trafo_unprepared_k2r_fails(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(moist_math_grid_3d_molecular_type), target :: mgrid
      type(moist_math_grid_3d_molecular_trafo_type) :: mtrafo
      type(mctc_error), allocatable :: err
      real(wp), allocatable :: f(:, :)
      complex(wp), allocatable :: fk(:, :)

      call molecular_trafo_fixture(mgrid, mtrafo, err)
      if (allocated(err)) return
      allocate (f(mgrid%ngrid, 1), fk(mgrid%npts_k, 1))
      fk = (0.0_wp, 0.0_wp)
      call mtrafo%fft_k2r(fk, f, err)
      call expect_error(error, err, "call prepare(nv)")
      call mtrafo%destroy()
   end subroutine test_trafo_unprepared_k2r_fails

   !> Molecular engine owns its plans through replacement, copies and finalization
   !>
   !> - Re-initializing a prepared engine leaves it unprepared; it prepares again
   !>   and reproduces the earlier spectrum
   !> - Assignment yields an unprepared copy bound to the same grid; preparing and
   !>   finalizing the copy leaves the source's plans intact
   !> - The allocatable factory replaces and deallocates a prepared engine cleanly
   subroutine test_molecular_trafo_lifecycle(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(moist_math_grid_3d_molecular_type), target :: mgrid
      type(moist_math_grid_3d_molecular_trafo_type) :: mtrafo
      type(moist_math_grid_3d_molecular_trafo_type), allocatable :: copy
      class(moist_math_grid_3d_trafo_type), allocatable :: owned
      type(mctc_error), allocatable :: merr
      real(wp), allocatable :: f(:, :)
      complex(wp), allocatable :: fk(:, :), ref(:, :)
      real(wp) :: thr
      integer :: i

      call molecular_trafo_fixture(mgrid, mtrafo, merr)
      if (.not. allocated(merr)) call mtrafo%prepare(1, merr)
      call check_moist_error(error, merr)
      if (allocated(error)) return
      allocate (f(mgrid%ngrid, 1), fk(mgrid%npts_k, 1), ref(mgrid%npts_k, 1))
      f(:, 1) = exp(-sum(mgrid%xyz**2, dim=1))
      call forward(error, mtrafo, f, ref)
      if (allocated(error)) return
      thr = 1.0e-12_wp*sum(abs(mgrid%w*f(:, 1)))

      ! Replacement releases the plans and leaves an unprepared engine
      call new_molecular_grid_trafo(mtrafo, mgrid)
      call check(error, mtrafo%get_ntrans() == 0, "replacement must release the plans")
      if (allocated(error)) return
      call mtrafo%fft_r2k(f, fk, merr)
      call require_message(error, merr, "call prepare(nv)", "replaced engine")
      if (allocated(error)) return
      call mtrafo%prepare(1, merr)
      call check_moist_error(error, merr)
      if (allocated(error)) return
      call forward(error, mtrafo, f, fk)
      if (allocated(error)) return
      do i = 1, mgrid%npts_k
         call check(error, complex_is_finite(ref(i, 1)), "reference spectrum is not finite")
         if (allocated(error)) return
         call check(error, fk(i, 1), ref(i, 1), thr=thr, more="re-prepared engine deviates from the first spectrum")
         if (allocated(error)) return
      end do

      ! Copy is unprepared and owns its own plans once prepared
      allocate (copy)
      copy = mtrafo
      call check(error, copy%get_ntrans() == 0 .and. copy%get_nufft_tol() == mtrafo%get_nufft_tol(), &
         & "copy must keep the settings and come out unprepared")
      if (allocated(error)) return
      call copy%fft_r2k(f, fk, merr)
      call require_message(error, merr, "call prepare(nv)", "copied engine")
      if (allocated(error)) return
      call copy%prepare(1, merr)
      call check_moist_error(error, merr)
      if (allocated(error)) return
      ! Same grid binding and transform family: the copy reproduces the source's spectrum
      call forward(error, copy, f, fk)
      if (allocated(error)) return
      do i = 1, mgrid%npts_k
         call check(error, complex_is_finite(ref(i, 1)), "reference spectrum is not finite")
         if (allocated(error)) return
         call check(error, fk(i, 1), ref(i, 1), thr=thr, more="prepared copy deviates from its source")
         if (allocated(error)) return
      end do
      call check(error, copy%uses_type12() .eqv. mtrafo%uses_type12(), "prepared copy deviates from its source")
      if (allocated(error)) return
      deallocate (copy)
      call forward(error, mtrafo, f, fk)
      if (allocated(error)) return
      do i = 1, mgrid%npts_k
         call check(error, complex_is_finite(ref(i, 1)), "reference spectrum is not finite")
         if (allocated(error)) return
         call check(error, fk(i, 1), ref(i, 1), thr=thr, more="source engine deviates after its copy was finalized")
         if (allocated(error)) return
      end do

      ! Factory replacement and deallocation finalize prepared engines
      call mgrid%new_trafo(owned, merr)
      if (.not. allocated(merr)) call owned%prepare(1, merr)
      if (.not. allocated(merr)) call mgrid%new_trafo(owned, merr)
      if (.not. allocated(merr)) call owned%prepare(1, merr)
      call check_moist_error(error, merr)
      if (allocated(error)) return
      call forward(error, owned, f, fk)
      if (allocated(error)) return
      do i = 1, mgrid%npts_k
         call check(error, complex_is_finite(ref(i, 1)), "reference spectrum is not finite")
         if (allocated(error)) return
         call check(error, fk(i, 1), ref(i, 1), thr=thr, more="factory-replaced engine deviates from the first spectrum")
         if (allocated(error)) return
      end do
      deallocate (owned)
      call mtrafo%destroy()
   end subroutine test_molecular_trafo_lifecycle

   !> Build the water probe grid and its single-column transform
   !>
   !> @param[out] mol    carrier water structure, bohr
   !> @param[out] mg     molecular grid with its reciprocal grid configured
   !> @param[out] trafo  transform prepared for one column
   !> @param[out] err    library error from the construction
   subroutine probe_fixture(mol, mg, trafo, err)
      !> Carrier structure
      type(structure_type), intent(out) :: mol
      !> Molecular grid
      type(moist_math_grid_3d_molecular_type), intent(inout), target :: mg
      !> Prepared transform
      type(moist_math_grid_3d_molecular_trafo_type), intent(inout) :: trafo
      !> Library error from the construction
      type(mctc_error), allocatable, intent(out) :: err

      type(moist_math_grid_atomic_recipe_type) :: recipe
      type(moist_math_grid_atomic_recipe_override_type), allocatable :: overrides(:)

      call new(mol, [8, 1, 1], reshape([ &
                                       0.0_wp, 0.0_wp, 0.0_wp, &
                                       1.43_wp, 0.0_wp, 1.11_wp, &
                                       -1.43_wp, 0.0_wp, 1.11_wp], [3, 3]))
      call get_uniform_recipe(recipe, overrides, probe_nrad, probe_nang, err, rmax=probe_rmax)
      if (.not. allocated(err)) call new_molecular_point_grid(mg, err, recipe=recipe, overrides=overrides, &
         & reciprocal=.false.)
      if (.not. allocated(err)) call mg%update(mol, err)
      if (.not. allocated(err)) call molecular_grid_set_kgrid(mg, probe_dr, err)
      if (allocated(err)) return
      call new_molecular_grid_trafo(trafo, mg)
      call trafo%prepare(1, err)
   end subroutine probe_fixture

   !> A prepared trafo refuses to transform after its grid was updated
   !>
   !> - Holds for a pure translation (`ngrid` unchanged), only the generation guard catches it
   subroutine test_trafo_stale_after_update(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(structure_type) :: mol
      type(moist_math_grid_3d_molecular_type), target :: mg
      type(moist_math_grid_3d_molecular_trafo_type) :: trafo
      type(mctc_error), allocatable :: merr
      real(wp), allocatable :: f(:, :)
      complex(wp), allocatable :: fk(:, :)

      call probe_fixture(mol, mg, trafo, merr)
      if (.not. allocated(merr)) then
         mol%xyz(1, 1) = mol%xyz(1, 1) + 0.05_wp
         call mg%update(mol, merr)
      end if
      if (allocated(merr)) then
         call trafo%destroy()
         call mg%destroy()
         return
      end if
      allocate (f(mg%ngrid, 1), fk(mg%npts_k, 1))
      f = 0.0_wp
      call trafo%fft_r2k(f, fk, merr)
      call expect_error(error, merr, "grid geometry changed since prepare")
      call trafo%destroy()
      call mg%destroy()
   end subroutine test_trafo_stale_after_update

   !> A prepared trafo refuses to transform after its grid was constructed again
   !>
   !> - Same recipe, same sequence of update and `molecular_grid_set_kgrid`, one
   !>   atom moved: the point count is unchanged, so only a generation
   !>   counter that survives `new_molecular_point_grid` catches the stale plans
   subroutine test_trafo_stale_after_reconstruction(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(structure_type) :: mol
      type(moist_math_grid_3d_molecular_type), target :: mg
      type(moist_math_grid_3d_molecular_trafo_type) :: trafo
      type(mctc_error), allocatable :: merr
      type(moist_math_grid_atomic_recipe_type) :: recipe
      type(moist_math_grid_atomic_recipe_override_type), allocatable :: overrides(:)
      real(wp), allocatable :: f(:, :)
      complex(wp), allocatable :: fk(:, :)

      call probe_fixture(mol, mg, trafo, merr)
      if (.not. allocated(merr)) then
         mol%xyz(1, 1) = mol%xyz(1, 1) + 0.05_wp
         call get_uniform_recipe(recipe, overrides, probe_nrad, probe_nang, merr, rmax=probe_rmax)
      end if
      if (.not. allocated(merr)) call new_molecular_point_grid(mg, merr, recipe=recipe, overrides=overrides, &
         & reciprocal=.false.)
      if (.not. allocated(merr)) call mg%update(mol, merr)
      if (.not. allocated(merr)) call molecular_grid_set_kgrid(mg, probe_dr, merr)
      if (allocated(merr)) then
         call trafo%destroy()
         call mg%destroy()
         return
      end if
      allocate (f(mg%ngrid, 1), fk(mg%npts_k, 1))
      f = 0.0_wp
      call trafo%fft_r2k(f, fk, merr)
      call expect_error(error, merr, "grid geometry changed since prepare")
      call trafo%destroy()
      call mg%destroy()
   end subroutine test_trafo_stale_after_reconstruction

   !> Every mis-sized add_weights call fails by name
   !>
   !> - Uninitialized accumulator, position weights of the wrong shape,
   !>   integration and Gaussian-width weights of the wrong length
   subroutine test_volume_adjoint_sizes(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(volume_adjoint_type) :: acc
      type(mctc_error), allocatable :: merr
      real(wp) :: w_xyz(2, 4)

      call acc%add_weights(merr, w_w=[1.0_wp])
      call require_message(error, merr, "not initialized", "uninitialized accumulator")
      if (allocated(error)) return
      call acc%init(4)
      w_xyz = 1.0_wp
      call acc%add_weights(merr, w_xyz=w_xyz)
      call require_message(error, merr, "xyz weight shape mismatch", "position weights")
      if (allocated(error)) return
      call acc%add_weights(merr, w_w=[1.0_wp, 2.0_wp, 3.0_wp])
      call require_message(error, merr, "integration weight size mismatch", "integration weights")
      if (allocated(error)) return
      call acc%add_weights(merr, w_xi=[1.0_wp, 2.0_wp, 3.0_wp, 4.0_wp, 5.0_wp])
      call require_message(error, merr, "xi weight size mismatch", "Gaussian-width weights")
   end subroutine test_volume_adjoint_sizes

   !> Point grid rejects a nonzero Gaussian-width adjoint and leaves the gradient unchanged
   subroutine test_point_grid_width_adjoint_fails(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(moist_math_grid_3d_cartesian_type) :: grid
      type(volume_adjoint_type) :: acc
      type(structure_type) :: mol
      type(mctc_error), allocatable :: merr
      real(wp) :: gradient(3, 3)

      gradient = 0.25_wp
      call make_molecule(mol)
      call new_cartesian_point_grid(grid, merr, nx=4, ny=4, nz=4, margin=0.0_wp)
      if (.not. allocated(merr)) call grid%update(mol, merr)
      if (allocated(merr)) return
      call acc%init(grid%ngrid)
      acc%w_xi = 1.0_wp
      call grid%get_volume_gradient(acc, gradient, merr)
      if (any(gradient /= 0.25_wp)) return
      call expect_error(error, merr, "point grids have no Gaussian-width channel")
   end subroutine test_point_grid_width_adjoint_fails

   !> Molecular Gaussian-grid Hessian rejects width adjoints and leaves the Hessian unchanged
   subroutine test_hessian_width_adjoint_fails(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(moist_math_grid_3d_molecular_type) :: grid
      type(moist_math_grid_atomic_recipe_type) :: recipe
      type(volume_adjoint_type) :: acc
      type(structure_type) :: mol
      type(mctc_error), allocatable :: merr
      real(wp) :: hessian(3, 3, 3, 3)

      call make_molecule(mol)
      hessian = 0.25_wp
      call get_qc_handymod_recipe(recipe, merr, nrad=3, degree=5)
      if (.not. allocated(merr)) call new_molecular_gaussian_grid(grid, merr, recipe=recipe, reciprocal=.false.)
      if (.not. allocated(merr)) call grid%update(mol, merr)
      if (allocated(merr)) return
      call acc%init(grid%ngrid)
      acc%w_xi = 0.3_wp
      call grid%get_volume_hessian(acc, hessian, merr)
      if (any(hessian /= 0.25_wp)) return
      call expect_error(error, merr, "Gaussian-width curvature is unsupported")
   end subroutine test_hessian_width_adjoint_fails

   !> Becke stiffness one rejects a seeded sample on a nucleus and leaves the product unchanged
   subroutine test_hessian_becke1_nucleus_fails(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(mctc_error), allocatable :: merr
      type(becke_partition_type) :: becke
      real(wp) :: xyz(3, 2), points(3, 1), direction(3, 2), hvp(3, 2), point_hvp(3, 1)

      xyz = 0.0_wp
      xyz(1, :) = [-1.0_wp, 1.0_wp]
      direction = reshape([0.1_wp, 0.3_wp, -0.2_wp, -0.4_wp, 0.2_wp, 0.5_wp], [3, 2])
      points(:, 1) = xyz(:, 1)
      hvp = 0.25_wp
      becke%k = 1
      call becke%owner_hessian_vector(1, points, xyz, [1, 1], [1.0_wp], direction, &
         & [0.0_wp, 0.0_wp, 0.0_wp], hvp, point_hvp, merr)
      if (any(hvp /= 0.25_wp)) return
      call expect_error(error, merr, "k=1 is not twice differentiable at a nucleus")
   end subroutine test_hessian_becke1_nucleus_fails

   !> Coincident elliptic centers are rejected and leave the product unchanged
   subroutine test_hessian_coincident_fails(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(mctc_error), allocatable :: merr
      type(becke_partition_type) :: becke
      real(wp) :: xyz(3, 2), points(3, 1), direction(3, 2), hvp(3, 2), point_hvp(3, 1)

      xyz = 0.0_wp
      xyz(1, :) = [-1.0_wp, 1.0_wp]
      direction = reshape([0.1_wp, 0.3_wp, -0.2_wp, -0.4_wp, 0.2_wp, 0.5_wp], [3, 2])
      points(:, 1) = xyz(:, 1)
      xyz(:, 2) = xyz(:, 1)
      hvp = 0.25_wp
      becke%k = 3
      call becke%owner_hessian_vector(1, points, xyz, [1, 1], [1.0_wp], direction, &
         & [0.0_wp, 0.0_wp, 0.0_wp], hvp, point_hvp, merr)
      if (any(hvp /= 0.25_wp)) return
      call expect_error(error, merr, "coincident elliptic centers are not differentiable")
   end subroutine test_hessian_coincident_fails

   !> Power radii of the wrong length are rejected before any geometry is read
   !>
   !> - The owner gradient must name the mismatch and leave its accumulator
   !>   unchanged before the whole-grid weights are tried
   !> - The test raises only when the whole-grid weights also name the mismatch
   subroutine test_power_radii_size_fails(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(pvoronoi_partition_type) :: pvoronoi
      type(mctc_error), allocatable :: merr
      real(wp) :: xyz(3, 2), points(3, 1), gradient(3, 2), point_gradient(3, 1), w(1), partition_w(1)
      character(len=*), parameter :: expected = "radii needs one entry per atom"

      xyz = 0.0_wp
      xyz(1, 2) = 2.0_wp
      points(:, 1) = [0.5_wp, 0.1_wp, 0.0_wp]
      pvoronoi%radii = [1.0_wp, 1.2_wp, 0.9_wp]
      gradient = 0.25_wp
      call pvoronoi%owner_gradient(1, points, xyz, [1, 1], [1.0_wp], gradient, point_gradient, merr)
      if (.not. allocated(merr)) return
      if (index(merr%message, expected) == 0 .or. any(gradient /= 0.25_wp)) return
      deallocate (merr)
      call pvoronoi%grid_weights(points, xyz, [1, 1], [1], [1.0_wp], w, partition_w, merr)
      call expect_error(error, merr, expected)
   end subroutine test_power_radii_size_fails
end module test_math_grid_3d_guards
