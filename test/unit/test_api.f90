!> Defensive-preamble tests for the C API entry points
module test_api
   use iso_c_binding, only: c_ptr, c_loc, c_null_ptr, c_int, c_double, c_bool, &
                            c_funptr, c_funloc, c_associated, c_f_pointer
   use mctc_env_error, only: moist_error_type => error_type
   use moist_api, only: vp_cavity, vp_error
   use testdrive, only: new_unittest, unittest_type, error_type, check, test_failed
   implicit none
   private

   public :: collect_api

   !> C bindings of the entry points under test
   interface
      subroutine moist_update_drop_cavity(verror, vcav, vmol, nleb) &
            & bind(C, name="moist_update_drop_cavity")
         import :: c_ptr
         type(c_ptr), value :: verror
         type(c_ptr), value :: vcav
         type(c_ptr), value :: vmol
         type(c_ptr), value :: nleb
      end subroutine moist_update_drop_cavity

      subroutine moist_get_drop_results(verror, vcav, ngrid_cap, nsph_cap, &
            & area, volume, ngrid, nmax, nsph, &
            & xyz, normal0, wleb, a, r_iI0, f, rho, owner, converged, radii, asph) &
            & bind(C, name="moist_get_drop_results")
         import :: c_ptr, c_int, c_double, c_bool
         type(c_ptr), value :: verror
         type(c_ptr), value :: vcav
         integer(c_int), value :: ngrid_cap
         integer(c_int), value :: nsph_cap
         real(c_double), intent(out) :: area
         real(c_double), intent(out) :: volume
         integer(c_int), intent(out) :: ngrid
         integer(c_int), intent(out) :: nmax
         integer(c_int), intent(out) :: nsph
         real(c_double), intent(out) :: xyz(3, ngrid_cap)
         real(c_double), intent(out) :: normal0(3, ngrid_cap)
         real(c_double), intent(out) :: wleb(ngrid_cap)
         real(c_double), intent(out) :: a(ngrid_cap)
         real(c_double), intent(out) :: r_iI0(ngrid_cap)
         real(c_double), intent(out) :: f(ngrid_cap)
         real(c_double), intent(out) :: rho(ngrid_cap)
         integer(c_int), intent(out) :: owner(ngrid_cap)
         logical(c_bool), intent(out) :: converged(ngrid_cap)
         real(c_double), intent(out) :: radii(nsph_cap)
         real(c_double), intent(out) :: asph(nsph_cap)
      end subroutine moist_get_drop_results

      subroutine moist_get_cavity_sizes(verror, vcav, ngrid, nsph) &
            & bind(C, name="moist_get_cavity_sizes")
         import :: c_ptr, c_int
         type(c_ptr), value :: verror
         type(c_ptr), value :: vcav
         integer(c_int), intent(out) :: ngrid
         integer(c_int), intent(out) :: nsph
      end subroutine moist_get_cavity_sizes

      function moist_new_structure(verror, natoms, numbers, positions, &
            & lattice, periodic) result(vmol) &
            & bind(C, name="moist_new_structure")
         import :: c_ptr, c_int, c_double, c_bool
         type(c_ptr), value :: verror
         integer(c_int), value :: natoms
         integer(c_int), intent(in) :: numbers(natoms)
         real(c_double), intent(in) :: positions(3, natoms)
         real(c_double), intent(in), optional :: lattice(3, 3)
         logical(c_bool), intent(in), optional :: periodic(3)
         type(c_ptr) :: vmol
      end function moist_new_structure

      subroutine moist_delete_structure(vmol) &
            & bind(C, name="moist_delete_structure")
         import :: c_ptr
         type(c_ptr), intent(inout) :: vmol
      end subroutine moist_delete_structure

      function moist_new_drop_cavity(verror, nleb, debug, verbose, blendk, &
            & blend1b, blend2b, blend3b, do_fine, tolerance, proj_maxiter, &
            & proj_level, branch_weight_s, rho_grid_h, wleb_prune_level) result(vcav) &
            & bind(C, name="moist_new_drop_cavity")
         import :: c_ptr
         type(c_ptr), value :: verror
         type(c_ptr), value :: nleb
         type(c_ptr), value :: debug
         type(c_ptr), value :: verbose
         type(c_ptr), value :: blendk
         type(c_ptr), value :: blend1b
         type(c_ptr), value :: blend2b
         type(c_ptr), value :: blend3b
         type(c_ptr), value :: do_fine
         type(c_ptr), value :: tolerance
         type(c_ptr), value :: proj_maxiter
         type(c_ptr), value :: proj_level
         type(c_ptr), value :: branch_weight_s
         type(c_ptr), value :: rho_grid_h
         type(c_ptr), value :: wleb_prune_level
         type(c_ptr) :: vcav
      end function moist_new_drop_cavity

      subroutine moist_update_cavity(verror, vcav, vmol) &
            & bind(C, name="moist_update_cavity")
         import :: c_ptr
         type(c_ptr), value :: verror
         type(c_ptr), value :: vcav
         type(c_ptr), value :: vmol
      end subroutine moist_update_cavity

      subroutine moist_delete_cavity(vcav) &
            & bind(C, name="moist_delete_cavity")
         import :: c_ptr
         type(c_ptr), intent(inout) :: vcav
      end subroutine moist_delete_cavity

      subroutine moist_get_cavity_results(verror, vcav, ngrid_cap, nsph_cap, &
            & area, volume, ngrid, nsph, xyz, a, owner, converged, radii, asph) &
            & bind(C, name="moist_get_cavity_results")
         import :: c_ptr, c_int, c_double, c_bool
         type(c_ptr), value :: verror
         type(c_ptr), value :: vcav
         integer(c_int), value :: ngrid_cap
         integer(c_int), value :: nsph_cap
         real(c_double), intent(out) :: area
         real(c_double), intent(out) :: volume
         integer(c_int), intent(out) :: ngrid
         integer(c_int), intent(out) :: nsph
         real(c_double), intent(out) :: xyz(3, ngrid_cap)
         real(c_double), intent(out) :: a(ngrid_cap)
         integer(c_int), intent(out) :: owner(ngrid_cap)
         logical(c_bool), intent(out) :: converged(ngrid_cap)
         real(c_double), intent(out) :: radii(nsph_cap)
         real(c_double), intent(out) :: asph(nsph_cap)
      end subroutine moist_get_cavity_results

      subroutine moist_get_drop_specific(verror, vcav, ngrid_cap, nmax, &
            & normal0, wleb, r_iI0, f, rho) &
            & bind(C, name="moist_get_drop_specific")
         import :: c_ptr, c_int, c_double
         type(c_ptr), value :: verror
         type(c_ptr), value :: vcav
         integer(c_int), value :: ngrid_cap
         integer(c_int), intent(out) :: nmax
         real(c_double), intent(out) :: normal0(3, ngrid_cap)
         real(c_double), intent(out) :: wleb(ngrid_cap)
         real(c_double), intent(out) :: r_iI0(ngrid_cap)
         real(c_double), intent(out) :: f(ngrid_cap)
         real(c_double), intent(out) :: rho(ngrid_cap)
      end subroutine moist_get_drop_specific

      subroutine moist_get_drop_numbering(verror, vcav, ngrid_cap, numbering) &
            & bind(C, name="moist_get_drop_numbering")
         import :: c_ptr, c_int
         type(c_ptr), value :: verror
         type(c_ptr), value :: vcav
         integer(c_int), value :: ngrid_cap
         integer(c_int), intent(out) :: numbering(ngrid_cap)
      end subroutine moist_get_drop_numbering

      subroutine moist_assemble_amat(verror, vcav, ngrid_cap, amat0, xi) &
            & bind(C, name="moist_assemble_amat")
         import :: c_ptr, c_int, c_double
         type(c_ptr), value :: verror
         type(c_ptr), value :: vcav
         integer(c_int), value :: ngrid_cap
         real(c_double), intent(out) :: amat0(ngrid_cap, ngrid_cap)
         real(c_double), intent(out) :: xi(ngrid_cap)
      end subroutine moist_assemble_amat

      function moist_new_drop_cavity_isodensity_callback(verror, callback, context, &
            & rho_iso, scale, nleb, debug, verbose, do_fine, wleb_prune_level, tolerance) &
            & result(vcav) &
            & bind(C, name="moist_new_drop_cavity_isodensity_callback")
         import :: c_ptr, c_funptr, c_double
         type(c_ptr), value :: verror
         type(c_funptr), value :: callback
         type(c_ptr), value :: context
         real(c_double), value :: rho_iso
         type(c_ptr), value :: scale
         type(c_ptr), value :: nleb
         type(c_ptr), value :: debug
         type(c_ptr), value :: verbose
         type(c_ptr), value :: do_fine
         type(c_ptr), value :: wleb_prune_level
         type(c_ptr), value :: tolerance
         type(c_ptr) :: vcav
      end function moist_new_drop_cavity_isodensity_callback
   end interface

   !* ================================================================================= *!
   !*                    Isodensity callback state (failure-channel tests)              *!
   !* ================================================================================= *!

   !> Water geometry the test callback's model density is centered on (Bohr)
   real(c_double), parameter :: cb_centers(3, 3) = reshape( &
      [0.0_c_double, 0.0_c_double, 0.1173_c_double, &
       0.0_c_double, 1.4309_c_double, -0.9370_c_double, &
       0.0_c_double, -1.4309_c_double, -0.9370_c_double], [3, 3])
   !> Exponent of the per-atom s-Gaussian density, 2*alpha for alpha = 0.3
   real(c_double), parameter :: cb_a = 0.6_c_double
   !> Prefactor of the per-atom s-Gaussian density, 2*((2*alpha/pi)**0.75)**2
   real(c_double), parameter :: cb_c = 0.166929_c_double
   !> Isovalue defining the surface (Bohr^-3); puts it ~2.9 Bohr off each atom
   real(c_double), parameter :: cb_rho_iso = 1.0e-3_c_double

   !> Per-test callback state.
   !>
   !> Lives in the callback's own `context` rather than in module variables:
   !> test-drive runs the tests of a suite in an OpenMP parallel loop, so two
   !> tests using the same callback would otherwise share one counter.
   type :: iso_cb_ctx
      !> Evaluations answered before the callback starts reporting failure
      integer :: fail_after = huge(1)
      !> Status the callback reports once it starts failing
      integer(c_int) :: fail_status = 0_c_int
      !> Evaluation counter, bumped atomically: moist calls the callback from
      !> several OpenMP threads within one build
      integer :: calls = 0
   end type iso_cb_ctx

contains

   !> Collect all C API guard tests.
   subroutine collect_api(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("update_drop_cavity_null_handle", test_update_drop_null_handle), &
                  new_unittest("update_drop_cavity_null_inner", test_update_drop_null_inner), &
                  new_unittest("get_drop_results_null_handle", test_drop_results_null_handle), &
                  new_unittest("get_drop_results_null_inner", test_drop_results_null_inner), &
                  new_unittest("get_cavity_sizes_null_handle", test_cavity_sizes_null_handle), &
                  new_unittest("get_cavity_sizes_null_inner", test_cavity_sizes_null_inner), &
                  new_unittest("array_capacity_too_small", test_capacity_too_small), &
                  new_unittest("array_capacity_oversized", test_capacity_oversized), &
                  new_unittest("isodensity_callback_fails_first_call", test_iso_callback_fails_first), &
                  new_unittest("isodensity_callback_fails_mid_loop", test_iso_callback_fails_mid_loop) &
                  ]

   end subroutine collect_api

   !> Assert that the API error handle carries a message containing `needle`.
   subroutine check_api_error(error, err, needle)
      !> Test-drive error
      type(error_type), allocatable, intent(out) :: error
      !> API error handle produced by the entry point
      type(vp_error), pointer, intent(in) :: err
      !> Expected message fragment
      character(len=*), intent(in) :: needle

      if (.not. allocated(err%ptr)) then
         call test_failed(error, "API entry point did not report an error")
         return
      end if
      if (index(err%ptr%message, needle) <= 0) then
         call test_failed(error, "Unexpected API error message: "//err%ptr%message)
         return
      end if

   end subroutine check_api_error

   !> A missing DROP cavity handle is reported, not dereferenced.
   subroutine test_update_drop_null_handle(error)
      type(error_type), allocatable, intent(out) :: error
      type(vp_error), pointer :: err

      allocate (err)
      call moist_update_drop_cavity(c_loc(err), c_null_ptr, c_null_ptr, c_null_ptr)
      call check_api_error(error, err, "DROP cavity handle is missing")
      deallocate (err)

   end subroutine test_update_drop_null_handle

   !> A DROP cavity handle with a null inner pointer is reported, not dereferenced.
   subroutine test_update_drop_null_inner(error)
      type(error_type), allocatable, intent(out) :: error
      type(vp_error), pointer :: err
      type(vp_cavity), pointer :: cav

      allocate (err)
      allocate (cav)
      call moist_update_drop_cavity(c_loc(err), c_loc(cav), c_null_ptr, c_null_ptr)
      call check_api_error(error, err, "Cavity is not initialized")
      deallocate (cav)
      deallocate (err)

   end subroutine test_update_drop_null_inner

   !> A missing DROP cavity handle is reported, not dereferenced.
   subroutine test_drop_results_null_handle(error)
      type(error_type), allocatable, intent(out) :: error
      type(vp_error), pointer :: err

      allocate (err)
      call call_get_drop_results(c_loc(err), c_null_ptr)
      call check_api_error(error, err, "DROP cavity handle is missing")
      deallocate (err)

   end subroutine test_drop_results_null_handle

   !> A DROP cavity handle with a null inner pointer is reported, not dereferenced.
   subroutine test_drop_results_null_inner(error)
      type(error_type), allocatable, intent(out) :: error
      type(vp_error), pointer :: err
      type(vp_cavity), pointer :: cav

      allocate (err)
      allocate (cav)
      call call_get_drop_results(c_loc(err), c_loc(cav))
      call check_api_error(error, err, "Cavity is not initialized")
      deallocate (cav)
      deallocate (err)

   end subroutine test_drop_results_null_inner

   !> A missing cavity handle is reported, not dereferenced.
   subroutine test_cavity_sizes_null_handle(error)
      type(error_type), allocatable, intent(out) :: error
      type(vp_error), pointer :: err
      integer(c_int) :: ngrid, nsph

      allocate (err)
      call moist_get_cavity_sizes(c_loc(err), c_null_ptr, ngrid, nsph)
      call check_api_error(error, err, "Cavity handle is missing")
      deallocate (err)

   end subroutine test_cavity_sizes_null_handle

   !> A cavity handle with a null inner pointer is reported, not dereferenced.
   subroutine test_cavity_sizes_null_inner(error)
      type(error_type), allocatable, intent(out) :: error
      type(vp_error), pointer :: err
      type(vp_cavity), pointer :: cav
      integer(c_int) :: ngrid, nsph

      allocate (err)
      allocate (cav)
      call moist_get_cavity_sizes(c_loc(err), c_loc(cav), ngrid, nsph)
      call check_api_error(error, err, "Cavity is not initialized")
      deallocate (cav)
      deallocate (err)

   end subroutine test_cavity_sizes_null_inner

   !> Provide the caller-owned result buffers moist_get_drop_results writes into.
   !> The guarded entry point returns before touching them; they only have to
   !> exist so the call is well formed.
   subroutine call_get_drop_results(verror, vcav)
      !> Opaque error handle
      type(c_ptr), intent(in) :: verror
      !> Opaque cavity handle under test
      type(c_ptr), intent(in) :: vcav
      !> Scalar result buffers
      real(c_double) :: area, volume
      integer(c_int) :: ngrid, nmax, nsph
      !> Per-point result buffers
      real(c_double) :: xyz(3, 1), normal0(3, 1)
      real(c_double) :: wleb(1), a(1), r_iI0(1), f(1), rho(1)
      integer(c_int) :: owner(1)
      logical(c_bool) :: converged(1)
      !> Per-sphere result buffers
      real(c_double) :: radii(1), asph(1)

      call moist_get_drop_results(verror, vcav, 1_c_int, 1_c_int, &
                                  area, volume, ngrid, nmax, nsph, &
                                  xyz, normal0, wleb, a, r_iI0, f, rho, owner, &
                                  converged, radii, asph)

   end subroutine call_get_drop_results

   !> Build a DROP cavity for water through the public C entry points, exactly
   !> as a host would.  The caller owns the returned handles.
   subroutine build_water_cavity(error, err, verror, vmol, vcav, ngrid, nsph)
      !> Test-drive error
      type(error_type), allocatable, intent(out) :: error
      !> API error handle (allocated here, deallocated by the caller)
      type(vp_error), pointer, intent(out) :: err
      !> Opaque view of `err`
      type(c_ptr), intent(out) :: verror
      !> Structure and cavity handles
      type(c_ptr), intent(out) :: vmol, vcav
      !> Sizes of the built cavity
      integer(c_int), intent(out) :: ngrid, nsph
      !> Water, coordinates in Bohr
      integer(c_int), parameter :: numbers(3) = [8, 1, 1]
      real(c_double), parameter :: positions(3, 3) = reshape( &
         [0.0_c_double, 0.0_c_double, 0.1173_c_double, &
          0.0_c_double, 1.4309_c_double, -0.9370_c_double, &
          0.0_c_double, -1.4309_c_double, -0.9370_c_double], [3, 3])
      !> Lebedev order, kept small so the build stays cheap
      integer(c_int), target :: nleb

      ngrid = 0
      nsph = 0
      nleb = 26_c_int

      allocate (err)
      verror = c_loc(err)

      vmol = moist_new_structure(verror, 3_c_int, numbers, positions)
      vcav = moist_new_drop_cavity(verror, c_loc(nleb), c_null_ptr, c_null_ptr, &
                                   c_null_ptr, c_null_ptr, c_null_ptr, c_null_ptr, &
                                   c_null_ptr, c_null_ptr, c_null_ptr, c_null_ptr, &
                                   c_null_ptr, c_null_ptr, c_null_ptr)
      call moist_update_cavity(verror, vcav, vmol)
      call moist_get_cavity_sizes(verror, vcav, ngrid, nsph)

      if (allocated(err%ptr)) then
         call test_failed(error, "Cavity setup failed: "//err%ptr%message)
         return
      end if
      if (ngrid <= 1 .or. nsph /= 3) then
         call test_failed(error, "Unexpected cavity sizes for water")
         return
      end if

   end subroutine build_water_cavity

   !> Release the handles built by build_water_cavity.
   subroutine drop_water_cavity(err, vmol, vcav)
      type(vp_error), pointer, intent(inout) :: err
      type(c_ptr), intent(inout) :: vmol, vcav

      call moist_delete_cavity(vcav)
      call moist_delete_structure(vmol)
      deallocate (err)

   end subroutine drop_water_cavity

   !> Assert that an entry point rejected an undersized capacity, then clear the
   !> error so the handle can be reused.
   subroutine expect_capacity_error(error, err, routine)
      type(error_type), allocatable, intent(out) :: error
      type(vp_error), pointer, intent(inout) :: err
      character(len=*), intent(in) :: routine

      if (.not. allocated(err%ptr)) then
         call test_failed(error, routine//" accepted an undersized array capacity")
         return
      end if
      if (index(err%ptr%message, "Array capacity too small") <= 0) then
         call test_failed(error, "Unexpected API error message: "//err%ptr%message)
         return
      end if
      deallocate (err%ptr)

   end subroutine expect_capacity_error

   !> A capacity below the cavity's own ngrid/nsph is refused before a single
   !> element is written.  Without the check these calls write ngrid values into
   !> buffers holding ngrid-1, which is silent heap corruption; the sentinels
   !> below detect any write at all.
   subroutine test_capacity_too_small(error)
      type(error_type), allocatable, intent(out) :: error
      type(vp_error), pointer :: err
      type(c_ptr) :: verror, vmol, vcav
      integer(c_int) :: ngrid, nsph, cap
      integer(c_int) :: out_ngrid, out_nsph, nmax
      real(c_double) :: area, volume
      real(c_double), allocatable :: xyz(:, :), a(:), radii(:), asph(:)
      real(c_double), allocatable :: normal0(:, :), wleb(:), r_iI0(:), f(:), rho(:)
      real(c_double), allocatable :: amat0(:, :), xi(:)
      integer(c_int), allocatable :: owner(:), numbering(:)
      logical(c_bool), allocatable :: converged(:)
      real(c_double), parameter :: sentinel = -12345.0_c_double

      call build_water_cavity(error, err, verror, vmol, vcav, ngrid, nsph)
      if (allocated(error)) return

      cap = ngrid - 1_c_int
      allocate (xyz(3, cap), a(cap), owner(cap), converged(cap))
      allocate (normal0(3, cap), wleb(cap), r_iI0(cap), f(cap), rho(cap))
      allocate (numbering(cap), amat0(cap, cap), xi(cap))
      allocate (radii(nsph), asph(nsph))
      xyz = sentinel
      a = sentinel
      owner = -1_c_int
      converged = .false._c_bool
      normal0 = sentinel
      wleb = sentinel
      r_iI0 = sentinel
      f = sentinel
      rho = sentinel
      numbering = -1_c_int
      amat0 = sentinel
      xi = sentinel
      radii = sentinel
      asph = sentinel

      call moist_get_cavity_results(verror, vcav, cap, nsph, area, volume, &
                                    out_ngrid, out_nsph, xyz, a, owner, converged, &
                                    radii, asph)
      call expect_capacity_error(error, err, "get_cavity_results")
      if (.not. allocated(error)) then
         if (any(a /= sentinel) .or. any(xyz /= sentinel) .or. any(owner /= -1_c_int)) &
            call test_failed(error, "get_cavity_results wrote into a rejected buffer")
      end if

      if (.not. allocated(error)) then
         call moist_get_drop_specific(verror, vcav, cap, nmax, normal0, wleb, &
                                      r_iI0, f, rho)
         call expect_capacity_error(error, err, "get_drop_specific")
      end if
      if (.not. allocated(error)) then
         if (any(wleb /= sentinel) .or. any(normal0 /= sentinel)) &
            call test_failed(error, "get_drop_specific wrote into a rejected buffer")
      end if

      if (.not. allocated(error)) then
         call moist_get_drop_numbering(verror, vcav, cap, numbering)
         call expect_capacity_error(error, err, "get_drop_numbering")
      end if
      if (.not. allocated(error)) then
         if (any(numbering /= -1_c_int)) &
            call test_failed(error, "get_drop_numbering wrote into a rejected buffer")
      end if

      if (.not. allocated(error)) then
         call moist_assemble_amat(verror, vcav, cap, amat0, xi)
         call expect_capacity_error(error, err, "assemble_amat")
      end if
      if (.not. allocated(error)) then
         if (any(amat0 /= sentinel) .or. any(xi /= sentinel)) &
            call test_failed(error, "assemble_amat wrote into a rejected buffer")
      end if

      ! A short per-sphere capacity has to be caught the same way
      if (.not. allocated(error)) then
         deallocate (radii, asph)
         allocate (radii(nsph - 1), asph(nsph - 1))
         radii = sentinel
         asph = sentinel
         call moist_get_cavity_results(verror, vcav, ngrid, nsph - 1_c_int, area, &
                                       volume, out_ngrid, out_nsph, xyz, a, owner, &
                                       converged, radii, asph)
         call expect_capacity_error(error, err, "get_cavity_results (nsph)")
      end if

      call drop_water_cavity(err, vmol, vcav)

   end subroutine test_capacity_too_small

   !> A capacity above the cavity's own sizes is accepted: a host may allocate
   !> one max-size buffer and reuse it.  Data lands in the leading entries, with
   !> the capacity setting the stride of the vector array.
   subroutine test_capacity_oversized(error)
      type(error_type), allocatable, intent(out) :: error
      type(vp_error), pointer :: err
      type(c_ptr) :: verror, vmol, vcav
      integer(c_int) :: ngrid, nsph, ngrid_cap, nsph_cap
      integer(c_int) :: out_ngrid, out_nsph
      real(c_double) :: area, volume
      real(c_double), allocatable :: xyz(:, :), a(:), radii(:), asph(:)
      integer(c_int), allocatable :: owner(:)
      logical(c_bool), allocatable :: converged(:)
      real(c_double), parameter :: sentinel = -12345.0_c_double

      call build_water_cavity(error, err, verror, vmol, vcav, ngrid, nsph)
      if (allocated(error)) return

      ngrid_cap = ngrid + 17_c_int
      nsph_cap = nsph + 5_c_int
      allocate (xyz(3, ngrid_cap), a(ngrid_cap), owner(ngrid_cap), converged(ngrid_cap))
      allocate (radii(nsph_cap), asph(nsph_cap))
      xyz = sentinel
      a = sentinel
      radii = sentinel
      asph = sentinel

      call moist_get_cavity_results(verror, vcav, ngrid_cap, nsph_cap, area, volume, &
                                    out_ngrid, out_nsph, xyz, a, owner, converged, &
                                    radii, asph)

      if (allocated(err%ptr)) then
         call test_failed(error, "Oversized capacity was rejected: "//err%ptr%message)
      else
         call check(error, out_ngrid, ngrid)
         if (.not. allocated(error)) call check(error, out_nsph, nsph)
         if (.not. allocated(error)) call check(error, area > 0.0_c_double)
         ! The leading entries carry the data ...
         if (.not. allocated(error)) call check(error, a(1) > 0.0_c_double)
         if (.not. allocated(error)) call check(error, radii(nsph) > 0.0_c_double)
         ! ... and the tail beyond the cavity sizes is left untouched.
         if (.not. allocated(error)) call check(error, a(ngrid + 1) == sentinel)
         if (.not. allocated(error)) call check(error, radii(nsph + 1) == sentinel)
         if (.not. allocated(error)) call check(error, xyz(1, ngrid + 1) == sentinel)
      end if

      call drop_water_cavity(err, vmol, vcav)

   end subroutine test_capacity_oversized

   !* ================================================================================= *!
   !*                        Isodensity callback failure channel                        *!
   !* ================================================================================= *!

   !> Isodensity LSF callback with a switchable failure channel.
   !>
   !> Answers the first `ctx%fail_after` evaluations with a smooth three-Gaussian
   !> model density (so the projection behaves like a real build), then reports
   !> `ctx%fail_status` for every evaluation after that.  The counter is bumped
   !> atomically because moist calls this from inside OpenMP loops.
   !>
   !> @param[in]  context Callback context, an `iso_cb_ctx`
   !> @param[in]  point   Evaluation point in Bohr
   !> @param[out] rho     Electron density
   !> @param[out] drho    Density spatial gradient
   !> @param[out] d2rho   Density spatial Hessian (Fortran (3,3), or NULL)
   !> @param[out] d3rho   Density third spatial derivative (Fortran (3,3,3), or NULL)
   !> @returns            0 while healthy, `ctx%fail_status` once failing
   function iso_switchable_callback(context, point, rho, drho, d2rho, d3rho) &
      result(status) bind(C)
      type(c_ptr), value :: context
      real(c_double), intent(in) :: point(3)
      real(c_double), intent(out) :: rho
      real(c_double), intent(out) :: drho(3)
      type(c_ptr), value :: d2rho
      type(c_ptr), value :: d3rho
      integer(c_int) :: status

      type(iso_cb_ctx), pointer :: ctx
      real(c_double), pointer :: hptr(:, :), tptr(:, :, :)
      real(c_double) :: d2rho_l(3, 3), d3rho_l(3, 3, 3)
      real(c_double) :: d(3), g
      integer :: iatom, i, j, k, ncalls
      logical :: want_hess, want_third

      if (.not. c_associated(context)) then
         status = -1_c_int
         return
      end if
      call c_f_pointer(context, ctx)

      !$omp atomic capture
      ctx%calls = ctx%calls + 1
      ncalls = ctx%calls
      !$omp end atomic

      if (ncalls > ctx%fail_after) then
         status = ctx%fail_status
         return
      end if

      status = 0_c_int
      want_hess = c_associated(d2rho)
      want_third = c_associated(d3rho)

      rho = 0.0_c_double
      drho = 0.0_c_double
      d2rho_l = 0.0_c_double
      d3rho_l = 0.0_c_double

      do iatom = 1, 3
         d = point - cb_centers(:, iatom)
         g = cb_c*exp(-cb_a*dot_product(d, d))
         rho = rho + g
         drho = drho - 2.0_c_double*cb_a*d*g
         if (want_hess) then
            do j = 1, 3
               do i = 1, 3
                  d2rho_l(i, j) = d2rho_l(i, j) + (4.0_c_double*cb_a*cb_a*d(i)*d(j) &
                                                   - 2.0_c_double*cb_a*delta(i, j))*g
               end do
            end do
         end if
         if (want_third) then
            do k = 1, 3
               do j = 1, 3
                  do i = 1, 3
                     d3rho_l(i, j, k) = d3rho_l(i, j, k) &
                                        + (-8.0_c_double*cb_a**3*d(i)*d(j)*d(k) &
                                           + 4.0_c_double*cb_a*cb_a*(d(i)*delta(j, k) &
                                                                     + d(j)*delta(i, k) &
                                                                     + d(k)*delta(i, j)))*g
                  end do
               end do
            end do
         end if
      end do

      ! The bare density: moist applies the isovalue and the DROP sign convention
      if (want_hess) then
         call c_f_pointer(d2rho, hptr, [3, 3])
         hptr = d2rho_l
      end if
      if (want_third) then
         call c_f_pointer(d3rho, tptr, [3, 3, 3])
         tptr = d3rho_l
      end if

   end function iso_switchable_callback

   !> Kronecker delta as a c_double, keeping the derivative expressions readable
   pure function delta(i, j) result(d)
      integer, intent(in) :: i, j
      real(c_double) :: d

      d = 0.0_c_double
      if (i == j) d = 1.0_c_double
   end function delta

   !> Build a water cavity backed by `iso_switchable_callback` and update it.
   !>
   !> Returns the handles so the caller can inspect the error and delete them.
   !>
   !> @param[out]   err        API error handle (caller deallocates)
   !> @param[out]   verror     Opaque view of `err`
   !> @param[out]   vmol,vcav  Structure and cavity handles (caller deletes)
   !> @param[inout] ctx        Callback state; must outlive the cavity handle
   subroutine run_iso_callback_update(err, verror, vmol, vcav, ctx)
      type(vp_error), pointer, intent(out) :: err
      type(c_ptr), intent(out) :: verror
      type(c_ptr), intent(out) :: vmol, vcav
      type(iso_cb_ctx), intent(inout), target :: ctx

      integer(c_int), parameter :: numbers(3) = [8, 1, 1]
      !> Lebedev order, kept small so the build stays cheap
      integer(c_int), target :: nleb

      nleb = 26_c_int
      ctx%calls = 0

      allocate (err)
      verror = c_loc(err)

      vmol = moist_new_structure(verror, 3_c_int, numbers, cb_centers)
      vcav = moist_new_drop_cavity_isodensity_callback(verror, &
                                                       c_funloc(iso_switchable_callback), &
                                                       c_loc(ctx), cb_rho_iso, c_null_ptr, &
                                                       c_loc(nleb), &
                                                       c_null_ptr, c_null_ptr, c_null_ptr, &
                                                       c_null_ptr, c_null_ptr)
      call moist_update_cavity(verror, vcav, vmol)

   end subroutine run_iso_callback_update

   !> Assert that the update failed with the callback's own status in the message
   !>
   !> @param[out]   error  Test-drive error
   !> @param[in]    err    API error handle produced by the update
   !> @param[in]    needle Status fragment the message has to name
   subroutine expect_iso_callback_error(error, err, needle)
      type(error_type), allocatable, intent(out) :: error
      type(vp_error), pointer, intent(in) :: err
      character(len=*), intent(in) :: needle

      if (.not. allocated(err%ptr)) then
         call test_failed(error, "A failing isodensity callback still produced a cavity")
         return
      end if
      if (index(err%ptr%message, "External LSF evaluation failed") <= 0) then
         call test_failed(error, "Unexpected API error message: "//err%ptr%message)
         return
      end if
      if (index(err%ptr%message, needle) <= 0) then
         call test_failed(error, "Error does not name the callback status: "//err%ptr%message)
         return
      end if

   end subroutine expect_iso_callback_error

   !> A callback that fails on its very first evaluation aborts the build.
   subroutine test_iso_callback_fails_first(error)
      type(error_type), allocatable, intent(out) :: error
      type(vp_error), pointer :: err
      type(c_ptr) :: verror, vmol, vcav
      type(iso_cb_ctx), target :: ctx

      ctx%fail_after = 0
      ctx%fail_status = 13_c_int

      call run_iso_callback_update(err, verror, vmol, vcav, ctx)
      call expect_iso_callback_error(error, err, "status 13")

      ! The callback is entered once and then not again: the failure latches and
      ! the build unwinds instead of grinding through the rest of the grid.
      if (.not. allocated(error)) call check(error, ctx%calls, 1)

      call drop_water_cavity(err, vmol, vcav)

   end subroutine test_iso_callback_fails_first

   !> A callback that fails partway through the projection aborts the build.
   !>
   !> This is the case that exercises the abort machinery: the projection runs
   !> its anchors in an OpenMP worksharing construct, which cannot be branched
   !> out of, so the failure has to travel out on a shared flag and become an
   !> error only after the region closes.  A first-evaluation failure would pass
   !> even if that were broken.
   !>
   !> The callback must also stop being called once it has reported failure, and
   !> a later build with a healthy callback must succeed -- the failure latch is
   !> per-build state, not a permanent poisoning of the cavity.
   subroutine test_iso_callback_fails_mid_loop(error)
      type(error_type), allocatable, intent(out) :: error
      type(vp_error), pointer :: err
      type(c_ptr) :: verror, vmol, vcav
      type(iso_cb_ctx), target :: ctx
      integer :: calls_at_abort
      integer(c_int) :: ngrid, nsph
      integer, parameter :: fail_after = 50

      ctx%fail_after = fail_after
      ctx%fail_status = 7_c_int

      call run_iso_callback_update(err, verror, vmol, vcav, ctx)
      calls_at_abort = ctx%calls
      call expect_iso_callback_error(error, err, "status 7")

      if (.not. allocated(error)) then
         if (calls_at_abort <= fail_after) then
            call test_failed(error, "The callback's failing branch was never reached")
         end if
      end if

      ! Once the failure is latched the host callback is not entered again, so
      ! the count stays near the point of failure instead of covering the grid.
      if (.not. allocated(error)) then
         if (calls_at_abort > fail_after + 4096) then
            call test_failed(error, "The failing callback kept being called after it reported failure")
         end if
      end if

      ! Recovery: the same cavity handle rebuilds cleanly once the callback
      ! stops failing.
      if (.not. allocated(error)) then
         deallocate (err%ptr)
         ctx%calls = 0
         ctx%fail_after = huge(1)
         call moist_update_cavity(verror, vcav, vmol)
         if (allocated(err%ptr)) then
            call test_failed(error, "Cavity stayed broken after the callback recovered: " &
                             //err%ptr%message)
         else
            call moist_get_cavity_sizes(verror, vcav, ngrid, nsph)
            if (.not. allocated(error)) call check(error, nsph, 3_c_int)
            if (.not. allocated(error)) call check(error, ngrid > 0)
         end if
      end if

      call drop_water_cavity(err, vmol, vcav)

   end subroutine test_iso_callback_fails_mid_loop

end module test_api
