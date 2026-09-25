!> Defensive-preamble tests for the C API entry points
module test_api
   use iso_c_binding, only: c_ptr, c_loc, c_null_ptr, c_int, c_double, c_bool, &
                            c_funptr, c_funloc, c_associated, c_f_pointer, &
                            c_char, c_null_char, c_size_t, c_sizeof
   use mctc_env_error, only: moist_error_type => error_type
   use moist_api, only: vp_cavity, vp_error, vp_response, response_get_api, &
      & next_response_item_api, response_item_name_api
   use moist_channels_response, only: density_response_type, &
      & potential_adjoint_response_type, response_accumulate
   use testdrive, only: new_unittest, unittest_type, error_type, check, test_failed
   implicit none
   private

   public :: collect_api

   !> Frozen isodensity options layout exercised through the C interface
   type, bind(C) :: iso_options
      !> Caller layout size
      integer(c_size_t) :: struct_size
      !> Density isovalue
      real(c_double) :: rho_iso
      !> Level-set scale
      real(c_double) :: scale
   end type iso_options

   !> C bindings of the entry points under test
   interface

      !> C binding for moist_new_svdw_lsf
      !>
      !> @param[in] error Native error
      !> @param[in] options Native options
      function moist_new_svdw_lsf(error, options) result(lsf) bind(C)
         import :: c_ptr
         implicit none(type, external)
         !> Native error
         type(c_ptr), value, intent(in) :: error
         !> Native options
         type(c_ptr), value, intent(in) :: options
         !> Owning result handle
         type(c_ptr) :: lsf
      end function moist_new_svdw_lsf

      !> C binding for moist_new_isodensity_callback_lsf
      !>
      !> @param[in] error Native error
      !> @param[in] callback Native callback
      !> @param[in] context Native context
      !> @param[in] options Native options
      function moist_new_isodensity_callback_lsf(error, callback, context, options) result(lsf) bind(C)
         import :: c_ptr, c_funptr
         implicit none(type, external)
         !> Native error
         type(c_ptr), value, intent(in) :: error
         !> Native callback
         type(c_funptr), value, intent(in) :: callback
         !> Native context
         type(c_ptr), value, intent(in) :: context
         !> Native options
         type(c_ptr), value, intent(in) :: options
         !> Owning result handle
         type(c_ptr) :: lsf
      end function moist_new_isodensity_callback_lsf

      !> C binding for moist_new_drop_cavity
      !>
      !> @param[in] error Native error
      !> @param[in] lsf Native lsf
      !> @param[in] radii Native radii
      !> @param[in] options Native options
      function moist_new_drop_cavity(error, lsf, radii, options) result(cavity) bind(C)
         import :: c_ptr
         implicit none(type, external)
         !> Native error
         type(c_ptr), value, intent(in) :: error
         !> Native lsf
         type(c_ptr), value, intent(in) :: lsf
         !> Native radii
         type(c_ptr), value, intent(in) :: radii
         !> Native options
         type(c_ptr), value, intent(in) :: options
         !> Owning result handle
         type(c_ptr) :: cavity
      end function moist_new_drop_cavity

      !> C binding for moist_delete_lsf
      !>
      !> @param[in,out] lsf Native lsf
      subroutine moist_delete_lsf(lsf) bind(C)
         import :: c_ptr
         implicit none(type, external)
         !> Native lsf
         type(c_ptr), intent(inout) :: lsf
      end subroutine moist_delete_lsf

      subroutine moist_get_cavity_sizes(verror, vcav, ngrid, nsph) &
            & bind(C, name="moist_get_cavity_sizes")
         import :: c_ptr, c_int
         type(c_ptr), value :: verror
         type(c_ptr), value :: vcav
         integer(c_int), intent(inout) :: ngrid
         integer(c_int), intent(inout) :: nsph
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
         real(c_double), intent(inout) :: area
         real(c_double), intent(inout) :: volume
         integer(c_int), intent(inout) :: ngrid
         integer(c_int), intent(inout) :: nsph
         real(c_double), intent(inout) :: xyz(3, ngrid_cap)
         real(c_double), intent(inout) :: a(ngrid_cap)
         integer(c_int), intent(inout) :: owner(ngrid_cap)
         logical(c_bool), intent(inout) :: converged(ngrid_cap)
         real(c_double), intent(inout) :: radii(nsph_cap)
         real(c_double), intent(inout) :: asph(nsph_cap)
      end subroutine moist_get_cavity_results

      subroutine moist_get_cavity_field_count(verror, vcav, nfield) &
            & bind(C, name="moist_get_cavity_field_count")
         import :: c_ptr, c_int
         type(c_ptr), value :: verror
         type(c_ptr), value :: vcav
         integer(c_int), intent(inout) :: nfield
      end subroutine moist_get_cavity_field_count

      subroutine moist_get_cavity_field_info(verror, vcav, ifield, name, &
            & dtype, rank, dims, count) &
            & bind(C, name="moist_get_cavity_field_info")
         import :: c_ptr, c_int, c_char
         type(c_ptr), value :: verror
         type(c_ptr), value :: vcav
         integer(c_int), value :: ifield
         character(kind=c_char), intent(inout) :: name(*)
         integer(c_int), intent(inout) :: dtype
         integer(c_int), intent(inout) :: rank
         integer(c_int), intent(inout) :: dims(2)
         integer(c_int), intent(inout) :: count
      end subroutine moist_get_cavity_field_info

      subroutine moist_get_cavity_field_real(verror, vcav, cname, values) &
            & bind(C, name="moist_get_cavity_field_real")
         import :: c_ptr, c_double
         type(c_ptr), value :: verror
         type(c_ptr), value :: vcav
         type(c_ptr), value :: cname
         real(c_double), intent(inout) :: values(*)
      end subroutine moist_get_cavity_field_real

      subroutine moist_get_cavity_field_int(verror, vcav, cname, values) &
            & bind(C, name="moist_get_cavity_field_int")
         import :: c_ptr, c_int
         type(c_ptr), value :: verror
         type(c_ptr), value :: vcav
         type(c_ptr), value :: cname
         integer(c_int), intent(inout) :: values(*)
      end subroutine moist_get_cavity_field_int

      subroutine moist_get_cavity_field_bool(verror, vcav, cname, values) &
            & bind(C, name="moist_get_cavity_field_bool")
         import :: c_ptr, c_bool
         type(c_ptr), value :: verror
         type(c_ptr), value :: vcav
         type(c_ptr), value :: cname
         logical(c_bool), intent(inout) :: values(*)
      end subroutine moist_get_cavity_field_bool

      subroutine moist_get_cavity_field_about(verror, vcav, cname, about, capacity, length) &
            & bind(C, name="moist_get_cavity_field_about")
         import :: c_ptr, c_size_t, c_char
         implicit none(type, external)
         type(c_ptr), value, intent(in) :: verror
         type(c_ptr), value, intent(in) :: vcav
         type(c_ptr), value, intent(in) :: cname
         character(kind=c_char), intent(inout), optional :: about(*)
         integer(c_size_t), value, intent(in) :: capacity
         integer(c_size_t), intent(inout), optional :: length
      end subroutine moist_get_cavity_field_about

      subroutine moist_assemble_amat(verror, vcav, ngrid_cap, amat0, xi) &
            & bind(C, name="moist_assemble_amat")
         import :: c_ptr, c_int, c_double
         type(c_ptr), value :: verror
         type(c_ptr), value :: vcav
         integer(c_int), value :: ngrid_cap
         real(c_double), intent(inout) :: amat0(ngrid_cap, ngrid_cap)
         real(c_double), intent(inout) :: xi(ngrid_cap)
      end subroutine moist_assemble_amat

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

   !> Per-test callback state
   !>
   !> Lives in the callback's own `context` rather than in module variables:
   !> test-drive runs the tests of a suite in an OpenMP parallel loop, so two
   !> tests using the same callback would otherwise share one counter
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

   !> Collect all C API guard tests
   subroutine collect_api(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("response_arrays", test_response_arrays), &
                  new_unittest("response_walk", test_response_walk), &
                  new_unittest("update_drop_cavity_null_handle", test_update_drop_null_handle), &
                  new_unittest("update_drop_cavity_null_inner", test_update_drop_null_inner), &
                  new_unittest("get_drop_results_null_handle", test_drop_results_null_handle), &
                  new_unittest("get_drop_results_null_inner", test_drop_results_null_inner), &
                  new_unittest("get_cavity_sizes_null_handle", test_cavity_sizes_null_handle), &
                  new_unittest("get_cavity_sizes_null_inner", test_cavity_sizes_null_inner), &
                  new_unittest("cavity_field_null_name", test_cavity_field_null_name), &
                  new_unittest("cavity_field_null_handle", test_cavity_field_null_handle), &
                  new_unittest("cavity_field_null_inner", test_cavity_field_null_inner), &
                  new_unittest("cavity_field_empty_name", test_cavity_field_empty_name), &
                  new_unittest("cavity_field_index_out_of_range", test_cavity_field_index_range), &
                  new_unittest("cavity_field_about_length_and_truncation", test_cavity_field_about_buffer), &
                  new_unittest("cavity_field_bool_accessor", test_cavity_field_bool), &
                  new_unittest("array_capacity_too_small", test_capacity_too_small), &
                  new_unittest("array_capacity_oversized", test_capacity_oversized), &
                  new_unittest("isodensity_callback_fails_first_call", test_iso_callback_fails_first), &
                  new_unittest("isodensity_callback_fails_mid_loop", test_iso_callback_fails_mid_loop) &
                  ]

   end subroutine collect_api

   !> Assert that the API error handle carries a message containing `needle`
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

   !> A missing DROP cavity handle is reported, not dereferenced
   subroutine test_update_drop_null_handle(error)
      type(error_type), allocatable, intent(out) :: error
      type(vp_error), pointer :: err

      allocate (err)
      call moist_update_cavity(c_loc(err), c_null_ptr, c_null_ptr)
      call check_api_error(error, err, "Cavity handle is missing")
      deallocate (err)

   end subroutine test_update_drop_null_handle

   !> A DROP cavity handle with a null inner pointer is reported, not dereferenced
   subroutine test_update_drop_null_inner(error)
      type(error_type), allocatable, intent(out) :: error
      type(vp_error), pointer :: err
      type(vp_cavity), pointer :: cav

      allocate (err)
      allocate (cav)
      call moist_update_cavity(c_loc(err), c_loc(cav), c_null_ptr)
      call check_api_error(error, err, "Cavity is not initialized")
      deallocate (cav)
      deallocate (err)

   end subroutine test_update_drop_null_inner

   !> A missing DROP cavity handle is reported, not dereferenced
   subroutine test_drop_results_null_handle(error)
      type(error_type), allocatable, intent(out) :: error
      type(vp_error), pointer :: err

      allocate (err)
      call call_get_drop_results(c_loc(err), c_null_ptr)
      call check_api_error(error, err, "Cavity handle is missing")
      deallocate (err)

   end subroutine test_drop_results_null_handle

   !> A DROP cavity handle with a null inner pointer is reported, not dereferenced
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

   !> A missing cavity handle is reported, not dereferenced
   subroutine test_cavity_sizes_null_handle(error)
      type(error_type), allocatable, intent(out) :: error
      type(vp_error), pointer :: err
      integer(c_int) :: ngrid, nsph

      allocate (err)
      call moist_get_cavity_sizes(c_loc(err), c_null_ptr, ngrid, nsph)
      call check_api_error(error, err, "Cavity handle is missing")
      deallocate (err)

   end subroutine test_cavity_sizes_null_handle

   !> A cavity handle with a null inner pointer is reported, not dereferenced
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

   !> Provide the caller-owned result buffers moist_get_cavity_results writes into
   !> The guarded entry point returns before touching them; they only have to
   !> exist so the call is well formed
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

      call moist_get_cavity_results(verror, vcav, 1_c_int, 1_c_int, &
                                    area, volume, ngrid, nsph, xyz, a, owner, converged, radii, asph)

   end subroutine call_get_drop_results

   !> Build a DROP cavity for water through the public C entry points, exactly
   !> as a host would.  The caller owns the returned handles
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
      block
         type(c_ptr) :: lsf
         lsf = moist_new_svdw_lsf(verror, c_null_ptr)
         vcav = moist_new_drop_cavity(verror, lsf, c_null_ptr, c_null_ptr)
         call moist_delete_lsf(lsf)
      end block
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

   !> Release the handles built by build_water_cavity
   subroutine drop_water_cavity(err, vmol, vcav)
      type(vp_error), pointer, intent(inout) :: err
      type(c_ptr), intent(inout) :: vmol, vcav

      call moist_delete_cavity(vcav)
      call moist_delete_structure(vmol)
      deallocate (err)

   end subroutine drop_water_cavity

   !> Null-terminate a Fortran string into a buffer c_loc can be taken of
   function c_string(text) result(buffer)
      !> Text to hand to a C entry point
      character(len=*), intent(in) :: text
      !> Null-terminated copy
      character(kind=c_char), allocatable :: buffer(:)
      integer :: ichar

      allocate (buffer(len(text) + 1))
      do ichar = 1, len(text)
         buffer(ichar) = text(ichar:ichar)
      end do
      buffer(len(text) + 1) = c_null_char

   end function c_string

   !> A null field name is reported rather than dereferenced
   !>
   !> The name crosses as a bare `const char*` that the entry point has to scan,
   !> so a caller that passes NULL must get an error instead of a walk off
   !> address zero
   subroutine test_cavity_field_null_name(error)
      type(error_type), allocatable, intent(out) :: error
      type(vp_error), pointer :: err
      type(c_ptr) :: verror, vmol, vcav
      integer(c_int) :: ngrid, nsph
      real(c_double) :: values(1)

      call build_water_cavity(error, err, verror, vmol, vcav, ngrid, nsph)
      if (allocated(error)) return

      values = -12345.0_c_double
      call moist_get_cavity_field_real(verror, vcav, c_null_ptr, values)
      call check_api_error(error, err, "Field name is missing")

      if (.not. allocated(error)) then
         if (any(values /= -12345.0_c_double)) &
            call test_failed(error, "get_cavity_field_real wrote for a null name")
      end if

      call drop_water_cavity(err, vmol, vcav)

   end subroutine test_cavity_field_null_name

   !* ================================================================================= *!
   !*                          Named field API defensive guards                         *!
   !* ================================================================================= *!

   !> Position of the null terminator in a C string buffer, i.e. its length
   pure integer function c_string_len(buffer) result(nchar)
      !> Null-terminated buffer written by a C entry point
      character(kind=c_char), intent(in) :: buffer(:)
      integer :: ichar

      nchar = size(buffer)
      do ichar = 1, size(buffer)
         if (buffer(ichar) == c_null_char) then
            nchar = ichar - 1
            return
         end if
      end do

   end function c_string_len

   !> A missing cavity handle is reported by the field API, not dereferenced
   !>
   !> `resolve_field_cavity` is the shared preamble of all six field entry
   !> points, so this single guard is what stands between every named read and a
   !> walk off a null handle.  Both an enumerating and a fetching entry point
   !> are exercised because they reach the guard by different routes
   subroutine test_cavity_field_null_handle(error)
      type(error_type), allocatable, intent(out) :: error
      type(vp_error), pointer :: err
      character(kind=c_char), allocatable, target :: name_wleb(:)
      integer(c_int) :: nfield
      real(c_double) :: values(1)

      allocate (err)
      name_wleb = c_string("wleb")

      nfield = -1_c_int
      call moist_get_cavity_field_count(c_loc(err), c_null_ptr, nfield)
      call check_api_error(error, err, "Cavity handle is missing")
      if (.not. allocated(error)) then
         if (nfield /= -1_c_int) &
            call test_failed(error, "get_cavity_field_count changed nfield on failure")
      end if
      if (allocated(err%ptr)) deallocate (err%ptr)

      if (.not. allocated(error)) then
         values = -12345.0_c_double
         call moist_get_cavity_field_real(c_loc(err), c_null_ptr, c_loc(name_wleb), values)
         call check_api_error(error, err, "Cavity handle is missing")
      end if
      if (.not. allocated(error)) then
         if (any(values /= -12345.0_c_double)) &
            call test_failed(error, "get_cavity_field_real wrote for a missing handle")
      end if

      deallocate (err)

   end subroutine test_cavity_field_null_handle

   !> A cavity handle with a null inner pointer is reported, not dereferenced
   !>
   !> This is the state a host holds between `moist_new_*_cavity` and the first
   !> `moist_update_cavity`: the handle is real, the cavity behind it is not
   subroutine test_cavity_field_null_inner(error)
      type(error_type), allocatable, intent(out) :: error
      type(vp_error), pointer :: err
      type(vp_cavity), pointer :: cav
      character(kind=c_char), allocatable, target :: name_wleb(:)
      integer(c_int) :: nfield
      real(c_double) :: values(1)

      allocate (err)
      allocate (cav)
      name_wleb = c_string("wleb")

      nfield = -1_c_int
      call moist_get_cavity_field_count(c_loc(err), c_loc(cav), nfield)
      call check_api_error(error, err, "Cavity is not initialized")
      if (.not. allocated(error)) then
         if (nfield /= -1_c_int) &
            call test_failed(error, "get_cavity_field_count changed nfield on failure")
      end if
      if (allocated(err%ptr)) deallocate (err%ptr)

      if (.not. allocated(error)) then
         values = -12345.0_c_double
         call moist_get_cavity_field_real(c_loc(err), c_loc(cav), c_loc(name_wleb), values)
         call check_api_error(error, err, "Cavity is not initialized")
      end if
      if (.not. allocated(error)) then
         if (any(values /= -12345.0_c_double)) &
            call test_failed(error, "get_cavity_field_real wrote for an empty cavity")
      end if

      deallocate (cav)
      deallocate (err)

   end subroutine test_cavity_field_null_inner

   !> An empty field name is rejected rather than matched against a declaration
   !>
   !> A bare `""` survives the null-pointer guard and reaches the name scan, so
   !> without its own check it would be compared against every declared name --
   !> harmless today, but only because no cavity declares an empty name
   subroutine test_cavity_field_empty_name(error)
      type(error_type), allocatable, intent(out) :: error
      type(vp_error), pointer :: err
      type(c_ptr) :: verror, vmol, vcav
      integer(c_int) :: ngrid, nsph
      character(kind=c_char), allocatable, target :: name_empty(:)
      real(c_double) :: values(1)

      call build_water_cavity(error, err, verror, vmol, vcav, ngrid, nsph)
      if (allocated(error)) return

      name_empty = c_string("")
      values = -12345.0_c_double
      call moist_get_cavity_field_real(verror, vcav, c_loc(name_empty), values)
      call check_api_error(error, err, "Field name is empty")

      if (.not. allocated(error)) then
         if (any(values /= -12345.0_c_double)) &
            call test_failed(error, "get_cavity_field_real wrote for an empty name")
      end if

      call drop_water_cavity(err, vmol, vcav)

   end subroutine test_cavity_field_empty_name

   !> A field index outside the reported count is refused, and the descriptor
   !> outputs remain unchanged
   !>
   !> The index is the one field-API argument the caller derives from an earlier
   !> call, so a stale count is the realistic way to get here.  Unchecked it
   !> would index `query%info` at 0 or past `nfield`; the cleared outputs matter
   !> because a caller that ignores the error must not read a plausible-looking
   !> descriptor out of the buffers
   subroutine test_cavity_field_index_range(error)
      type(error_type), allocatable, intent(out) :: error
      type(vp_error), pointer :: err
      type(c_ptr) :: verror, vmol, vcav
      integer(c_int) :: ngrid, nsph, nfield
      integer(c_int) :: dtype, rank, dims(2), count
      character(kind=c_char), target :: name(64)

      call build_water_cavity(error, err, verror, vmol, vcav, ngrid, nsph)
      if (allocated(error)) return

      call moist_get_cavity_field_count(verror, vcav, nfield)
      if (allocated(err%ptr)) then
         call test_failed(error, "get_cavity_field_count failed: "//err%ptr%message)
      else if (nfield <= 0_c_int) then
         call test_failed(error, "A built cavity declares no fields")
      end if

      if (.not. allocated(error)) then
         call probe_field_info(verror, vcav, -1_c_int, name, dtype, rank, dims, count)
         call check_api_error(error, err, "Field index out of range")
      end if
      if (.not. allocated(error)) &
         call check_unchanged_info(error, "negative index", name, dtype, rank, dims, count)
      if (allocated(err%ptr)) deallocate (err%ptr)

      if (.not. allocated(error)) then
         call probe_field_info(verror, vcav, nfield, name, dtype, rank, dims, count)
         call check_api_error(error, err, "Field index out of range")
      end if
      if (.not. allocated(error)) &
         call check_unchanged_info(error, "one past the end", name, dtype, rank, dims, count)

      call drop_water_cavity(err, vmol, vcav)

   end subroutine test_cavity_field_index_range

   !> Seed the descriptor buffers with detectable values, then describe a field
   subroutine probe_field_info(verror, vcav, ifield, name, dtype, rank, dims, count)
      !> Opaque error and cavity handles
      type(c_ptr), intent(in) :: verror, vcav
      !> Field position to describe
      integer(c_int), intent(in) :: ifield
      !> Name buffer, seeded so any write is visible
      character(kind=c_char), intent(inout) :: name(:)
      !> Descriptor outputs, seeded so a cleared value is distinguishable
      integer(c_int), intent(out) :: dtype, rank, dims(2), count

      name = "Z"
      dtype = -1_c_int
      rank = -1_c_int
      dims = -1_c_int
      count = -1_c_int
      call moist_get_cavity_field_info(verror, vcav, ifield, name, dtype, rank, dims, count)

   end subroutine probe_field_info

   !> Assert a rejected describe left the descriptor and name unchanged
   subroutine check_unchanged_info(error, what, name, dtype, rank, dims, count)
      type(error_type), allocatable, intent(out) :: error
      !> Case being checked, used in the failure message
      character(len=*), intent(in) :: what
      !> Name buffer that must still carry its seed
      character(kind=c_char), intent(in) :: name(:)
      !> Descriptor outputs that must retain their sentinel values
      integer(c_int), intent(in) :: dtype, rank, dims(2), count

      if (dtype /= -1_c_int .or. rank /= -1_c_int .or. count /= -1_c_int) then
         call test_failed(error, "get_cavity_field_info changed the descriptor for "//what)
         return
      end if
      if (any(dims /= -1_c_int)) then
         call test_failed(error, "get_cavity_field_info changed extents for "//what)
         return
      end if
      if (any(name /= "Z")) &
         call test_failed(error, "get_cavity_field_info wrote a name for "//what)

   end subroutine check_unchanged_info

   !> Description queries report full length while copying only what fits
   subroutine test_cavity_field_about_buffer(error)
      !> Test diagnostic
      type(error_type), allocatable, intent(out) :: error
      type(vp_error), pointer :: err
      type(c_ptr) :: verror, vmol, vcav
      integer(c_int) :: ngrid, nsph
      character(kind=c_char), allocatable, target :: name_ngrid(:)
      character(kind=c_char), allocatable :: about(:)
      integer(c_size_t) :: nchar, length

      call build_water_cavity(error, err, verror, vmol, vcav, ngrid, nsph)
      if (allocated(error)) return
      name_ngrid = c_string("ngrid")
      description_checks: block
         call moist_get_cavity_field_about(verror, vcav, c_loc(name_ngrid), &
            & capacity=0_c_size_t, length=nchar)
         call check(error, .not. allocated(err%ptr), "Length-only query failed")
         if (allocated(error)) exit description_checks
         call check(error, nchar > 0, "Expected a nonempty description")
         if (allocated(error)) exit description_checks
         allocate(about(nchar + 1))
         about = "Z"
         call moist_get_cavity_field_about(verror, vcav, c_loc(name_ngrid), about, nchar, length)
         call check(error, .not. allocated(err%ptr), "Truncation should succeed")
         if (allocated(error)) exit description_checks
         call check(error, length == nchar .and. int(c_string_len(about), c_size_t) == nchar - 1)
         if (allocated(error)) exit description_checks
         call check(error, about(nchar + 1) == "Z", "Wrote past supplied capacity")
         if (allocated(error)) exit description_checks
         call moist_get_cavity_field_about(verror, vcav, c_loc(name_ngrid), about, nchar + 1, length)
         call check(error, .not. allocated(err%ptr), "Exact allocation should succeed")
         if (allocated(error)) exit description_checks
         call check(error, length == nchar .and. int(c_string_len(about), c_size_t) == nchar)
      end block description_checks
      call drop_water_cavity(err, vmol, vcav)
   end subroutine test_cavity_field_about_buffer

   !> The logical accessor honours the same contract as the other two
   !>
   !> `converged` is the only boolean a cavity declares, so the bool accessor is
   !> the one branch of `check_field_payload` a C host reaches without a Python
   !> layer in between.  Its payload is cross-checked against the independent
   !> `get_cavity_results` path, since every value it can return is also a legal
   !> seed and a buffer left untouched would otherwise pass unnoticed
   subroutine test_cavity_field_bool(error)
      type(error_type), allocatable, intent(out) :: error
      type(vp_error), pointer :: err
      type(c_ptr) :: verror, vmol, vcav
      integer(c_int) :: ngrid, nsph, out_ngrid, out_nsph
      real(c_double) :: area, volume
      real(c_double), allocatable :: xyz(:, :), a(:), radii(:), asph(:)
      integer(c_int), allocatable :: owner(:)
      logical(c_bool), allocatable :: reference(:), values(:)
      character(kind=c_char), allocatable, target :: name_conv(:), name_wleb(:)

      call build_water_cavity(error, err, verror, vmol, vcav, ngrid, nsph)
      if (allocated(error)) return

      name_conv = c_string("converged")
      name_wleb = c_string("wleb")
      allocate (xyz(3, ngrid), a(ngrid), owner(ngrid), reference(ngrid))
      allocate (radii(nsph), asph(nsph), values(ngrid))

      call moist_get_cavity_results(verror, vcav, ngrid, nsph, area, volume, &
                                    out_ngrid, out_nsph, xyz, a, owner, reference, &
                                    radii, asph)
      if (allocated(err%ptr)) &
         call test_failed(error, "get_cavity_results failed: "//err%ptr%message)

      !> The declared field has to agree with the fixed getter, point for point
      if (.not. allocated(error)) then
         values = .not. reference
         call moist_get_cavity_field_bool(verror, vcav, c_loc(name_conv), values)
         if (allocated(err%ptr)) then
            call test_failed(error, "get_cavity_field_bool failed: "//err%ptr%message)
         else if (.not. all(values .eqv. reference)) then
            call test_failed(error, "get_cavity_field_bool disagrees with get_cavity_results")
         end if
      end if

      !> Reading a real field through the logical accessor is refused too
      if (.not. allocated(error)) then
         values = .false._c_bool
         call moist_get_cavity_field_bool(verror, vcav, c_loc(name_wleb), values)
         call check_api_error(error, err, "has a different element type")
      end if
      if (.not. allocated(error)) then
         if (any(values)) &
            call test_failed(error, "get_cavity_field_bool wrote for a real-valued field")
      end if

      call drop_water_cavity(err, vmol, vcav)

   end subroutine test_cavity_field_bool

   !> Assert that an entry point rejected an undersized capacity, then clear the
   !> error so the handle can be reused
   subroutine expect_capacity_error(error, err, routine)
      type(error_type), allocatable, intent(out) :: error
      type(vp_error), pointer, intent(inout) :: err
      character(len=*), intent(in) :: routine

      if (.not. allocated(err%ptr)) then
         call test_failed(error, routine//" accepted an undersized array capacity")
         return
      end if
      if (index(err%ptr%message, "Array capacity") <= 0) then
         call test_failed(error, "Unexpected API error message: "//err%ptr%message)
         return
      end if
      deallocate (err%ptr)

   end subroutine expect_capacity_error

   !> A capacity below the cavity's own ngrid/nsph is refused before a single
   !> element is written.  Without the check these calls write ngrid values into
   !> buffers holding ngrid-1, which is silent heap corruption; the sentinels
   !> below detect any write at all
   subroutine test_capacity_too_small(error)
      type(error_type), allocatable, intent(out) :: error
      type(vp_error), pointer :: err
      type(c_ptr) :: verror, vmol, vcav
      integer(c_int) :: ngrid, nsph, cap
      integer(c_int) :: out_ngrid, out_nsph
      real(c_double) :: area, volume
      real(c_double), allocatable :: xyz(:, :), a(:), radii(:), asph(:)
      real(c_double), allocatable :: amat0(:, :), xi(:)
      integer(c_int), allocatable :: owner(:)
      logical(c_bool), allocatable :: converged(:)
      real(c_double), parameter :: sentinel = -12345.0_c_double

      call build_water_cavity(error, err, verror, vmol, vcav, ngrid, nsph)
      if (allocated(error)) return

      cap = ngrid - 1_c_int
      allocate (xyz(3, cap), a(cap), owner(cap), converged(cap))
      allocate (amat0(cap, cap), xi(cap))
      allocate (radii(nsph), asph(nsph))
      xyz = sentinel
      a = sentinel
      owner = -1_c_int
      converged = .false._c_bool
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

   !> Accept oversized capacities and preserve padding outside the logical extent
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

      call check(error, .not. allocated(err%ptr), "Oversized capacity must be accepted")
      if (.not. allocated(error)) call check(error, out_ngrid, ngrid)
      if (.not. allocated(error)) call check(error, out_nsph, nsph)
      if (.not. allocated(error)) call check(error, all(a(:ngrid) /= sentinel))
      if (.not. allocated(error)) call check(error, all(xyz(:, :ngrid) /= sentinel))
      if (.not. allocated(error)) call check(error, all(radii(:nsph) /= sentinel))
      if (.not. allocated(error)) call check(error, all(a(ngrid + 1:) == sentinel))
      if (.not. allocated(error)) call check(error, all(xyz(:, ngrid + 1:) == sentinel))
      if (.not. allocated(error)) call check(error, all(radii(nsph + 1:) == sentinel))

      call drop_water_cavity(err, vmol, vcav)

   end subroutine test_capacity_oversized

   !* ================================================================================= *!
   !*                        Isodensity callback failure channel                        *!
   !* ================================================================================= *!

   !> Isodensity LSF callback with a switchable failure channel
   !>
   !> Answers the first `ctx%fail_after` evaluations with a smooth three-Gaussian
   !> model density (so the projection behaves like a real build), then reports
   !> `ctx%fail_status` for every evaluation after that.  The counter is bumped
   !> atomically because moist calls this from inside OpenMP loops
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

   !> Build a water cavity backed by `iso_switchable_callback` and update it
   !>
   !> Returns the handles so the caller can inspect the error and delete them
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
      block
         type(c_ptr) :: lsf
         type(iso_options), target :: options
         options = iso_options(c_sizeof(options), cb_rho_iso, 1000.0_c_double)
         lsf = moist_new_isodensity_callback_lsf(verror, c_funloc(iso_switchable_callback), &
                                               c_loc(ctx), c_loc(options))
         vcav = moist_new_drop_cavity(verror, lsf, c_null_ptr, c_null_ptr)
         call moist_delete_lsf(lsf)
      end block
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

   !> A callback that fails on its very first evaluation aborts the build
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
      ! the build unwinds instead of grinding through the rest of the grid
      if (.not. allocated(error)) call check(error, ctx%calls, 1)

      call drop_water_cavity(err, vmol, vcav)

   end subroutine test_iso_callback_fails_first

   !> A callback that fails partway through the projection aborts the build
   !>
   !> This is the case that exercises the abort machinery: the projection runs
   !> its anchors in an OpenMP worksharing construct, which cannot be branched
   !> out of, so the failure has to travel out on a shared flag and become an
   !> error only after the region closes.  A first-evaluation failure would pass
   !> even if that were broken
   !>
   !> The callback must also stop being called once it has reported failure, and
   !> a later build with a healthy callback must succeed -- the failure latch is
   !> per-build state, not a permanent poisoning of the cavity
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
      ! the count stays near the point of failure instead of covering the grid
      if (.not. allocated(error)) then
         if (calls_at_abort > fail_after + 4096) then
            call test_failed(error, "The failing callback kept being called after it reported failure")
         end if
      end if

      ! Recovery: the same cavity handle rebuilds cleanly once the callback
      ! stops failing
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

   !> Response copies write exactly the logical shape at all ranks, leave the
   !> rest of a larger buffer untouched, and name a missing current item, an
   !> array the current item does not have and a NULL buffer
   !>
   !> @param[out] error Test failure
   subroutine test_response_arrays(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Native wrappers
      type(vp_error), target :: err
      type(vp_response), target :: response
      !> Distinct values expose layout and truncation errors
      type(density_response_type) :: density
      !> Two logical points in four-point buffers
      real(c_double), target :: rho(4), grad(3, 4), hess(3, 3, 4)
      !> NUL-terminated array names
      character(kind=c_char), allocatable, target :: w_rho(:), w_grad_rho(:), &
         & w_hess_rho(:), w_phi(:), bogus(:)
      !> Padding marker
      real(c_double), parameter :: sentinel = -12345.0_c_double
      !> Value index
      integer :: i

      allocate (density%w_rho(2), density%w_grad_rho(3, 2), density%w_hess_rho(3, 3, 2))
      density%w_rho = [1.0_c_double, 2.0_c_double]
      density%w_grad_rho = reshape([(real(i, c_double), i=3, 8)], [3, 2])
      density%w_hess_rho = reshape([(real(i, c_double), i=9, 26)], [3, 3, 2])
      call response_accumulate(response%ptr, density, err%ptr)
      call check(error, .not. allocated(err%ptr))
      if (allocated(error)) return
      w_rho = c_string("w_rho")
      w_grad_rho = c_string("w_grad_rho")
      w_hess_rho = c_string("w_hess_rho")
      w_phi = c_string("w_phi")
      bogus = c_string("w_bogus")

      ! Arrays are read from the item the walk stopped at, and none is current yet
      call response_get_api(c_loc(err), c_loc(response), c_loc(w_rho), c_loc(rho))
      call check_api_error(error, err, "[moist_get_response_array] No current response item - call next() first")
      if (allocated(error)) return
      call check(error, logical(next_response_item_api(c_loc(err), c_loc(response))))
      if (allocated(error)) return

      ! The logical points are written and the rest of the buffer is left alone
      rho = sentinel
      grad = sentinel
      hess = sentinel
      call response_get_api(c_loc(err), c_loc(response), c_loc(w_rho), c_loc(rho))
      call check(error, .not. allocated(err%ptr))
      if (allocated(error)) return
      call response_get_api(c_loc(err), c_loc(response), c_loc(w_grad_rho), c_loc(grad))
      call check(error, .not. allocated(err%ptr))
      if (allocated(error)) return
      call response_get_api(c_loc(err), c_loc(response), c_loc(w_hess_rho), c_loc(hess))
      call check(error, .not. allocated(err%ptr))
      if (allocated(error)) return
      call check(error, all(rho(:2) == density%w_rho))
      if (allocated(error)) return
      call check(error, all(grad(:, :2) == density%w_grad_rho))
      if (allocated(error)) return
      call check(error, all(hess(:, :, :2) == density%w_hess_rho))
      if (allocated(error)) return
      call check(error, all(rho(3:) == sentinel) .and. all(grad(:, 3:) == sentinel) &
         & .and. all(hess(:, :, 3:) == sentinel))
      if (allocated(error)) return

      ! An array the current item does not have and a NULL buffer are named
      call response_get_api(c_loc(err), c_loc(response), c_loc(bogus), c_loc(rho))
      call check_api_error(error, err, "[moist_get_response_array] density has no array 'w_bogus'")
      if (allocated(error)) return
      call response_get_api(c_loc(err), c_loc(response), c_loc(w_phi), c_loc(rho))
      call check_api_error(error, err, "[moist_get_response_array] density has no array 'w_phi'")
      if (allocated(error)) return
      call response_get_api(c_loc(err), c_loc(response), c_loc(w_grad_rho), c_null_ptr)
      call check_api_error(error, err, "Null array pointer provided for 'w_grad_rho'")
      if (allocated(error)) return

      ! The pass ends after the only item, and nothing is current any more
      call check(error, .not. logical(next_response_item_api(c_loc(err), c_loc(response))))
      if (allocated(error)) return
      call check(error, .not. allocated(err%ptr))
      if (allocated(error)) return
      call response_get_api(c_loc(err), c_loc(response), c_loc(w_rho), c_loc(rho))
      call check_api_error(error, err, "No current response item - call next() first")
   end subroutine test_response_arrays

   !> The C walk visits every item once per pass in accumulation order, rewinds
   !> after the last one, names the current item and fails by name outside a pass
   !>
   !> @param[out] error Test failure
   subroutine test_response_walk(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Native wrappers
      type(vp_error), target :: err
      type(vp_response), target :: response
      !> Items accumulated in this order
      type(potential_adjoint_response_type) :: adjoint
      type(density_response_type) :: density
      !> Pass index
      integer :: pass

      ! A missing response handle ends the loop with an error
      call check(error, .not. logical(next_response_item_api(c_loc(err), c_null_ptr)))
      if (allocated(error)) return
      call check_api_error(error, err, "[moist_next_response_item] Response handle is missing")
      if (allocated(error)) return
      ! A response nothing has filled has nothing to visit, which is no error
      call check(error, .not. logical(next_response_item_api(c_loc(err), c_loc(response))))
      if (allocated(error)) return
      call check(error, .not. allocated(err%ptr))
      if (allocated(error)) return
      call check(error, current_item_name(err, response), "")
      if (allocated(error)) return
      call check_api_error(error, err, &
         & "[moist_get_response_item_name] No current response item - call next() first")
      if (allocated(error)) return

      adjoint%w_phi = [1.0_c_double, 2.0_c_double]
      allocate (density%w_rho(2), source=3.0_c_double)
      allocate (density%w_grad_rho(3, 2), source=4.0_c_double)
      allocate (density%w_hess_rho(3, 3, 2), source=5.0_c_double)
      call response_accumulate(response%ptr, adjoint, err%ptr)
      call response_accumulate(response%ptr, density, err%ptr)
      call check(error, .not. allocated(err%ptr))
      if (allocated(error)) return

      ! Two passes: the first ends by rewinding, so the second starts over
      do pass = 1, 2
         call check(error, logical(next_response_item_api(c_loc(err), c_loc(response))))
         if (allocated(error)) return
         call check(error, current_item_name(err, response), "potential_adjoint")
         if (allocated(error)) return
         call check(error, logical(next_response_item_api(c_loc(err), c_loc(response))))
         if (allocated(error)) return
         call check(error, current_item_name(err, response), "density")
         if (allocated(error)) return
         call check(error, .not. logical(next_response_item_api(c_loc(err), c_loc(response))), &
            & more="the pass did not end after the last item")
         if (allocated(error)) return
         call check(error, .not. allocated(err%ptr))
         if (allocated(error)) return
      end do

      ! A missing name buffer is refused
      call check(error, logical(next_response_item_api(c_loc(err), c_loc(response))))
      if (allocated(error)) return
      call response_item_name_api(c_loc(err), c_loc(response))
      call check_api_error(error, err, "[moist_get_response_item_name] Output pointer is missing")
   end subroutine test_response_walk

   !> Name of the current response item through the C entry point, empty on error
   function current_item_name(err, response) result(name)
      !> Error wrapper receiving a failure
      type(vp_error), target, intent(inout) :: err
      !> Response wrapper under test
      type(vp_response), target, intent(inout) :: response
      !> Decoded name
      character(len=:), allocatable :: name
      !> MOIST_NAME_MAX + 1 fits every name
      character(kind=c_char) :: buffer(33)
      integer :: ichar

      buffer = c_null_char
      call response_item_name_api(c_loc(err), c_loc(response), buffer)
      name = ""
      do ichar = 1, c_string_len(buffer)
         name = name//buffer(ichar)
      end do
   end function current_item_name


end module test_api
