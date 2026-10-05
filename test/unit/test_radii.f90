module test_radii
   use moist_cavity_drop_lsf_svdw_param, only: moist_cavity_drop_lsf_svdw_param_type
   use moist_cavity_drop_parameters, only: moist_cavity_drop_parameters_type
   use mctc_env, only: wp
   use mctc_env_error, only: moist_error_type => error_type
   use mctc_io, only: structure_type
   use mstore, only: get_structure
   use moist_cavity_drop, only: cavity_type_drop, new_cavity_drop
   use moist_cavity_drop_lsf_svdw, only: moist_cavity_drop_lsf_svdw_type
   use moist_data_radii_legacy, only: get_radius_func, rad_type
   use moist_radii, only: radius_type, radius_type_static, radius_type_custom, default_cpcm_radii
   use moist_radii_static, only: new_rahm_radii, new_gauss_radii
   use moist_radii, only: new_cpcm_radii, new_smd_radii, new_d3_radii
   use moist_radii, only: new_cosmo_radii, new_bondi_radii
   use moist_radii, only: new_radii, new_radii_custom_atoms, new_radii_custom_elements
   use testdrive, only: new_unittest, unittest_type, error_type, check, test_failed
   use moist_context, only: moist_context_type, new_context
   implicit none(type, external)
   private

   public :: collect_radii

   real(wp), parameter :: thr = 10*epsilon(1.0_wp)

contains

   !> Collect all radii tests
   subroutine collect_radii(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("radii_factory_contracts", test_factory_contracts), &
                  new_unittest("radii_cache_lifecycle", test_cache_lifecycle), &
                  new_unittest("radii_guard_contracts", test_guard_contracts), &
                  new_unittest("static_radii_cpcm", test_static_radii_cpcm), &
                  new_unittest("static_radii_constructors", test_static_constructors), &
                  new_unittest("radii_constructor_verbosity", test_constructor_verbosity), &
                  new_unittest("static_radii_zero_gradient", test_static_zero_gradient), &
                  new_unittest("static_radii_needs_update", test_static_requires_update), &
                  new_unittest("custom_radii_atoms_works", test_custom_atoms_dropcess), &
                  new_unittest("custom_radii_elements_works", test_custom_elements_dropcess), &
                  new_unittest("custom_radii_cavity_integration", test_custom_radii_cavity_integration), &
                  new_unittest("custom_radii_atoms_bad_empty", test_custom_atoms_empty_fails, should_fail=.true.), &
                  new_unittest("custom_radii_atoms_bad_value", test_custom_atoms_nonpositive_fails, should_fail=.true.), &
                  new_unittest("custom_radii_atoms_bad_nat", test_custom_atoms_size_mismatch_fails, should_fail=.true.), &
                  new_unittest("custom_radii_elements_bad_empty", test_custom_elements_empty_fails, should_fail=.true.), &
                  new_unittest("custom_radii_elements_bad_size", test_custom_elements_size_mismatch_fails, should_fail=.true.), &
                  new_unittest("custom_radii_elements_bad_z", test_custom_elements_invalid_z_fails, should_fail=.true.), &
                  new_unittest("custom_radii_elements_bad_value", test_custom_elements_nonpositive_fails, should_fail=.true.), &
                  new_unittest("custom_radii_elements_bad_duplicate", test_custom_elements_duplicate_fails, should_fail=.true.), &
                  new_unittest("custom_radii_elements_missing", test_custom_elements_missing_for_molecule_fails, &
                               should_fail=.true.), &
                  new_unittest("custom_radii_string_guidance", test_custom_string_guidance, should_fail=.true.) &
                  ]
   end subroutine collect_radii

   !> Fetch a small reference molecule for the radii suite from mstore
   !> AlH3 from MB16-43 has composition [Al, H, H, H] (Z = [13, 1, 1, 1],
   !> nat = 4) - enough atomic diversity to exercise both per-element and
   !> per-atom code paths without being large enough to slow the suite
   !>
   !> @param[out] mol  structure populated from MB16-43/AlH3
   subroutine make_test_molecule(mol)
      !> Structure to populate with AlH3 from MB16-43
      type(structure_type), intent(out) :: mol

      call get_structure(mol, "MB16-43", "AlH3")
   end subroutine make_test_molecule

   subroutine test_static_radii_cpcm(error)
      type(error_type), allocatable, intent(out) :: error

      type(structure_type) :: mol
      type(radius_type_static) :: model
      type(moist_error_type), allocatable :: err
      real(wp) :: ref
      integer :: iat

      call make_test_molecule(mol)
      call new_cpcm_radii(model)

      call model%update(mol, err)
      if (allocated(err)) then
         call test_failed(error, "CPCM static update failed")
         return
      end if

      if (.not. allocated(model%f0)) then
         call test_failed(error, "CPCM static update did not cache f0")
         return
      end if

      do iat = 1, mol%nat
         ref = get_radius_func(mol%num(mol%id(iat)), "cpcm", err)
         if (allocated(err)) then
            call test_failed(error, "CPCM reference lookup failed: "//trim(err%message))
            return
         end if
         call check(error, model%f0(iat), ref, thr=thr, &
                    more="CPCM static radii mismatch")
         if (allocated(error)) return
      end do
   end subroutine test_static_radii_cpcm

   subroutine test_static_constructors(error)
      type(error_type), allocatable, intent(out) :: error

      type(structure_type) :: mol
      type(radius_type_static) :: model
      type(moist_error_type), allocatable :: err
      real(wp) :: ref

      call make_test_molecule(mol)

      call new_smd_radii(model)
      call model%update(mol, err)
      if (allocated(err)) then
         call test_failed(error, "SMD static update failed")
         return
      end if
      if (.not. allocated(model%f0)) then
         call test_failed(error, "SMD static update did not cache f0")
         return
      end if
      ref = get_radius_func(mol%num(mol%id(1)), "smd", err)
      if (allocated(err)) then
         call test_failed(error, "smd reference lookup failed: "//trim(err%message))
         return
      end if
      call check(error, model%f0(1), ref, thr=thr, more="SMD radius mismatch")
      if (allocated(error)) return

      call new_d3_radii(model)
      call model%update(mol, err)
      if (allocated(err)) then
         call test_failed(error, "D3 static update failed")
         return
      end if
      if (.not. allocated(model%f0)) then
         call test_failed(error, "D3 static update did not cache f0")
         return
      end if
      ref = get_radius_func(mol%num(mol%id(1)), "d3", err)
      if (allocated(err)) then
         call test_failed(error, "d3 reference lookup failed: "//trim(err%message))
         return
      end if
      call check(error, model%f0(1), ref, thr=thr, more="D3 radius mismatch")
      if (allocated(error)) return

      call new_cosmo_radii(model)
      call model%update(mol, err)
      if (allocated(err)) then
         call test_failed(error, "COSMO static update failed")
         return
      end if
      if (.not. allocated(model%f0)) then
         call test_failed(error, "COSMO static update did not cache f0")
         return
      end if
      ref = get_radius_func(mol%num(mol%id(1)), "cosmo", err)
      if (allocated(err)) then
         call test_failed(error, "cosmo reference lookup failed: "//trim(err%message))
         return
      end if
      call check(error, model%f0(1), ref, thr=thr, more="COSMO radius mismatch")
      if (allocated(error)) return

      call new_bondi_radii(model)
      call model%update(mol, err)
      if (allocated(err)) then
         call test_failed(error, "Bondi static update failed")
         return
      end if
      if (.not. allocated(model%f0)) then
         call test_failed(error, "Bondi static update did not cache f0")
         return
      end if
      ref = get_radius_func(mol%num(mol%id(1)), "bondi", err)
      if (allocated(err)) then
         call test_failed(error, "bondi reference lookup failed: "//trim(err%message))
         return
      end if
      call check(error, model%f0(1), ref, thr=thr, more="Bondi radius mismatch")
   end subroutine test_static_constructors

   subroutine test_constructor_verbosity(error)
      type(error_type), allocatable, intent(out) :: error

      type(radius_type_static) :: static_model
      class(radius_type), allocatable :: model
      type(moist_error_type), allocatable :: err

      call new_cpcm_radii(static_model)
      call check(error, static_model%verbosity, 0, &
                 more="default static constructor verbosity must be zero")
      if (allocated(error)) return

      call new_cpcm_radii(static_model, verbosity=3)
      call check(error, static_model%verbosity, 3, &
                 more="explicit static constructor verbosity mismatch")
      if (allocated(error)) return

      call new_radii("cpcm", model, err, verbosity=4)
      if (allocated(err)) then
         call test_failed(error, "new_radii with verbosity failed: "//trim(err%message))
         return
      end if
      call check(error, model%verbosity, 4, &
                 more="new_radii string constructor verbosity mismatch")
      if (allocated(error)) return

      call new_radii_custom_atoms([2.10_wp, 1.25_wp, 1.25_wp], model, err, verbosity=5)
      if (allocated(err)) then
         call test_failed(error, "new_radii_custom_atoms with verbosity failed: "//trim(err%message))
         return
      end if
      call check(error, model%verbosity, 5, &
                 more="custom atom constructor verbosity mismatch")
   end subroutine test_constructor_verbosity

   subroutine test_static_zero_gradient(error)
      type(error_type), allocatable, intent(out) :: error

      type(structure_type) :: mol
      type(radius_type_static) :: model
      type(moist_error_type), allocatable :: err

      call make_test_molecule(mol)
      call new_cpcm_radii(model)

      call model%update(mol, err)
      if (allocated(err)) then
         call test_failed(error, "Static update failed in zero-gradient test")
         return
      end if

      if (.not. allocated(model%f1_rA)) then
         call test_failed(error, "Static update did not cache f1_rA")
         return
      end if

      call check(error, size(model%f1_rA, 1), 3, more="Gradient first dimension must be Cartesian")
      if (allocated(error)) return
      call check(error, size(model%f1_rA, 2), mol%nat, more="Gradient second dimension must be nat")
      if (allocated(error)) return
      call check(error, size(model%f1_rA, 3), mol%nat, more="Gradient third dimension must be nat")
      if (allocated(error)) return
      call check(error, maxval(abs(model%f1_rA)), 0.0_wp, thr=0.0_wp, more="Static gradient must be zero")
   end subroutine test_static_zero_gradient

   subroutine test_static_requires_update(error)
      type(error_type), allocatable, intent(out) :: error

      type(radius_type_static) :: model

      call new_cpcm_radii(model)

      if (allocated(model%f0)) then
         call test_failed(error, "f0 should not be allocated before update")
         return
      end if

      if (allocated(model%f1_rA)) then
         call test_failed(error, "f1_rA should not be allocated before update")
      end if
   end subroutine test_static_requires_update


   subroutine test_custom_atoms_dropcess(error)
      type(error_type), allocatable, intent(out) :: error

      type(structure_type) :: mol
      class(radius_type), allocatable :: model
      type(moist_error_type), allocatable :: err
      real(wp), parameter :: radii(4) = [2.10_wp, 1.25_wp, 1.25_wp, 1.25_wp]

      call make_test_molecule(mol)
      call new_radii_custom_atoms(radii, model, err)
      if (allocated(err)) then
         call test_failed(error, "new_radii_custom_atoms failed: "//trim(err%message))
         return
      end if

      call model%update(mol, err)
      if (allocated(err)) then
         call test_failed(error, "custom atom model update failed: "//trim(err%message))
         return
      end if

      call check(error, size(model%f1_rA, 1), 3, more="custom gradient Cartesian extent")
      if (allocated(error)) return
      call check(error, size(model%f1_rA, 2), mol%nat, more="custom gradient radius extent")
      if (allocated(error)) return
      call check(error, size(model%f1_rA, 3), mol%nat, more="custom gradient atom extent")
      if (allocated(error)) return
      call check(error, size(model%f0), mol%nat, more="custom atom radii size mismatch")
      if (allocated(error)) return
      call check(error, maxval(abs(model%f0 - radii)), 0.0_wp, thr=thr, more="custom atom radii mismatch")
      if (allocated(error)) return
      call check(error, maxval(abs(model%f1_rA)), 0.0_wp, thr=0.0_wp, more="custom atom radii derivative must be zero")
   end subroutine test_custom_atoms_dropcess

   subroutine test_custom_elements_dropcess(error)
      type(error_type), allocatable, intent(out) :: error

      type(structure_type) :: mol
      class(radius_type), allocatable :: model
      type(moist_error_type), allocatable :: err
      !> Test molecule is AlH3 (Z = [13, 1, 1, 1]); element table maps
      !> H -> 1.30, Al -> 2.30; so the per-atom expected values track the
      !> [Al, H, H, H] ordering returned by mstore
      integer, parameter :: atomic_numbers(2) = [1, 13]
      real(wp), parameter :: element_radii(2) = [1.30_wp, 2.30_wp]
      real(wp), parameter :: expected(4) = [2.30_wp, 1.30_wp, 1.30_wp, 1.30_wp]

      call make_test_molecule(mol)
      call new_radii_custom_elements(atomic_numbers, element_radii, model, err)
      if (allocated(err)) then
         call test_failed(error, "new_radii_custom_elements failed: "//trim(err%message))
         return
      end if

      call model%update(mol, err)
      if (allocated(err)) then
         call test_failed(error, "custom element model update failed: "//trim(err%message))
         return
      end if

      call check(error, size(model%f0), mol%nat, more="custom element radii size mismatch")
      if (allocated(error)) return
      call check(error, maxval(abs(model%f0 - expected)), 0.0_wp, thr=thr, more="custom element radii mismatch")
      if (allocated(error)) return
      call check(error, maxval(abs(model%f1_rA)), 0.0_wp, thr=0.0_wp, more="custom element radii derivative must be zero")
   end subroutine test_custom_elements_dropcess

   subroutine test_custom_atoms_empty_fails(error)
      type(error_type), allocatable, intent(out) :: error

      class(radius_type), allocatable :: model
      type(moist_error_type), allocatable :: err
      real(wp) :: radii(0)

      call new_radii_custom_atoms(radii, model, err)
      if (allocated(err)) call test_failed(error, trim(err%message))
   end subroutine test_custom_atoms_empty_fails

   subroutine test_custom_atoms_nonpositive_fails(error)
      type(error_type), allocatable, intent(out) :: error

      class(radius_type), allocatable :: model
      type(moist_error_type), allocatable :: err

      call new_radii_custom_atoms([1.20_wp, 0.0_wp], model, err)
      if (allocated(err)) call test_failed(error, trim(err%message))
   end subroutine test_custom_atoms_nonpositive_fails

   subroutine test_custom_atoms_size_mismatch_fails(error)
      type(error_type), allocatable, intent(out) :: error

      type(structure_type) :: mol
      class(radius_type), allocatable :: model
      type(moist_error_type), allocatable :: err

      call make_test_molecule(mol)
      call new_radii_custom_atoms([2.00_wp, 1.10_wp], model, err)
      if (allocated(err)) then
         call test_failed(error, "unexpected constructor failure: "//trim(err%message))
         return
      end if

      call model%update(mol, err)
      if (allocated(err)) call test_failed(error, trim(err%message))
   end subroutine test_custom_atoms_size_mismatch_fails

   subroutine test_custom_elements_empty_fails(error)
      type(error_type), allocatable, intent(out) :: error

      class(radius_type), allocatable :: model
      type(moist_error_type), allocatable :: err
      integer :: atomic_numbers(0)
      real(wp) :: element_radii(0)

      call new_radii_custom_elements(atomic_numbers, element_radii, model, err)
      if (allocated(err)) call test_failed(error, trim(err%message))
   end subroutine test_custom_elements_empty_fails

   subroutine test_custom_elements_size_mismatch_fails(error)
      type(error_type), allocatable, intent(out) :: error

      class(radius_type), allocatable :: model
      type(moist_error_type), allocatable :: err

      call new_radii_custom_elements([1, 8], [1.20_wp], model, err)
      if (allocated(err)) call test_failed(error, trim(err%message))
   end subroutine test_custom_elements_size_mismatch_fails

   subroutine test_custom_elements_invalid_z_fails(error)
      type(error_type), allocatable, intent(out) :: error

      class(radius_type), allocatable :: model
      type(moist_error_type), allocatable :: err

      call new_radii_custom_elements([0], [1.20_wp], model, err)
      if (allocated(err)) call test_failed(error, trim(err%message))
   end subroutine test_custom_elements_invalid_z_fails

   subroutine test_custom_elements_nonpositive_fails(error)
      type(error_type), allocatable, intent(out) :: error

      class(radius_type), allocatable :: model
      type(moist_error_type), allocatable :: err

      call new_radii_custom_elements([1], [0.0_wp], model, err)
      if (allocated(err)) call test_failed(error, trim(err%message))
   end subroutine test_custom_elements_nonpositive_fails

   subroutine test_custom_elements_duplicate_fails(error)
      type(error_type), allocatable, intent(out) :: error

      class(radius_type), allocatable :: model
      type(moist_error_type), allocatable :: err

      call new_radii_custom_elements([1, 1], [1.20_wp, 1.30_wp], model, err)
      if (allocated(err)) call test_failed(error, trim(err%message))
   end subroutine test_custom_elements_duplicate_fails

   subroutine test_custom_elements_missing_for_molecule_fails(error)
      type(error_type), allocatable, intent(out) :: error

      type(structure_type) :: mol
      class(radius_type), allocatable :: model
      type(moist_error_type), allocatable :: err

      call make_test_molecule(mol)
      call new_radii_custom_elements([1], [1.20_wp], model, err)
      if (allocated(err)) then
         call test_failed(error, "unexpected constructor failure: "//trim(err%message))
         return
      end if

      call model%update(mol, err)
      if (allocated(err)) call test_failed(error, trim(err%message))
   end subroutine test_custom_elements_missing_for_molecule_fails

   subroutine test_custom_string_guidance(error)
      type(error_type), allocatable, intent(out) :: error

      class(radius_type), allocatable :: model
      type(moist_error_type), allocatable :: err

      call new_radii("custom", model, err)
      if (allocated(err)) call test_failed(error, trim(err%message))
   end subroutine test_custom_string_guidance

   subroutine test_custom_radii_cavity_integration(error)
      type(error_type), allocatable, intent(out) :: error

      type(structure_type) :: mol
      type(cavity_type_drop) :: cavity
      class(radius_type), allocatable :: model
      type(moist_error_type), allocatable :: err
      real(wp), parameter :: radii(4) = [2.15_wp, 1.35_wp, 1.35_wp, 1.35_wp]
      !> Local run context borrowed by the cavities built here
      type(moist_context_type), target :: ctx

      call new_context(ctx, verbosity=0)

      call make_test_molecule(mol)

      call new_radii_custom_atoms(radii, model, err)
      if (allocated(err)) then
         call test_failed(error, "new_radii_custom_atoms failed in cavity test: "//trim(err%message))
         return
      end if

      block
         type(moist_cavity_drop_lsf_svdw_type) :: svdw_template
         call svdw_template%new(param=moist_cavity_drop_lsf_svdw_param_type(blend_k=2.5_wp, &
            blend_3b=1.0_wp))
         call new_cavity_drop(cavity, ctx, radius_model=model, lsf_model=svdw_template, error=err, &
            param=moist_cavity_drop_parameters_type(num_leb=110))
      end block
      if (allocated(err)) then
         call test_failed(error, "new_cavity_drop failed with custom radii model: "//trim(err%message))
         return
      end if

      call cavity%update(mol, err)
      if (allocated(err)) then
         call test_failed(error, "cavity update failed with custom radii: "//trim(err%message))
         return
      end if

      if (.not. allocated(cavity%radii)) then
         call test_failed(error, "cavity radii not allocated after update")
         return
      end if

      call check(error, size(cavity%radii), mol%nat, more="cavity radii size mismatch")
      if (allocated(error)) return
      call check(error, maxval(abs(cavity%radii - radii)), 0.0_wp, thr=thr, &
                 more="cavity radii do not match custom radii")
   end subroutine test_custom_radii_cavity_integration

   !> Constructor selectors and public model factory contracts
   subroutine test_factory_contracts(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(radius_type_static) :: static_model
      class(radius_type), allocatable :: model
      type(moist_error_type), allocatable :: err
      integer :: i
      integer, parameter :: tags(5) = [rad_type%cpcm, rad_type%smd, rad_type%d3, rad_type%cosmo, rad_type%bondi]
      character(len=8), parameter :: names(5) = [character(len=8) :: " CPCM ", " SMD ", " D3 ", " COSMO ", " BONDI "]

      do i = 1, size(tags)
         select case (i)
         case (1)
            call new_cpcm_radii(static_model)
         case (2)
            call new_smd_radii(static_model)
         case (3)
            call new_d3_radii(static_model)
         case (4)
            call new_cosmo_radii(static_model)
         case (5)
            call new_bondi_radii(static_model)
         case default
            call test_failed(error, "invalid static constructor index")
            return
         end select
         call check(error, static_model%model_tag, tags(i), more="direct constructor model selector")
         if (allocated(error)) return
         call new_radii(tags(i), model, err, verbosity=3)
         if (allocated(err)) then
            call test_failed(error, "integer factory rejected valid tag")
            return
         end if
         call check_selector()
         if (allocated(error)) return
         call new_radii(names(i), model, err, verbosity=3)
         if (allocated(err)) then
            call test_failed(error, "string factory rejected normalized model name")
            return
         end if
         call check_selector()
         if (allocated(error)) return
      end do
      static_model = default_cpcm_radii(verbosity=4)
      call check(error, static_model%model_tag, rad_type%cpcm, more="default CPCM selector")
      if (allocated(error)) return
      call check(error, static_model%verbosity, 4, more="default CPCM verbosity")
      if (allocated(error)) return
      call new_rahm_radii(static_model)
      call check(error, static_model%model_tag, rad_type%rahm, more="Rahm selector")
      if (allocated(error)) return
      call new_gauss_radii(static_model)
      call check(error, static_model%model_tag, rad_type%gauss, more="Gaussian selector")
      if (allocated(error)) return
      call new_radii(-99, model, err)
      call check(error, allocated(err), more="unknown integer tag must fail")
      if (allocated(error)) return
      call new_radii("unknown", model, err)
      call check(error, allocated(err), more="unknown string name must fail")
      if (allocated(error)) return
      call new_radii_custom_elements([1], [1.3_wp], model, err, verbosity=4)
      if (allocated(err)) then
         call test_failed(error, "custom element constructor rejected valid input")
         return
      end if
      call check(error, model%verbosity, 4, more="custom element explicit verbosity")
      if (allocated(error)) return
      call new_radii_custom_elements([1], [1.3_wp], model, err)
      call check(error, model%verbosity, 0, more="custom element default verbosity")
      if (allocated(error)) return
      call new_radii_custom_atoms([1.3_wp], model, err)
      call check(error, model%verbosity, 0, more="custom atom default verbosity")
   contains
      !> Validate static selector and verbosity in a factory result
      subroutine check_selector()
         select type (model)
         type is (radius_type_static)
            call check(error, model%model_tag, tags(i), more="factory model selector")
         class default
            call test_failed(error, "factory must construct static model")
         end select
         if (allocated(error)) return
         call check(error, model%verbosity, 3, more="factory verbosity")
      end subroutine check_selector
   end subroutine test_factory_contracts

   !> Repeated updates rebuild every cache, and a failed lookup clears them
   subroutine test_cache_lifecycle(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(structure_type) :: mol
      type(radius_type_static) :: static_model
      class(radius_type), allocatable :: model
      type(moist_error_type), allocatable :: err
      real(wp), allocatable :: f0_first(:)
      integer :: i

      call make_test_molecule(mol)
      call new_radii_custom_atoms([2.1_wp, 1.2_wp, 1.3_wp, 1.4_wp], model, err)
      if (allocated(err)) then
         call test_failed(error, "custom constructor failed")
         return
      end if
      ! Poison the caches after each pass; the next update must overwrite them
      do i = 1, 2
         call model%update(mol, err)
         if (allocated(err)) then
            call test_failed(error, "custom repeat update failed")
            return
         end if
         call check(error, model%nat, mol%nat, more="cached custom nat")
         if (allocated(error)) return
         call check(error, maxval(abs(model%f0 - [2.1_wp, 1.2_wp, 1.3_wp, 1.4_wp])), 0.0_wp, thr=thr, &
                    more="custom radii not refreshed")
         if (allocated(error)) return
         call check(error, maxval(abs(model%f1_rA)), 0.0_wp, thr=0.0_wp, more="custom gradient not zeroed")
         if (allocated(error)) return
         model%f0 = 8.0_wp
         model%f1_rA = 8.0_wp
      end do
      call new_cpcm_radii(static_model)
      do i = 1, 2
         call static_model%update(mol, err)
         if (allocated(err)) then
            call test_failed(error, "static repeat update failed")
            return
         end if
         call check(error, static_model%nat, mol%nat, more="cached static nat")
         if (allocated(error)) return
         call check(error, all(static_model%atomic_numbers == mol%num(mol%id)), &
                    more="static element cache not refreshed")
         if (allocated(error)) return
         if (i == 1) f0_first = static_model%f0
         call check(error, maxval(abs(static_model%f0 - f0_first)), 0.0_wp, thr=0.0_wp, &
                    more="static radii not refreshed")
         if (allocated(error)) return
         call check(error, maxval(abs(static_model%f1_rA)), 0.0_wp, thr=0.0_wp, &
                    more="static gradient not zeroed")
         if (allocated(error)) return
         static_model%atomic_numbers = 0
         static_model%f0 = 8.0_wp
         static_model%f1_rA = 8.0_wp
      end do
      mol%num = 0
      call static_model%update(mol, err)
      call check(error, allocated(err), more="static invalid element lookup")
      if (allocated(error)) return
      call check(error, .not. allocated(static_model%f0), more="failed lookup clears radius cache")
      if (allocated(error)) return
      call check(error, .not. allocated(static_model%f1_rA), more="failed lookup clears derivative cache")
      if (allocated(error)) return
      call check(error, .not. allocated(static_model%atomic_numbers), more="failed lookup clears element cache")
      if (allocated(error)) return
      call check(error, static_model%nat, 0, more="failed lookup resets atom count")
   end subroutine test_cache_lifecycle

   !> Invalid empty structures, missing storage, and sparse element overrides
   subroutine test_guard_contracts(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      type(structure_type) :: mol
      type(radius_type_static) :: static_model
      type(radius_type_custom) :: custom_model
      class(radius_type), allocatable :: model
      type(moist_error_type), allocatable :: err

      call make_test_molecule(mol)
      mol%nat = 0
      call new_cpcm_radii(static_model)
      call static_model%update(mol, err)
      call check(error, allocated(err), more="static empty structure guard")
      if (allocated(error)) return
      ! Element mode has no per-atom size check, so only the atom-count guard
      ! can reject an empty structure here
      call new_radii_custom_elements([1], [1.2_wp], model, err)
      if (allocated(err)) then
         call test_failed(error, "custom element constructor failed: "//err%message)
         return
      end if
      call model%update(mol, err)
      call check(error, allocated(err), more="custom empty structure guard")
      if (allocated(error)) return
      call make_test_molecule(mol)
      custom_model%has_atom_radii = .true.
      call custom_model%update(mol, err)
      call check(error, allocated(err), more="custom missing atom storage")
      if (allocated(error)) return
      call check(error, index(err%message, "not allocated") > 0, &
                 more="custom missing atom storage diagnostic")
      if (allocated(error)) return
      custom_model%has_atom_radii = .false.
      custom_model%has_element_radii = .true.
      call custom_model%update(mol, err)
      call check(error, allocated(err), more="custom missing element storage")
      if (allocated(error)) return
      call check(error, index(err%message, "not allocated") > 0, &
                 more="custom missing element storage diagnostic")
      if (allocated(error)) return
      ! With both storages present, only the mode guard rejects either flag pattern
      allocate (custom_model%atom_radii(mol%nat), source=1.2_wp)
      allocate (custom_model%element_radii(14), source=1.2_wp)
      custom_model%has_element_radii = .false.
      call custom_model%update(mol, err)
      call check(error, allocated(err), more="custom unset mode guard")
      if (allocated(error)) return
      custom_model%has_atom_radii = .true.
      custom_model%has_element_radii = .true.
      call custom_model%update(mol, err)
      call check(error, allocated(err), more="custom conflicting mode guard")
      if (allocated(error)) return
      call new_radii_custom_elements([0], [1.2_wp], model, err)
      call check(error, allocated(err), more="invalid atomic number guard")
      if (allocated(error)) return
      call check(error, index(err%message, ">= 1") > 0, more="invalid element diagnostic")
      if (allocated(error)) return
      call new_radii_custom_elements([1, 14], [1.2_wp, 2.4_wp], model, err)
      if (allocated(err)) then
         call test_failed(error, "custom element constructor failed: "//err%message)
         return
      end if
      call model%update(mol, err)
      call check(error, allocated(err), more="custom missing interior element")
   end subroutine test_guard_contracts

end module test_radii
