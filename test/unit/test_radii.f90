module test_radii
   use mctc_env, only : wp
   use mctc_env_error, only : moist_error_type => error_type
   use mctc_io, only : structure_type
   use mstore, only : get_structure
   use moist_cavity_drop, only : cavity_type_drop, new_cavity_drop
   use moist_cavity_drop_lsf_svdw, only : moist_cavity_drop_lsf_svdw_type
   use moist_data_radii_legacy, only : get_radius_func
   use moist_radii, only : radius_type, static_radius_type, draco_radius_type
   use moist_radii, only : new_cpcm_radii, new_smd_radii, new_d3_radii
   use moist_radii, only : new_cosmo_radii, new_bondi_radii, new_draco_radii
   use moist_radii, only : new_radii, new_radii_custom_atoms, new_radii_custom_elements
   use testdrive, only : new_unittest, unittest_type, error_type, check, test_failed
   implicit none
   private

   public :: collect_radii

   real(wp), parameter :: thr = 10*epsilon(1.0_wp)

contains

   !> Collect all radii tests.
   subroutine collect_radii(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
         new_unittest("StaticRadiiCPCM", test_static_radii_cpcm), &
         new_unittest("StaticRadiiConstructors", test_static_constructors), &
         new_unittest("RadiiConstructorVerbosity", test_constructor_verbosity), &
         new_unittest("StaticRadiiZeroGradient", test_static_zero_gradient), &
         new_unittest("StaticRadiiNeedsUpdate", test_static_requires_update), &
         new_unittest("DracoRadiiDummy", test_draco_dummy), &
         new_unittest("CustomRadiiAtomsWorks", test_custom_atoms_dropcess), &
         new_unittest("CustomRadiiElementsWorks", test_custom_elements_dropcess), &
         new_unittest("CustomRadiiCavityIntegration", test_custom_radii_cavity_integration), &
         new_unittest("CustomRadiiAtomsBadEmpty", test_custom_atoms_empty_fails, should_fail=.true.), &
         new_unittest("CustomRadiiAtomsBadValue", test_custom_atoms_nonpositive_fails, should_fail=.true.), &
         new_unittest("CustomRadiiAtomsBadNat", test_custom_atoms_size_mismatch_fails, should_fail=.true.), &
         new_unittest("CustomRadiiElementsBadEmpty", test_custom_elements_empty_fails, should_fail=.true.), &
         new_unittest("CustomRadiiElementsBadSize", test_custom_elements_size_mismatch_fails, should_fail=.true.), &
         new_unittest("CustomRadiiElementsBadZ", test_custom_elements_invalid_z_fails, should_fail=.true.), &
         new_unittest("CustomRadiiElementsBadValue", test_custom_elements_nonpositive_fails, should_fail=.true.), &
         new_unittest("CustomRadiiElementsBadDuplicate", test_custom_elements_duplicate_fails, should_fail=.true.), &
         new_unittest("CustomRadiiElementsMissing", test_custom_elements_missing_for_molecule_fails, &
            should_fail=.true.), &
         new_unittest("CustomRadiiStringGuidance", test_custom_string_guidance, should_fail=.true.) &
         ]
   end subroutine collect_radii

   !> Fetch a small reference molecule for the radii suite from mstore.
   !> AlH3 from MB16-43 has composition [Al, H, H, H] (Z = [13, 1, 1, 1],
   !> nat = 4) - enough atomic diversity to exercise both per-element and
   !> per-atom code paths without being large enough to slow the suite.
   !> @param[out] mol  structure populated from MB16-43/AlH3
   subroutine make_test_molecule(mol)
      !> Structure to populate with AlH3 from MB16-43.
      type(structure_type), intent(out) :: mol

      call get_structure(mol, "MB16-43", "AlH3")
   end subroutine make_test_molecule

   subroutine test_static_radii_cpcm(error)
      type(error_type), allocatable, intent(out) :: error

      type(structure_type) :: mol
      type(static_radius_type) :: model
      type(moist_error_type), allocatable :: err
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
         call check(error, model%f0(iat), get_radius_func(mol%num(mol%id(iat)), "cpcm"), thr=thr, &
            more="CPCM static radii mismatch")
         if (allocated(error)) return
      end do
   end subroutine test_static_radii_cpcm

   subroutine test_static_constructors(error)
      type(error_type), allocatable, intent(out) :: error

      type(structure_type) :: mol
      type(static_radius_type) :: model
      type(moist_error_type), allocatable :: err

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
      call check(error, model%f0(1), get_radius_func(mol%num(mol%id(1)), "smd"), thr=thr, &
         more="SMD radius mismatch")
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
      call check(error, model%f0(1), get_radius_func(mol%num(mol%id(1)), "d3"), thr=thr, &
         more="D3 radius mismatch")
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
      call check(error, model%f0(1), get_radius_func(mol%num(mol%id(1)), "cosmo"), thr=thr, &
         more="COSMO radius mismatch")
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
      call check(error, model%f0(1), get_radius_func(mol%num(mol%id(1)), "bondi"), thr=thr, &
         more="Bondi radius mismatch")
   end subroutine test_static_constructors

   subroutine test_constructor_verbosity(error)
      type(error_type), allocatable, intent(out) :: error

      type(static_radius_type) :: static_model
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
      type(static_radius_type) :: model
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

      type(static_radius_type) :: model

      call new_cpcm_radii(model)

      if (allocated(model%f0)) then
         call test_failed(error, "f0 should not be allocated before update")
         return
      end if

      if (allocated(model%f1_rA)) then
         call test_failed(error, "f1_rA should not be allocated before update")
      end if
   end subroutine test_static_requires_update

   subroutine test_draco_dummy(error)
      type(error_type), allocatable, intent(out) :: error

      type(structure_type) :: mol
      type(draco_radius_type) :: model
      type(moist_error_type), allocatable :: err

      call make_test_molecule(mol)
      call new_draco_radii(model)

      call model%update(mol, err)
      if (.not. allocated(err)) then
         call test_failed(error, "DRACO dummy update should report not implemented")
      end if
   end subroutine test_draco_dummy

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
      !> [Al, H, H, H] ordering returned by mstore.
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

      call make_test_molecule(mol)

      call new_radii_custom_atoms(radii, model, err)
      if (allocated(err)) then
         call test_failed(error, "new_radii_custom_atoms failed in cavity test: "//trim(err%message))
         return
      end if

      block
         type(moist_cavity_drop_lsf_svdw_type) :: svdw_template
         call svdw_template%new(blend_k=2.5_wp, blend_3b=1.0_wp)
         call new_cavity_drop(cavity, nleb=110, &
            debug=.false., verbose=0, radius_model=model, &
            lsf_model=svdw_template, error=err)
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

end module test_radii
