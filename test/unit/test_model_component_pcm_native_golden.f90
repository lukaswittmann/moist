!> CPCM/COSMO numerical reference from the pre-pcmgen target implementation
module test_model_component_pcm_native_golden
   use mctc_env, only: wp
   use mctc_env_error, only: moist_error_type => error_type
   use mctc_io, only: structure_type, new
   use mctc_io_convert, only: aatoau
   use mstore, only: get_structure
   use testdrive, only: new_unittest, unittest_type, error_type, test_failed, check
   use moist_context, only: moist_context_type, new_context
   use moist_cavity_type, only: cavity_type
   use moist_cavity_drop, only: cavity_type_drop, new_cavity_drop
   use moist_cavity_drop_parameters, only: moist_cavity_drop_parameters_type
   use moist_cavity_drop_lsf_svdw, only: moist_cavity_drop_lsf_svdw_type
   use moist_cavity_drop_lsf_svdw_param, only: moist_cavity_drop_lsf_svdw_param_type
   use moist_cavity_iswig, only: cavity_type_iswig, new_cavity_iswig, moist_cavity_iswig_parameters_type
   use moist_radii, only: default_cpcm_radii
   use moist_model_continuum_component_pcm_conductor, only: model_continuum_component_cpcm, new_component_cpcm, &
      & model_continuum_component_cosmo, new_component_cosmo
   use moist_channels_coupling, only: coupling_type
   use test_helpers, only: stage_point_charge_energy, component_view
   implicit none(type, external)
   private

   public :: collect_model_component_pcm_native_golden

   real(wp), parameter :: dielectric = 78.355_wp
   !> Frozen from the native target CPCM/COSMO implementation with A-matrix assembly.
   !> Rows: iSwiG/CPCM, iSwiG/COSMO, DROP/CPCM, DROP/COSMO.
   real(wp), parameter :: energy_ref(8) = [ &
                          -5.1479237409116717e-2_wp, -5.1152820330877501e-2_wp, &
                          -3.3650890639200508e-2_wp, -3.3437518686634346e-2_wp, &
                          -1.4341658757774461e-1_wp, -1.4250721856133600e-1_wp, &
                          -1.1051084963947831e-1_wp, -1.0981012774714757e-1_wp]
   !> Sum, absolute sum, sin projection and cos projection of the surface charges.
   real(wp), parameter :: q_ref(4, 8) = reshape([ &
                                                2.3323510907060155e-5_wp, 1.3277834083699780e+0_wp, &
                                                1.1480364912609561e-2_wp, -2.6256547470154389e-2_wp, &
                                                2.3175622308890895e-5_wp, 1.3193642630502787e+0_wp, &
                                                1.1407570765677607e-2_wp, -2.6090061213922686e-2_wp, &
                                                1.3380545566955805e-4_wp, 1.0834448979646505e+0_wp, &
                                                1.5608844324624076e-2_wp, 9.4324166436354218e-3_wp, &
                                                1.3295702845694780e-4_wp, 1.0765750425467024e+0_wp, &
                                                1.5509872513549428e-2_wp, 9.3726080288125454e-3_wp, &
                                                1.6551025782662139e-3_wp, 1.8783039402401367e+0_wp, &
                                                -7.1917754561314939e-2_wp, -5.1809307653944044e-2_wp, &
                                                1.6446079832605175e-3_wp, 1.8663940807496766e+0_wp, &
                                                -7.1461741914296106e-2_wp, -5.1480797682135176e-2_wp, &
                                                1.7998905813674214e-3_wp, 1.7146912113731476e+0_wp, &
                                                -6.4606889890538380e-2_wp, -7.2453728219245297e-2_wp, &
                                                1.7884779215400111e-3_wp, 1.7038187796226330e+0_wp, &
                                                -6.4197233623399186e-2_wp, -7.1994317096176069e-2_wp], [4, 8])
   !> Cavity fingerprints for the chloroglycine test, iSwiG then DROP.
   real(wp), parameter :: xyz_ref(4, 2) = reshape([ &
                                                  -2.3098537341174193e+1_wp, 5.2269582477566591e+3_wp, &
                                                  -4.6965463272005962e+1_wp, -3.3973194568477830e+1_wp, &
                                                  -1.9060932221531669e+2_wp, 6.0469792014353961e+3_wp, &
                                                  -6.7810904138394832e+1_wp, -4.6230162296330910e+1_wp], [4, 2])
   real(wp), parameter :: area_ref(4, 2) = reshape([ &
                                                   4.3558591928526522e+2_wp, 4.3558591928526522e+2_wp, &
                                                   2.4177013710857427e+1_wp, 5.7976936831370516e+0_wp, &
                                                   4.7726049232040492e+2_wp, 4.7726049232040492e+2_wp, &
                                                   3.1073454240033904e+1_wp, 8.5212392029092587e+0_wp], [4, 2])
   real(wp), parameter :: switch_ref(4, 2) = reshape([ &
                                                     3.7866512766592115e+2_wp, 3.7866512766592115e+2_wp, &
                                                     1.6819629720478318e+1_wp, 8.5408515280464066e+0_wp, &
                                                     3.7866466285780655e+2_wp, 3.7866466285780655e+2_wp, &
                                                     1.7397781360936126e+1_wp, 9.9669495405682689e+0_wp], [4, 2])

contains

   subroutine collect_model_component_pcm_native_golden(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
         & new_unittest("iswig_cpcm", test_iswig_cpcm), &
         & new_unittest("iswig_cosmo", test_iswig_cosmo), &
         & new_unittest("drop_cpcm", test_drop_cpcm), &
         & new_unittest("drop_cosmo", test_drop_cosmo), &
         & new_unittest("chloroglycine_iswig_cpcm", test_chloroglycine_iswig_cpcm), &
         & new_unittest("chloroglycine_iswig_cosmo", test_chloroglycine_iswig_cosmo), &
         & new_unittest("chloroglycine_drop_cpcm", test_chloroglycine_drop_cpcm), &
         & new_unittest("chloroglycine_drop_cosmo", test_chloroglycine_drop_cosmo)]
   end subroutine collect_model_component_pcm_native_golden

   subroutine test_iswig_cpcm(error)
      type(error_type), allocatable, intent(out) :: error
      call run_case(error, "iswig", "cpcm")
   end subroutine test_iswig_cpcm

   subroutine test_iswig_cosmo(error)
      type(error_type), allocatable, intent(out) :: error
      call run_case(error, "iswig", "cosmo")
   end subroutine test_iswig_cosmo

   subroutine test_drop_cpcm(error)
      type(error_type), allocatable, intent(out) :: error
      call run_case(error, "drop", "cpcm")
   end subroutine test_drop_cpcm

   subroutine test_drop_cosmo(error)
      type(error_type), allocatable, intent(out) :: error
      call run_case(error, "drop", "cosmo")
   end subroutine test_drop_cosmo

   subroutine test_chloroglycine_iswig_cpcm(error)
      type(error_type), allocatable, intent(out) :: error
      call run_case(error, "iswig", "cpcm", chloroglycine=.true.)
   end subroutine test_chloroglycine_iswig_cpcm

   subroutine test_chloroglycine_iswig_cosmo(error)
      type(error_type), allocatable, intent(out) :: error
      call run_case(error, "iswig", "cosmo", chloroglycine=.true.)
   end subroutine test_chloroglycine_iswig_cosmo

   subroutine test_chloroglycine_drop_cpcm(error)
      type(error_type), allocatable, intent(out) :: error
      call run_case(error, "drop", "cpcm", chloroglycine=.true.)
   end subroutine test_chloroglycine_drop_cpcm

   subroutine test_chloroglycine_drop_cosmo(error)
      type(error_type), allocatable, intent(out) :: error
      call run_case(error, "drop", "cosmo", chloroglycine=.true.)
   end subroutine test_chloroglycine_drop_cosmo

   subroutine run_case(error, cavity_kind, model_kind, chloroglycine)
      type(error_type), allocatable, intent(out) :: error
      character(len=*), intent(in) :: cavity_kind, model_kind
      logical, optional, intent(in) :: chloroglycine

      type(moist_context_type), target :: ctx
      type(structure_type) :: mol
      class(cavity_type), allocatable :: cavity
      type(cavity_type_drop), allocatable :: drop
      type(cavity_type_iswig), allocatable :: iswig
      type(moist_cavity_drop_lsf_svdw_type) :: svdw
      type(model_continuum_component_cpcm) :: cpcm
      type(model_continuum_component_cosmo) :: cosmo
      type(coupling_type) :: coupling
      type(moist_error_type), allocatable :: err
      real(wp), allocatable :: qat(:), q(:)
      real(wp) :: energy, q_fingerprint(4)
      integer :: i, icase, ngrid_ref, icavity

      call new_context(ctx, verbosity=0)
      if (present(chloroglycine)) then
         if (chloroglycine) then
            ! Approximate 2-chloroglycine zwitterion, NH3+-CH(Cl)-COO-.
            ! Atom order: alpha C, carboxyl C, N, O, O, Cl, alpha H, N-H x3.
            call new(mol, [6, 6, 7, 8, 8, 17, 1, 1, 1, 1], reshape([ &
                                                                   0.00_wp, 0.00_wp, 0.00_wp, &
                                                                   1.52_wp, 0.00_wp, 0.00_wp, &
                                                                   -1.45_wp, 0.00_wp, 0.00_wp, &
                                                                   2.15_wp, 1.05_wp, 0.00_wp, &
                                                                   2.15_wp, -1.05_wp, 0.00_wp, &
                                                                   -0.10_wp, 0.00_wp, 1.78_wp, &
                                                                   -0.15_wp, -0.93_wp, -0.35_wp, &
                                                                   -1.80_wp, 0.90_wp, 0.35_wp, &
                                                                   -1.80_wp, -0.70_wp, 0.70_wp, &
                                                                   -1.80_wp, -0.20_wp, -0.90_wp], [3, 10])*aatoau)
            qat = [0.15_wp, 0.70_wp, 0.20_wp, -0.85_wp, -0.85_wp, &
                   -0.15_wp, 0.05_wp, 0.25_wp, 0.25_wp, 0.25_wp]
         end if
      end if
      if (.not. allocated(qat)) then
         call get_structure(mol, "But14diol", "1")
         allocate (qat(mol%nat))
         do i = 1, mol%nat
            qat(i) = 0.3_wp*sin(real(i, wp))
         end do
         qat = qat - sum(qat)/real(mol%nat, wp)
      end if

      if (cavity_kind == "drop") then
         allocate (drop)
         call svdw%new(param=moist_cavity_drop_lsf_svdw_param_type( &
                       blend_k=2.5_wp, blend_2b=1.0_wp, blend_3b=1.0_wp))
         call new_cavity_drop(drop, ctx=ctx, radius_model=default_cpcm_radii(), &
                              lsf_model=svdw, error=err, param=moist_cavity_drop_parameters_type( &
                              num_leb=110, tolerance=1.0e-10_wp, proj_maxiter=1000, proj_level=1))
         if (allocated(err)) then
            call test_failed(error, err%message)
            return
         end if
         call move_alloc(drop, cavity)
      else
         allocate (iswig)
         call new_cavity_iswig(iswig, ctx=ctx, radius_model=default_cpcm_radii(), &
                               error=err, param=moist_cavity_iswig_parameters_type(num_leb=110))
         if (allocated(err)) then
            call test_failed(error, err%message)
            return
         end if
         call move_alloc(iswig, cavity)
      end if
      call cavity%update(mol, error=err)
      if (allocated(err)) then
         call test_failed(error, err%message)
         return
      end if

      if (model_kind == "cpcm") then
         call new_component_cpcm(cpcm, epsilon=dielectric, error=err, ctx=ctx)
         if (.not. allocated(err)) call cpcm%update(mol, cavity, err)
         if (allocated(err)) then
            call test_failed(error, err%message)
            return
         end if
         call stage_point_charge_energy(error, cpcm, cavity, qat, mol, coupling)
         if (allocated(error)) return
         energy = 0.0_wp
         call cpcm%get_energy(component_view(coupling), cavity, energy, err)
         if (.not. allocated(err)) q = cpcm%q
      else
         call new_component_cosmo(cosmo, epsilon=dielectric, error=err, ctx=ctx)
         if (.not. allocated(err)) call cosmo%update(mol, cavity, err)
         if (allocated(err)) then
            call test_failed(error, err%message)
            return
         end if
         call stage_point_charge_energy(error, cosmo, cavity, qat, mol, coupling)
         if (allocated(error)) return
         energy = 0.0_wp
         call cosmo%get_energy(component_view(coupling), cavity, energy, err)
         if (.not. allocated(err)) q = cosmo%q
      end if
      if (allocated(err)) then
         call test_failed(error, err%message)
         return
      end if

      icase = 1
      if (model_kind == "cosmo") icase = icase + 1
      if (cavity_kind == "drop") icase = icase + 2
      icavity = 1
      if (cavity_kind == "drop") icavity = 2
      if (present(chloroglycine)) then
         if (chloroglycine) then
            icase = icase + 4
            call check(error, maxval(abs(fingerprint(reshape(cavity%xyz(:, 1:cavity%ngrid), &
                                                             [3*cavity%ngrid])) - xyz_ref(:, icavity))), 0.0_wp, thr=1.0e-8_wp, &
                       message="native golden cavity positions changed")
            if (allocated(error)) return
            call check(error, maxval(abs(fingerprint(cavity%a(1:cavity%ngrid)) &
                                         - area_ref(:, icavity))), 0.0_wp, thr=1.0e-8_wp, &
                       message="native golden cavity areas changed")
            if (allocated(error)) return
            call check(error, maxval(abs(fingerprint(cavity%f(1:cavity%ngrid)) &
                                         - switch_ref(:, icavity))), 0.0_wp, thr=1.0e-8_wp, &
                       message="native golden cavity switching changed")
            if (allocated(error)) return
         end if
      end if
      ngrid_ref = 997
      if (cavity_kind == "drop") ngrid_ref = 1001
      if (icase > 4) ngrid_ref = 669
      if (icase > 6) ngrid_ref = 670
      call check(error, cavity%ngrid, ngrid_ref, message="native golden cavity grid changed")
      if (allocated(error)) return
      call check(error, energy, energy_ref(icase), thr=1.0e-10_wp, &
                 message="native CPCM/COSMO energy changed")
      if (allocated(error)) return

      q_fingerprint = [sum(q), sum(abs(q)), &
                       sum(q*sin(0.37_wp*[(real(i, wp), i=1, size(q))] + 0.1_wp)), &
                       sum(q*cos(1.13_wp*[(real(i, wp), i=1, size(q))]))]
      call check(error, maxval(abs(q_fingerprint - q_ref(:, icase))), 0.0_wp, &
                 thr=1.0e-10_wp, message="native CPCM/COSMO surface charges changed")
   end subroutine run_case

   function fingerprint(values) result(summary)
      real(wp), intent(in) :: values(:)
      real(wp) :: summary(4)
      integer :: i
      summary = [sum(values), sum(abs(values)), &
                 sum(values*sin(0.37_wp*[(real(i, wp), i=1, size(values))] + 0.1_wp)), &
                 sum(values*cos(1.13_wp*[(real(i, wp), i=1, size(values))]))]
   end function fingerprint

end module test_model_component_pcm_native_golden
