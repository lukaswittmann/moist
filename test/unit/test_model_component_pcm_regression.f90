!> Change detector: CPCM/COSMO energies, charges and cavities frozen from the
!> pre-pcmgen implementation
!>
!> - Not a reference: re-pin when a deliberate numerical change moves them
!> - Fingerprints: sum, absolute sum, sin and cos projections of each array
module test_model_component_pcm_regression
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

   public :: collect_model_component_pcm_regression

   real(wp), parameter :: dielectric = 78.355_wp
   !> Frozen from the native target CPCM/COSMO implementation with A-matrix assembly.
   !> Rows: iSwiG/CPCM, iSwiG/COSMO, DROP/CPCM, DROP/COSMO.
   !> Re-pinned 2026-10-15 for the shared full-precision Born zeta table (110
   !> points: iSwiG 4.901, DROP 4.900490, now both 4.90101060987); energies
   !> moved 4e-7 (iSwiG) and 1.8e-5 (DROP) relative, positions unchanged.
   real(wp), parameter :: energy_ref(8) = [ &
                          -5.1479216499987563e-2_wp, -5.1152799554327810e-2_wp, &
                          -3.3650283393618854e-2_wp, -3.3436915291446425e-2_wp, &
                          -1.4341653734303421e-1_wp, -1.4250716864515180e-1_wp, &
                          -1.1050896797640573e-1_wp, -1.0980825801523386e-1_wp]
   !> Sum, absolute sum, sin projection and cos projection of the surface charges.
   real(wp), parameter :: q_ref(4, 8) = reshape([ &
                                                2.3321643339626883e-5_wp, 1.3277827942042677e+0_wp, &
                                                1.1480374893988068e-2_wp, -2.6256526865012586e-2_wp, &
                                                2.3173766583087403e-5_wp, 1.3193636527788433e+0_wp, &
                                                1.1407580683766851e-2_wp, -2.6090040739433893e-2_wp, &
                                                1.3377133170037581e-4_wp, 1.0834232266813133e+0_wp, &
                                                1.5609160106153848e-2_wp, 9.4320963195475542e-3_wp, &
                                                1.3292312085976330e-4_wp, 1.0765535086755971e+0_wp, &
                                                1.5510186292786601e-2_wp, 9.3722897358207824e-3_wp, &
                                                1.6550689601669459e-3_wp, 1.8783032443252865e+0_wp, &
                                                -7.1917654405905500e-2_wp, -5.1809274083117209e-2_wp, &
                                                1.6445745783256134e-3_wp, 1.8663933892474505e+0_wp, &
                                                -7.1461642393947296e-2_wp, -5.1480764324172776e-2_wp, &
                                                1.7989591517178052e-3_wp, 1.7146602107797559e+0_wp, &
                                                -6.4603290326013587e-2_wp, -7.2450521074466394e-2_wp, &
                                                1.7875523978540295e-3_wp, 1.7037879755963237e+0_wp, &
                                                -6.4193656882820251e-2_wp, -7.1991130287106868e-2_wp], [4, 8])
   !> Cavity fingerprints for the chloroglycine test, iSwiG then DROP.
   real(wp), parameter :: xyz_ref(4, 2) = reshape([ &
                                                  -2.3098537341174193e+1_wp, 5.2269582477566591e+3_wp, &
                                                  -4.6965463272005962e+1_wp, -3.3973194568477830e+1_wp, &
                                                  -1.9060932221531669e+2_wp, 6.0469792014353961e+3_wp, &
                                                  -6.7810904138394832e+1_wp, -4.6230162296330910e+1_wp], [4, 2])
   real(wp), parameter :: area_ref(4, 2) = reshape([ &
                                                   4.3558592783437673e+2_wp, 4.3558592783437673e+2_wp, &
                                                   2.4177011102698205e+1_wp, 5.7976948118787712e+0_wp, &
                                                   4.7726077146771750e+2_wp, 4.7726077146771750e+2_wp, &
                                                   3.1073264444545163e+1_wp, 8.5213883359103360e+0_wp], [4, 2])
   real(wp), parameter :: switch_ref(4, 2) = reshape([ &
                                                     3.7866513733389091e+2_wp, 3.7866513733389091e+2_wp, &
                                                     1.6819628403573632e+1_wp, 8.5408541941912635e+0_wp, &
                                                     3.7866513733397711e+2_wp, 3.7866513733397711e+2_wp, &
                                                     1.7397685745955502e+1_wp, 9.9670596543947578e+0_wp], [4, 2])

contains

   subroutine collect_model_component_pcm_regression(testsuite)
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
   end subroutine collect_model_component_pcm_regression

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

end module test_model_component_pcm_regression
