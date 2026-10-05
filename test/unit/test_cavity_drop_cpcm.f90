module test_cavity_drop_cpcm
   use moist_cavity_drop_lsf_svdw_param, only: moist_cavity_drop_lsf_svdw_param_type
   use moist_cavity_drop_parameters, only: moist_cavity_drop_parameters_type
   use mctc_env_accuracy, only: wp
   use mctc_env_error, only: mctc_error => error_type
   use mctc_io, only: structure_type, new
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use testdrive, only: to_string, test_failed
   use moist_cavity_drop, only: cavity_type_drop, new_cavity_drop
   use moist_cavity_drop_lsf_svdw, only: moist_cavity_drop_lsf_svdw_type
   use moist_radii, only: default_cpcm_radii
   use mstore, only: get_structure
   use moist_math_lapack, only: getrf, getri
   use moist_model_continuum_component_pcm_amat, only: &
      & assemble_pcm_amat_with_gradient, pcm_amat_surface_weights, &
      & pcm_amat_nuclear_gradient
   use moist_model_continuum_component_pcm_electrostatics, only: &
      & pcm_electrostatic_nuclear_gradient
   use moist_context, only: moist_context_type, new_context
   use test_helpers, only: fd6_scalar, fd6_offsets
   use, intrinsic :: iso_fortran_env, only: error_unit
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   implicit none(type, external)
   private

   public :: collect_cavity_drop_cpcm

   integer, parameter :: ndim = 3

   !> SvdW blending steepness
   real(wp), parameter :: k = 1.5_wp
   !> Three-body SvdW blending weight
   real(wp), parameter :: gamma = 1.0_wp
   !> Lebedev rule size
   integer, parameter :: NUM_LEB = 26

   !> Fourth-order point-charge finite-difference step in bohr
   real(wp), parameter :: POINTCHARGE_STEP = 4.0e-3_wp
   !> Point-charge finite-difference absolute tolerance
   real(wp), parameter :: POINTCHARGE_ATOL = 1.0e-10_wp
   !> Point-charge finite-difference relative tolerance
   real(wp), parameter :: POINTCHARGE_RTOL = 1.0e-9_wp

   !> Precision of the independent CPCM reference
   integer, parameter :: rk = wp
   ! Quadruple precision, about 20x slower (!?); swap in with MATRIX_FLOOR and REFINE_TOL below
   ! integer, parameter :: rk = selected_real_kind(30, 300)

   !> Sixth-order matrix finite-difference step in bohr
   real(wp), parameter :: MATRIX_STEP = 1.0e-4_wp
   !> Matrix finite-difference absolute tolerance
   real(wp), parameter :: MATRIX_ATOL = 1.0e-10_wp
   !> Matrix finite-difference relative tolerance
   real(wp), parameter :: MATRIX_RTOL = 1.0e-9_wp
   !> Stencil round-off floor in units of epsilon(wp)*|A_ij|/step
   real(wp), parameter :: MATRIX_FLOOR = 1.0e4_wp
   ! Quadruple precision: only the rounding of the samples to wp for fd6_scalar remains
   ! real(wp), parameter :: MATRIX_FLOOR = 4.0_wp

   !> Absolute tolerance of production against reference matrix values
   real(wp), parameter :: VALUE_ATOL = 1.0e-10_wp
   !> Relative tolerance of production against reference matrix values
   real(wp), parameter :: VALUE_RTOL = 1.0e-10_wp

   !> KKT residual at which a reference branch counts as refined
   real(rk), parameter :: REFINE_TOL = 1.0e-13_rk
   ! real(rk), parameter :: REFINE_TOL = 1.0e-27_rk

   !> Largest accepted entry of A A^-1 - 1
   real(wp), parameter :: INVERSE_THR = 1.0e-10_wp

   !> DROP projection settings
   real(wp), parameter :: PROJ_TOL = 1.0e-14_wp
   integer, parameter :: PROJ_MAXITER = 150
   integer, parameter :: PROJ_LEVEL = 2

contains

   subroutine collect_cavity_drop_cpcm(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("contract_amat1_q1q2", test_contract_amat1_q1q2_rA), &
                  new_unittest("contract_nuc_elec_pointcharge_fd", &
                               test_contract_nuc_elec_pointcharge_fd), &
                  new_unittest("single_atom", test_single_atom), &
                  new_unittest("dimer", test_dimer), &
                  new_unittest("ar5_blendk_09", test_ar5_blendk_09), &
                  new_unittest("bih3_h2o", test_bih3_h2o), &
                  new_unittest("mb16_43_01", test_mb16_43_01), &
                  new_unittest("mb16_43_19", test_mb16_43_19), &
                  new_unittest("but14diol_32", test_but14diol_32), &
                  new_unittest("il16_008", test_il16_008), &
                  new_unittest("mb16_43_h2", test_mb16_43_h2), &
                  new_unittest("heavy28_pbh4", test_heavy28_pbh4), &
                  new_unittest("negative_weight_lebedev_rejected", test_negative_weight_lebedev) &
                  ]
   end subroutine collect_cavity_drop_cpcm

   !> DROP refuses the 74-point Lebedev rule, which has negative weights
   !>
   !> - At construction, where no fitted Born zeta exists for that size
   !> - In `update`, if the size is changed on a constructed cavity: the
   !>   Lebedev cache must not hand negative weights to the surface
   subroutine test_negative_weight_lebedev(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      type(structure_type) :: mol
      type(cavity_type_drop), allocatable :: cavity
      type(moist_cavity_drop_lsf_svdw_type) :: svdw_template
      type(mctc_error), allocatable :: cavity_error
      !> Local run context borrowed by the cavities built here
      type(moist_context_type), target :: ctx

      call new_context(ctx, verbosity=0)
      call new(mol, [1, 1], reshape([0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 1.4_wp], [3, 2]))
      call svdw_template%new(param=moist_cavity_drop_lsf_svdw_param_type(blend_k=k, blend_3b=gamma))

      allocate (cavity)
      call new_cavity_drop(cavity, ctx, radius_model=default_cpcm_radii(), lsf_model=svdw_template, &
         error=cavity_error, param=moist_cavity_drop_parameters_type(num_leb=74))
      call check(error, allocated(cavity_error), "DROP was constructed with the 74-point Lebedev rule")
      if (allocated(error)) return
      deallocate (cavity, cavity_error)

      allocate (cavity)
      call new_cavity_drop(cavity, ctx, radius_model=default_cpcm_radii(), lsf_model=svdw_template, &
         error=cavity_error, param=moist_cavity_drop_parameters_type(num_leb=NUM_LEB))
      if (allocated(cavity_error)) then
         call test_failed(error, cavity_error%message)
         return
      end if
      cavity%param%num_leb = 74
      call cavity%update(mol, error=cavity_error)
      if (.not. allocated(cavity_error)) then
         call test_failed(error, "DROP update accepted the 74-point Lebedev rule")
         return
      end if
      call check(error, index(cavity_error%message, "negative weights") > 0, &
         & "unexpected error message: "//cavity_error%message)
   end subroutine test_negative_weight_lebedev

   !> Test the contracted A-matrix gradient against the explicit tensor
   !> contraction of the dense derivative built by
   !> `assemble_pcm_amat_with_gradient`
   subroutine test_contract_amat1_q1q2_rA(error)
      type(error_type), allocatable, intent(out) :: error

      type(structure_type) :: mol
      type(cavity_type_drop), allocatable :: cavity
      real(wp), allocatable :: Amat0(:, :), Amat1_rA(:, :, :, :)
      real(wp), allocatable :: q1(:), q2(:)
      real(wp), allocatable :: w_xi(:), w_f(:), w_xyz(:, :)
      real(wp), allocatable :: grad_ref(:, :), grad_ctr(:, :)
      integer :: iat, iaxis, igrid, jgrid, ngrid, nsph
      type(mctc_error), allocatable :: cavity_error
      !> Local run context borrowed by the cavities built here
      type(moist_context_type), target :: ctx

      call new_context(ctx, verbosity=0)

      call get_structure(mol, "MB16-43", "04")

      allocate (cavity)
      block
         type(moist_cavity_drop_lsf_svdw_type) :: svdw_template
         call svdw_template%new(param=moist_cavity_drop_lsf_svdw_param_type(blend_k=k, blend_3b=gamma))
         call new_cavity_drop(cavity, ctx, radius_model=default_cpcm_radii(), lsf_model=svdw_template, &
            error=cavity_error, param=moist_cavity_drop_parameters_type(num_leb=NUM_LEB, tolerance=PROJ_TOL, &
            proj_maxiter=PROJ_MAXITER, proj_level=PROJ_LEVEL))
      end block
      if (allocated(cavity_error)) then
         call test_failed(error, cavity_error%message)
         return
      end if

      call cavity%update(mol, error=cavity_error)
      if (allocated(cavity_error)) then
         call test_failed(error, cavity_error%message)
         return
      end if

      call check(error, cavity%nsph, mol%nat, "Updated cavity sphere count")
      if (allocated(error)) return

      call cavity%get_gradient(cavity_error)
      if (allocated(cavity_error)) then
         call test_failed(error, cavity_error%message)
         return
      end if

      ngrid = cavity%ngrid
      nsph = cavity%nsph

      allocate (Amat0(ngrid, ngrid), Amat1_rA(3, nsph, ngrid, ngrid))
      call assemble_pcm_amat_with_gradient(cavity%xi0, cavity%f, cavity%xyz, &
                                           cavity%xi1_rA, cavity%f1_rA, cavity%xyz1_rA, &
                                           Amat0, Amat1_rA, cavity_error)
      if (allocated(cavity_error)) then
         call test_failed(error, cavity_error%message)
         return
      end if

      allocate (q1(ngrid), q2(ngrid))
      allocate (grad_ref(3, nsph), grad_ctr(3, nsph))

      do igrid = 1, ngrid
         q1(igrid) = real(igrid, wp)/(real(ngrid, wp) + 1.0_wp)
         if (mod(igrid, 2) == 0) then
            q2(igrid) = -1.0_wp/(real(igrid, wp) + 0.5_wp)
         else
            q2(igrid) = 1.0_wp/(real(igrid, wp) + 0.25_wp)
         end if
      end do

      grad_ref = 0.0_wp
      do iat = 1, nsph
         do iaxis = 1, 3
            do igrid = 1, ngrid
               do jgrid = 1, ngrid
                  grad_ref(iaxis, iat) = grad_ref(iaxis, iat) &
                                         + q1(igrid)*Amat1_rA(iaxis, iat, igrid, jgrid)*q2(jgrid)
               end do
            end do
         end do
      end do
      ! The dense contraction is the expected value below; check() only screens `actual`
      call check(error, all(ieee_is_finite(grad_ref)), &
                 "Non-finite dense A-matrix gradient contraction")
      if (allocated(error)) return

      allocate (w_xi(ngrid), w_f(ngrid), w_xyz(3, ngrid))
      call pcm_amat_surface_weights(cavity%xi0, cavity%f, cavity%xyz, q1, q2, &
                                    w_xi, w_f, w_xyz, cavity_error)
      if (allocated(cavity_error)) then
         call test_failed(error, cavity_error%message)
         return
      end if
      call pcm_amat_nuclear_gradient(cavity%xi1_rA, cavity%f1_rA, cavity%xyz1_rA, &
                                     w_xi, w_f, w_xyz, grad_ctr, cavity_error)
      if (allocated(cavity_error)) then
         call test_failed(error, cavity_error%message)
         return
      end if

      do iat = 1, nsph
         do iaxis = 1, 3
            call check(error, &
                       grad_ctr(iaxis, iat), &
                       grad_ref(iaxis, iat), &
                       thr_abs=5.0e-11_wp, thr_rel=5.0e-11_wp, &
                       more="contracted A-matrix gradient mismatch")
            if (allocated(error)) return
         end do
      end do
   end subroutine test_contract_amat1_q1q2_rA

   !> Test the fused nuclear/electronic contraction against fourth-order finite differences
   !>
   !> This checks the point-charge case: the host total position weight is
   !> `w_xyz(:, i) = q_i grad phi_nuc(r_i)` of the nuclear potential
   subroutine test_contract_nuc_elec_pointcharge_fd(error)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error

      type(structure_type) :: mol, mol_fd
      type(cavity_type_drop), allocatable :: cavity
      real(wp), allocatable :: surface_q(:), w_xyz(:, :), za(:)
      real(wp), allocatable :: grad_ctr(:, :), grad_num(:, :)
      integer, allocatable :: numbering_ref(:)
      integer :: iat, iaxis, igrid, ngrid, istencil
      real(wp) :: energies(-2:2), r_vec(3), r_dist
      type(mctc_error), allocatable :: cavity_error
      !> Local run context borrowed by the cavities built here
      type(moist_context_type), target :: ctx

      call new_context(ctx, verbosity=0)

      call get_structure(mol, "MB16-43", "15")

      allocate (cavity)
      block
         type(moist_cavity_drop_lsf_svdw_type) :: svdw_template
         call svdw_template%new(param=moist_cavity_drop_lsf_svdw_param_type(blend_k=k, blend_3b=gamma))
         call new_cavity_drop(cavity, ctx, radius_model=default_cpcm_radii(), lsf_model=svdw_template, &
            error=cavity_error, param=moist_cavity_drop_parameters_type(num_leb=NUM_LEB, tolerance=PROJ_TOL, &
            proj_maxiter=PROJ_MAXITER, proj_level=PROJ_LEVEL, wleb_prune_level=3))
      end block
      if (allocated(cavity_error)) then
         call test_failed(error, cavity_error%message)
         return
      end if

      call cavity%update(mol, error=cavity_error)
      if (allocated(cavity_error)) then
         call test_failed(error, cavity_error%message)
         return
      end if

      ngrid = cavity%ngrid
      allocate (numbering_ref(ngrid))
      numbering_ref = cavity%numbering(1:ngrid)

      call cavity%get_gradient(cavity_error)
      if (allocated(cavity_error)) then
         call test_failed(error, cavity_error%message)
         return
      end if

      allocate (surface_q(ngrid), w_xyz(3, ngrid), za(mol%nat))
      allocate (grad_ctr(3, mol%nat), grad_num(3, mol%nat))

      do igrid = 1, ngrid
         if (mod(igrid, 2) == 0) then
            surface_q(igrid) = -0.07_wp/(real(igrid, wp) + 0.5_wp)
         else
            surface_q(igrid) = 0.09_wp/(real(igrid, wp) + 0.25_wp)
         end if
      end do

      do iat = 1, mol%nat
         za(iat) = real(mol%num(mol%id(iat)), wp)
      end do
      ! Host total position weight of the nuclear point-charge potential
      w_xyz = 0.0_wp
      do igrid = 1, ngrid
         do iat = 1, mol%nat
            r_vec = cavity%xyz(:, igrid) - mol%xyz(:, iat)
            r_dist = norm2(r_vec)
            w_xyz(:, igrid) = w_xyz(:, igrid) - surface_q(igrid)*za(iat)*r_vec/(r_dist*r_dist*r_dist)
         end do
      end do

      call pcm_electrostatic_nuclear_gradient(cavity%xyz, cavity%sphxyz, &
                                              cavity%xyz1_rA, surface_q, w_xyz, za, &
                                              grad_ctr, cavity_error)
      if (allocated(cavity_error)) then
         call test_failed(error, cavity_error%message)
         return
      end if

      grad_num = 0.0_wp
      do iat = 1, mol%nat
         do iaxis = 1, 3
            do istencil = -2, 2
               if (istencil == 0) cycle
               mol_fd = mol
               mol_fd%xyz(iaxis, iat) = mol_fd%xyz(iaxis, iat) + real(istencil, wp)*POINTCHARGE_STEP
               call cavity%update(mol_fd, error=cavity_error)
               if (allocated(cavity_error)) then
                  call test_failed(error, cavity_error%message)
                  return
               end if
               if (cavity%ngrid /= ngrid) then
                  call test_failed(error, "contract_nuc_elec FD: ngrid changed for stencil point")
                  return
               end if
               if (any(cavity%numbering(1:ngrid) /= numbering_ref)) then
                  call test_failed(error, "contract_nuc_elec FD: numbering changed for stencil point")
                  return
               end if
               energies(istencil) = weighted_nuclear_potential(cavity, mol_fd, surface_q, za)
            end do
            ! Pair opposite energy values before combining the fourth-order stencil
            grad_num(iaxis, iat) = (8.0_wp*(energies(1) - energies(-1)) &
               & - (energies(2) - energies(-2)))/(12.0_wp*POINTCHARGE_STEP)
         end do
      end do

      if (.not. all(ieee_is_finite(grad_ctr)) .or. .not. all(ieee_is_finite(grad_num))) then
         call test_failed(error, "contract_nuc_elec FD: nonfinite gradient")
         return
      end if

      ! Restore reference geometry
      call cavity%update(mol, error=cavity_error)
      if (allocated(cavity_error)) then
         call test_failed(error, cavity_error%message)
         return
      end if

      do iat = 1, mol%nat
         do iaxis = 1, 3
            call check(error, &
                       grad_ctr(iaxis, iat), &
                       grad_num(iaxis, iat), &
                       thr_abs=POINTCHARGE_ATOL, thr_rel=POINTCHARGE_RTOL, &
                       more="pcm_electrostatic_nuclear_gradient point-charge FD mismatch")
            if (allocated(error)) return
         end do
      end do
   end subroutine test_contract_nuc_elec_pointcharge_fd

   pure function weighted_nuclear_potential(cavity, mol, surface_q, za) result(value)
      type(cavity_type_drop), intent(in) :: cavity
      type(structure_type), intent(in) :: mol
      real(wp), intent(in) :: surface_q(:)
      real(wp), intent(in) :: za(:)
      real(wp) :: value

      integer :: igrid, katom
      real(wp) :: r_vec(3), r

      value = 0.0_wp
      do igrid = 1, cavity%ngrid
         do katom = 1, cavity%nsph
            r_vec(:) = cavity%xyz(:, igrid) - mol%xyz(:, katom)
            r = sqrt(dot_product(r_vec, r_vec))
            if (r > 1.0e-12_wp) then
               value = value + surface_q(igrid)*za(katom)/r
            end if
         end do
      end do
   end function weighted_nuclear_potential

   !> Test A matrix gradient for a single atom
   subroutine test_single_atom(error)
      type(error_type), allocatable, intent(out) :: error
      type(structure_type) :: mol

      ! Create single oxygen atom
      call new(mol, [8], reshape([0.0_wp, 0.0_wp, 0.0_wp], [3, 1]))

      call do_test(error, mol)
   end subroutine test_single_atom

   !> Test A matrix gradient for a dimer
   subroutine test_dimer(error)
      type(error_type), allocatable, intent(out) :: error
      type(structure_type) :: mol

      ! Create dimer (two oxygen atoms)
      call new(mol, [8, 8], reshape([0.0_wp, 0.0_wp, 0.0_wp, &
                                     3.0_wp, 0.0_wp, 0.0_wp], [3, 2]))

      call do_test(error, mol)
   end subroutine test_dimer

   !> Test A matrix gradient for 5-argon geometry with custom blend-k
   subroutine test_ar5_blendk_09(error)
      type(error_type), allocatable, intent(out) :: error
      type(structure_type) :: mol

      call new(mol, [18, 18, 18, 18, 18], reshape([ &
                                                  0.2_wp, 0.0_wp, 5.1_wp, &
                                                  -2.2_wp, -2.2_wp, 0.0_wp, &
                                                  2.2_wp, -2.2_wp, 0.0_wp, &
                                                  -2.2_wp, 2.2_wp, 0.0_wp, &
                                                  2.2_wp, 2.2_wp, 0.0_wp], [3, 5]))

      call do_test(error, mol, blend_k_override=0.8_wp)
   end subroutine test_ar5_blendk_09

   !> Test A matrix gradient for bih3_h2o system
   subroutine test_bih3_h2o(error)
      type(error_type), allocatable, intent(out) :: error
      type(structure_type) :: mol

      call get_structure(mol, "Heavy28", "bih3_h2o")

      call do_test(error, mol)
   end subroutine test_bih3_h2o

   !> Test A matrix gradient for MB16-43 01
   subroutine test_mb16_43_01(error)
      type(error_type), allocatable, intent(out) :: error
      type(structure_type) :: mol

      call get_structure(mol, "MB16-43", "01")

      call do_test(error, mol)
   end subroutine test_mb16_43_01

   !> Test A matrix gradient for MB16-43 19
   subroutine test_mb16_43_19(error)
      type(error_type), allocatable, intent(out) :: error
      type(structure_type) :: mol

      call get_structure(mol, "MB16-43", "19")

      call do_test(error, mol)
   end subroutine test_mb16_43_19

   !> Test A matrix gradient for But14diol 32
   subroutine test_but14diol_32(error)
      type(error_type), allocatable, intent(out) :: error
      type(structure_type) :: mol

      call get_structure(mol, "But14diol", "32")

      call do_test(error, mol)
   end subroutine test_but14diol_32

   !> Test A matrix gradient for IL16 008
   subroutine test_il16_008(error)
      type(error_type), allocatable, intent(out) :: error
      type(structure_type) :: mol

      call get_structure(mol, "IL16", "008")

      call do_test(error, mol)
   end subroutine test_il16_008

   !> Test A matrix gradient for MB16-43 H2
   subroutine test_mb16_43_h2(error)
      type(error_type), allocatable, intent(out) :: error
      type(structure_type) :: mol

      call get_structure(mol, "MB16-43", "H2")

      call do_test(error, mol)
   end subroutine test_mb16_43_h2

   !> Test A matrix gradient for Heavy28 PbH4, displaced off Td
   subroutine test_heavy28_pbh4(error)
      type(error_type), allocatable, intent(out) :: error
      type(structure_type) :: mol

      call get_structure(mol, "Heavy28", "pbh4")
      mol%xyz = mol%xyz + 1.0e-2_wp*reshape([ &
                                            1.0_wp, -0.7_wp, 0.5_wp, &
                                            0.4_wp, -0.3_wp, -0.2_wp, &
                                            -0.8_wp, 1.0_wp, -0.8_wp, &
                                            -0.8_wp, -0.7_wp, -1.0_wp, &
                                            0.3_wp, -0.3_wp, -0.8_wp], [3, 5])

      call do_test(error, mol)
   end subroutine test_heavy28_pbh4

   !> Test the A matrix and its nuclear gradient against the independent reference
   !>
   !> - Production matrix values match the reference at the reference geometry
   !> - A A^-1 recovers the identity
   !> - Analytic derivatives match sixth-order finite differences of the
   !>   reference, on points that stay converged across the whole stencil
   !> - Sixth order because xi ~ f_foc^(-1/2) is stiff deep in the f_foc tail:
   !>   at f_foc = 4e-9 (MB16-43 19, A_302,302) the fourth-order stencil
   !>   truncates at 1.2e-9 relative
   !>
   !> @param[out] error             Test failure
   !> @param[in]  mol               Molecular structure
   !> @param[in]  blend_k_override  SvdW blending steepness override (optional)
   subroutine do_test(error, mol, blend_k_override)
      !> Test failure
      type(error_type), allocatable, intent(out) :: error
      !> Molecular structure
      type(structure_type), intent(in) :: mol
      !> SvdW blending steepness override
      real(wp), intent(in), optional :: blend_k_override

      type(structure_type) :: mol_fd
      type(cavity_type_drop), allocatable :: cavity
      type(mctc_error), allocatable :: cavity_error
      !> Local run context borrowed by the cavity built here
      type(moist_context_type), target :: ctx
      !> Production matrix and its analytic nuclear derivative
      real(wp), allocatable :: Amat0(:, :), Amat1_rA(:, :, :, :)
      !> Production inverse and the A A^-1 - 1 residual
      real(wp), allocatable :: Amat_inv(:, :), residual(:, :)
      integer, allocatable :: ipiv(:)
      !> Reference matrix at the current geometry
      real(rk), allocatable :: ref_values(:, :)
      !> Reference matrix at the six stencil geometries, on reference grid indices
      real(rk), allocatable :: samples(:, :, :)
      !> Anchor offsets from their owner nuclei, at the reference and current grid
      real(rk), allocatable :: anchor_offset(:, :), displaced_anchors(:, :)
      !> Displaced nuclear positions
      real(rk), allocatable :: centers(:, :)
      !> Persistent grid numbering -> reference grid index, and its current image
      integer, allocatable :: numbering_to_idx(:), idx(:)
      !> Points converged at the reference geometry, and at every stencil geometry
      logical, allocatable :: valid_ref(:), valid(:), seen(:)
      integer :: iat, idir, istep, igrid, jgrid, ngrid, info
      integer :: worst_loc(2)
      real(wp) :: blend_k_local, analytic, numeric, tol
      !> Worst |deviation|/tolerance seen, and where it occurred
      real(wp) :: worst_ratio, worst_a, worst_n, worst_tol
      integer :: worst_iat, worst_idir, worst_i, worst_j

      blend_k_local = k
      if (present(blend_k_override)) blend_k_local = blend_k_override
      worst_ratio = 0.0_wp; worst_a = 0.0_wp; worst_n = 0.0_wp; worst_tol = 0.0_wp
      worst_iat = 0; worst_idir = 0; worst_i = 0; worst_j = 0

      call new_context(ctx, verbosity=0, debug=.false.)
      allocate (cavity)
      block
         type(moist_cavity_drop_lsf_svdw_type) :: svdw_template
         call svdw_template%new(param=moist_cavity_drop_lsf_svdw_param_type(blend_k=blend_k_local, &
            blend_3b=gamma))
         call new_cavity_drop(cavity, ctx, radius_model=default_cpcm_radii(), lsf_model=svdw_template, &
            error=cavity_error, param=moist_cavity_drop_parameters_type(num_leb=NUM_LEB, tolerance=PROJ_TOL, &
            proj_maxiter=PROJ_MAXITER, proj_level=PROJ_LEVEL))
      end block
      if (.not. allocated(cavity_error)) call cavity%update(mol, error=cavity_error)
      if (.not. allocated(cavity_error)) call cavity%get_gradient(cavity_error)
      if (allocated(cavity_error)) then
         call test_failed(error, cavity_error%message)
         return
      end if
      ngrid = cavity%ngrid

      ! Reference grid indices follow the persistent numbering
      allocate (numbering_to_idx(max(1, maxval([0, cavity%numbering(1:ngrid)]))), source=0)
      do igrid = 1, ngrid
         if (cavity%numbering(igrid) > 0) numbering_to_idx(cavity%numbering(igrid)) = igrid
      end do
      valid_ref = cavity%converged(1:ngrid)
      allocate (anchor_offset(3, ngrid))
      do igrid = 1, ngrid
         anchor_offset(:, igrid) = real(cavity%anchorxyz(:, igrid), rk) &
                                   - real(mol%xyz(:, cavity%owner(igrid)), rk)
      end do

      allocate (Amat0(ngrid, ngrid), Amat1_rA(ndim, mol%nat, ngrid, ngrid))
      call assemble_pcm_amat_with_gradient(cavity%xi0, cavity%f, cavity%xyz, &
                                           cavity%xi1_rA, cavity%f1_rA, cavity%xyz1_rA, &
                                           Amat0, Amat1_rA, cavity_error)
      if (allocated(cavity_error)) then
         call test_failed(error, "assemble_pcm_amat_with_gradient failed: "//cavity_error%message)
         return
      end if

      !> Production values against the reference at the reference geometry
      allocate (ref_values(ngrid, ngrid))
      call cpcm_reference_matrix(cavity, real(mol%xyz, rk), anchor_offset, real(blend_k_local, rk), &
                                 ref_values, error)
      if (allocated(error)) return
      do jgrid = 1, ngrid
         if (.not. valid_ref(jgrid)) cycle
         do igrid = 1, ngrid
            if (.not. valid_ref(igrid)) cycle
            numeric = real(ref_values(igrid, jgrid), wp)
            analytic = Amat0(igrid, jgrid)
            call check(error, ieee_is_finite(analytic) .and. ieee_is_finite(numeric), &
                       "Non-finite CPCM reference matrix value")
            if (allocated(error)) return
            tol = max(VALUE_ATOL, VALUE_RTOL*abs(numeric))
            if (abs(analytic - numeric) > tol) then
               call test_failed(error, "CPCM reference value mismatch at ("//to_string(igrid)//", "// &
                                to_string(jgrid)//")", "production "//to_string(analytic)// &
                                ", reference "//to_string(numeric))
               return
            end if
         end do
      end do

      !> A A^-1 must recover the identity
      allocate (Amat_inv, source=Amat0)
      allocate (ipiv(ngrid))
      call getrf(Amat_inv, ipiv, info)
      if (info == 0) call getri(Amat_inv, ipiv, info)
      if (info /= 0) then
         call test_failed(error, "LAPACK inversion failed with info = "//to_string(info))
         return
      end if
      residual = matmul(Amat0, Amat_inv)
      do igrid = 1, ngrid
         residual(igrid, igrid) = residual(igrid, igrid) - 1.0_wp
      end do
      if (ngrid > 0) then
         worst_loc = maxloc(abs(residual))
         if (abs(residual(worst_loc(1), worst_loc(2))) >= INVERSE_THR) then
            call test_failed(error, "A A^-1 differs from identity by "// &
                             to_string(abs(residual(worst_loc(1), worst_loc(2))))//" at ("// &
                             to_string(worst_loc(1))//", "//to_string(worst_loc(2))//")")
            return
         end if
      end if

      !> Analytic derivatives against finite differences of the reference
      allocate (samples(ngrid, ngrid, size(fd6_offsets)), valid(ngrid), seen(ngrid))
      allocate (centers(3, mol%nat))
      do iat = 1, mol%nat
         do idir = 1, ndim
            valid = valid_ref
            samples = 0.0_rk
            do istep = 1, size(fd6_offsets)
               mol_fd = mol
               mol_fd%xyz(idir, iat) = mol_fd%xyz(idir, iat) + fd6_offsets(istep)*MATRIX_STEP
               call cavity%update(mol_fd, error=cavity_error)
               if (allocated(cavity_error)) then
                  call test_failed(error, cavity_error%message)
                  return
               end if
               centers = real(mol%xyz, rk)
               centers(idir, iat) = centers(idir, iat) + real(fd6_offsets(istep), rk)*real(MATRIX_STEP, rk)

               ! Surviving points keep their rigid anchor offset; new points use the production one
               call reference_indices(cavity%numbering(1:cavity%ngrid), idx)
               allocate (displaced_anchors(3, cavity%ngrid))
               do igrid = 1, cavity%ngrid
                  if (idx(igrid) > 0) then
                     displaced_anchors(:, igrid) = anchor_offset(:, idx(igrid))
                  else
                     displaced_anchors(:, igrid) = real(cavity%anchorxyz(:, igrid), rk) &
                                                   - real(mol_fd%xyz(:, cavity%owner(igrid)), rk)
                  end if
               end do
               deallocate (ref_values)
               allocate (ref_values(cavity%ngrid, cavity%ngrid))
               call cpcm_reference_matrix(cavity, centers, displaced_anchors, real(blend_k_local, rk), &
                                          ref_values, error)
               if (allocated(error)) return
               deallocate (displaced_anchors)

               seen = .false.
               do jgrid = 1, cavity%ngrid
                  if (idx(jgrid) == 0) cycle
                  seen(idx(jgrid)) = cavity%converged(jgrid)
                  do igrid = 1, cavity%ngrid
                     if (idx(igrid) == 0) cycle
                     samples(idx(igrid), idx(jgrid), istep) = ref_values(igrid, jgrid)
                  end do
               end do
               valid = valid .and. seen
            end do
            call check(error, any(valid), "No grid point stays converged across the stencil")
            if (allocated(error)) return

            do jgrid = 1, ngrid
               if (.not. valid(jgrid)) cycle
               do igrid = 1, ngrid
                  if (.not. valid(igrid)) cycle
                  call fd6_scalar(real(samples(igrid, jgrid, 1), wp), real(samples(igrid, jgrid, 2), wp), &
                                  real(samples(igrid, jgrid, 3), wp), real(samples(igrid, jgrid, 4), wp), &
                                  real(samples(igrid, jgrid, 5), wp), real(samples(igrid, jgrid, 6), wp), &
                                  MATRIX_STEP, numeric, error)
                  if (allocated(error)) return
                  analytic = Amat1_rA(idir, iat, igrid, jgrid)
                  call check(error, ieee_is_finite(analytic), "Non-finite CPCM matrix derivative")
                  if (allocated(error)) return
                  ! Stencil round-off scales with the differenced entry, not with its derivative
                  tol = max(MATRIX_ATOL, MATRIX_RTOL*abs(numeric), MATRIX_FLOOR*epsilon(1.0_wp) &
                            *real(maxval(abs(samples(igrid, jgrid, :))), wp)/MATRIX_STEP)
                  if (abs(analytic - numeric)/tol <= worst_ratio) cycle
                  worst_ratio = abs(analytic - numeric)/tol
                  worst_a = analytic; worst_n = numeric; worst_tol = tol
                  worst_iat = iat; worst_idir = idir; worst_i = igrid; worst_j = jgrid
               end do
            end do
         end do
      end do

      ! One assertion at the end rather than a check per entry, so the message
      ! can name the entry that actually decided the outcome
      if (worst_ratio > 1.0_wp) then
         write (error_unit, "(a)") "A-matrix gradient exceeds its tolerance:"
         write (error_unit, "(2x,a,i0,a,i0)") "atom/axis : ", worst_iat, " / ", worst_idir
         write (error_unit, "(2x,a,i0,a,i0)") "entry     : ", worst_i, " , ", worst_j
         write (error_unit, "(2x,a,es15.6)") "analytic  : ", worst_a
         write (error_unit, "(2x,a,es15.6)") "numeric   : ", worst_n
         write (error_unit, "(2x,a,es15.6)") "deviation : ", abs(worst_a - worst_n)
         write (error_unit, "(2x,a,es15.6)") "tolerance : ", worst_tol
         call test_failed(error, "A-matrix gradient off by "// &
                          to_string(worst_ratio)//" times its tolerance")
      end if

   contains

      !> Reference grid index of each current point, zero for points absent at the reference
      !>
      !> @param[in]  numbering  Persistent numbering of the current grid
      !> @param[out] idx        Reference grid indices
      subroutine reference_indices(numbering, idx)
         !> Persistent numbering of the current grid
         integer, intent(in) :: numbering(:)
         !> Reference grid indices
         integer, allocatable, intent(out) :: idx(:)
         integer :: i

         allocate (idx(size(numbering)), source=0)
         do i = 1, size(numbering)
            if (numbering(i) > 0 .and. numbering(i) <= size(numbering_to_idx)) then
               idx(i) = numbering_to_idx(numbering(i))
            end if
         end do
      end subroutine reference_indices

   end subroutine do_test

   !> Independent SvdW value, spatial gradient, and spatial Hessian
   !>
   !> @param[in] point Projected position in bohr, shape (3)
   !> @param[in] centers Nuclear positions in bohr, shape (3, nat)
   !> @param[in] radii Sphere radii in bohr, shape (nat)
   !> @param[in] blend SvdW blending sharpness
   !> @param[in] gamma_r Three-body SvdW blending weight
   !> @param[out] s Level-set value
   !> @param[out] g Spatial level-set gradient, shape (3)
   !> @param[out] hess Spatial level-set Hessian, shape (3, 3)
   subroutine ref_svdw(point, centers, radii, blend, gamma_r, s, g, hess)
      !> Projected position in bohr, shape (3)
      real(rk), intent(in) :: point(3)
      !> Nuclear positions in bohr, shape (3, nat)
      real(rk), intent(in) :: centers(:, :)
      !> Sphere radii in bohr, shape (nat)
      real(rk), intent(in) :: radii(:)
      !> SvdW blending sharpness
      real(rk), intent(in) :: blend
      !> Three-body SvdW blending weight
      real(rk), intent(in) :: gamma_r
      !> Level-set value
      real(rk), intent(out) :: s
      !> Spatial level-set gradient, shape (3)
      real(rk), intent(out) :: g(3)
      !> Spatial level-set Hessian, shape (3, 3)
      real(rk), intent(out) :: hess(3, 3)
      real(rk) :: ps(3), pg(3, 3), ph(3, 3, 3), r, dr(3), n(3), c, v
      real(rk) :: z, zg(3), zh(3, 3), unit(3, 3), outer(3, 3), product_hessian
      integer :: i, j, a, m
      unit = 0.0_rk
      do i = 1, 3
         unit(i, i) = 1.0_rk
      end do
      ps = 0.0_rk; pg = 0.0_rk; ph = 0.0_rk
      do a = 1, size(radii)
         dr = point - centers(:, a)
         r = sqrt(sum(dr*dr)); n = dr/r
         do j = 1, 3
            do i = 1, 3
               outer(i, j) = n(i)*n(j)
            end do
         end do
         do m = 1, 3
            c = blend*real(m, rk)/3.0_rk
            v = exp(-c*(r - radii(a)))
            ps(m) = ps(m) + v
            pg(:, m) = pg(:, m) - c*v*n
            ph(:, :, m) = ph(:, :, m) + v*(c*c*outer - c/r*(unit - outer))
         end do
      end do
      z = ps(3) + gamma_r/6.0_rk*(ps(1)**3 - 3.0_rk*ps(1)*ps(2) + 2.0_rk*ps(3))
      zg = pg(:, 3) + gamma_r/6.0_rk*(3.0_rk*ps(1)**2*pg(:, 1) &
                                      - 3.0_rk*(pg(:, 1)*ps(2) + ps(1)*pg(:, 2)) + 2.0_rk*pg(:, 3))
      do j = 1, 3
         do i = 1, 3
            product_hessian = ph(i, j, 1)*ps(2) + pg(i, 1)*pg(j, 2) &
                              + pg(j, 1)*pg(i, 2) + ps(1)*ph(i, j, 2)
            zh(i, j) = ph(i, j, 3) + gamma_r/6.0_rk*(6.0_rk*ps(1)*pg(i, 1)*pg(j, 1) &
                                                     + 3.0_rk*ps(1)**2*ph(i, j, 1) - 3.0_rk*product_hessian + 2.0_rk*ph(i, j, 3))
            hess(i, j) = -(zh(i, j)/z - zg(i)*zg(j)/(z*z))/blend
         end do
      end do
      s = -log(z)/blend
      g = -zg/(blend*z)
   end subroutine ref_svdw

   !> Independent critical or focal switching value
   !>
   !> @param[in] x Switching input
   !> @param[in] sw Critical or focal switching parameters
   function ref_bump(x, sw) result(value)
      use moist_cavity_drop_switching, only: moist_cavity_drop_swif_sigmoid_bump_type
      !> Switching input
      real(rk), intent(in) :: x
      !> Critical or focal switching parameters
      type(moist_cavity_drop_swif_sigmoid_bump_type), intent(in) :: sw
      !> Switching value
      real(rk) :: value
      real(rk) :: lo, hi, width, u, t, exponent
      lo = real(min(sw%from, sw%to), rk); hi = real(max(sw%from, sw%to), rk)
      width = hi - lo
      if (x <= lo) then
         value = 0.0_rk
      else if (x >= hi) then
         value = 1.0_rk
      else
         u = (hi - x)/width; t = 1.0_rk - u
         exponent = -real(sw%a_hi, rk)*u**(-real(sw%p_hi, rk)) &
                    + real(sw%a_lo, rk)*t**(-real(sw%p_lo, rk))
         if (exponent >= 50.0_rk) then
            value = 0.0_rk
         else if (exponent <= -50.0_rk) then
            value = 1.0_rk
         else
            value = 1.0_rk/(1.0_rk + exp(exponent))
         end if
      end if
      if (sw%from > sw%to) value = 1.0_rk - value
   end function ref_bump

   !> Independent closest-point area and focus weight
   !>
   !> @param[in] cavity Production branch seeds, identities, and fixed cavity parameters
   !> @param[in] igrid Production grid index
   !> @param[in] anchor Unprojected anchor in bohr, shape (3)
   !> @param[in] centers Nuclear positions in bohr, shape (3, nat)
   !> @param[in] blend SvdW blending sharpness
   !> @param[out] refined_point Refined branch position in bohr, shape (3)
   !> @param[out] error Reference failure
   function ref_weight(cavity, igrid, anchor, centers, blend, refined_point, error) result(weight)
      !> Production branch seeds, identities, and fixed cavity parameters
      type(cavity_type_drop), intent(in) :: cavity
      !> Production grid index
      integer, intent(in) :: igrid
      !> Unprojected anchor in bohr, shape (3)
      real(rk), intent(in) :: anchor(3)
      !> Nuclear positions in bohr, shape (3, nat)
      real(rk), intent(in) :: centers(:, :)
      !> SvdW blending sharpness
      real(rk), intent(in) :: blend
      !> Refined branch position in bohr, shape (3)
      real(rk), intent(out), optional :: refined_point(3)
      !> Reference failure
      type(error_type), allocatable, intent(out) :: error
      !> Area, focus, and branch weight before normalization
      real(rk) :: weight
      real(rk) :: point(3), g(3), hess(3, 3), s, gnorm, n(3), sphere(3)
      real(rk) :: amat(3, 3), adj(3, 3), alpha, lambda, detb, trb, beta2, area, focus
      real(rk) :: proj(3, 3), dev(3, 3)
      integer :: i, j
      point = real(cavity%xyz(:, igrid), rk)
      lambda = real(cavity%lambda0(igrid), rk)
      alpha = real(cavity%param%phi_alpha, rk)
      if (cavity%converged(igrid)) then
         call ref_refine(point, lambda, anchor, centers, real(cavity%iswig%radii, rk), blend, alpha, error)
      end if
      if (allocated(error)) then
         weight = 0.0_rk
         return
      end if
      call ref_svdw(point, centers, real(cavity%iswig%radii, rk), blend, real(gamma, rk), s, g, hess)
      if (present(refined_point)) refined_point = point
      gnorm = sqrt(sum(g*g)); n = g/gnorm
      amat = -lambda*hess
      do i = 1, 3
         amat(i, i) = amat(i, i) + alpha
      end do
      adj(1, 1) = amat(2, 2)*amat(3, 3) - amat(2, 3)**2
      adj(2, 2) = amat(1, 1)*amat(3, 3) - amat(1, 3)**2
      adj(3, 3) = amat(1, 1)*amat(2, 2) - amat(1, 2)**2
      adj(1, 2) = amat(1, 3)*amat(2, 3) - amat(1, 2)*amat(3, 3)
      adj(1, 3) = amat(1, 2)*amat(2, 3) - amat(1, 3)*amat(2, 2)
      adj(2, 3) = amat(1, 2)*amat(1, 3) - amat(1, 1)*amat(2, 3)
      adj(2, 1) = adj(1, 2); adj(3, 1) = adj(1, 3); adj(3, 2) = adj(2, 3)
      ! Tangent determinant n^T adj(A) n, independent of tangent-frame choices
      detb = dot_product(n, matmul(adj, n))
      trb = amat(1, 1) + amat(2, 2) + amat(3, 3) - dot_product(n, matmul(amat, n))
      ! Eigenvalue gap as a sum of squares, |P A P - trb/2 P|_F^2/2: trb^2/4 - detb
      ! cancels at umbilics and its square root puts sqrt(eps) noise on beta2
      do j = 1, 3
         do i = 1, 3
            proj(i, j) = -n(i)*n(j)
         end do
         proj(j, j) = proj(j, j) + 1.0_rk
      end do
      dev = matmul(proj, matmul(amat, proj)) - 0.5_rk*trb*proj
      beta2 = 0.5_rk*trb - sqrt(0.5_rk*sum(dev*dev))
      sphere = anchor - centers(:, cavity%owner(igrid)); sphere = sphere/sqrt(sum(sphere*sphere))
      ! Area of the closest-point map between sphere and surface tangent planes
      area = alpha*alpha*abs(dot_product(n, sphere))/abs(detb)
      focus = ref_bump(gnorm, cavity%f_crit)*ref_bump(beta2, cavity%f_foc)
      weight = real(cavity%anchor_wleb0(igrid), rk)*area*focus
   end function ref_weight

   !> Refine one converged production branch with a spatial KKT solve
   !>
   !> @param[in,out] point Projected position in bohr, shape (3)
   !> @param[in,out] lambda Closest-point Lagrange multiplier
   !> @param[in] anchor Unprojected anchor in bohr, shape (3)
   !> @param[in] centers Nuclear positions in bohr, shape (3, nat)
   !> @param[in] radii Sphere radii in bohr, shape (nat)
   !> @param[in] blend SvdW blending sharpness
   !> @param[in] alpha Closest-point objective coefficient
   !> @param[out] error Reference failure
   subroutine ref_refine(point, lambda, anchor, centers, radii, blend, alpha, error)
      !> Projected position in bohr, shape (3)
      real(rk), intent(inout) :: point(3)
      !> Closest-point Lagrange multiplier
      real(rk), intent(inout) :: lambda
      !> Unprojected anchor in bohr, shape (3)
      real(rk), intent(in) :: anchor(3)
      !> Nuclear positions in bohr, shape (3, nat)
      real(rk), intent(in) :: centers(:, :)
      !> Sphere radii in bohr, shape (nat)
      real(rk), intent(in) :: radii(:)
      !> SvdW blending sharpness
      real(rk), intent(in) :: blend
      !> Closest-point objective coefficient
      real(rk), intent(in) :: alpha
      !> Reference failure
      type(error_type), allocatable, intent(out) :: error
      real(rk) :: s, g(3), hess(3, 3), jac(4, 4), rhs(4), row(4), factor, delta(4)
      integer :: iter, i, j, pivot
      do iter = 1, 8
         call ref_svdw(point, centers, radii, blend, real(gamma, rk), s, g, hess)
         rhs(1:3) = -alpha*(point - anchor) + lambda*g
         rhs(4) = s
         if (maxval(abs(rhs)) < REFINE_TOL) return
         jac(1:3, 1:3) = -lambda*hess
         do i = 1, 3
            jac(i, i) = jac(i, i) + alpha
         end do
         jac(1:3, 4) = -g; jac(4, 1:3) = -g; jac(4, 4) = 0.0_rk
         do i = 1, 3
            pivot = maxloc(abs(jac(i:4, i)), dim=1) + i - 1
            if (pivot /= i) then
               row = jac(i, :); jac(i, :) = jac(pivot, :); jac(pivot, :) = row
               factor = rhs(i); rhs(i) = rhs(pivot); rhs(pivot) = factor
            end if
            do j = i + 1, 4
               factor = jac(j, i)/jac(i, i)
               jac(j, i:4) = jac(j, i:4) - factor*jac(i, i:4)
               rhs(j) = rhs(j) - factor*rhs(i)
            end do
         end do
         do i = 4, 1, -1
            delta(i) = (rhs(i) - dot_product(jac(i, i + 1:4), delta(i + 1:4)))/jac(i, i)
         end do
         point = point + delta(1:3); lambda = lambda + delta(4)
      end do
      call test_failed(error, "CPCM reference could not refine a converged production branch")
   end subroutine ref_refine

   !> Independent surface values for the retained branches
   !>
   !> @param[in] cavity Production branch seeds, identities, and fixed cavity parameters
   !> @param[in] centers Nuclear positions in bohr, shape (3, nat)
   !> @param[in] anchors Anchor offsets from their owner nuclei in bohr, shape (3, ngrid)
   !> @param[in] blend SvdW blending sharpness
   !> @param[out] points Refined branch positions in bohr, shape (3, ngrid)
   !> @param[out] widths Gaussian widths, shape (ngrid)
   !> @param[out] switches Anchor switching values, shape (ngrid)
   !> @param[out] error Reference failure
   subroutine cpcm_reference_surface(cavity, centers, anchors, blend, points, widths, switches, error)
      !> Production branch seeds, identities, and fixed cavity parameters
      type(cavity_type_drop), intent(in) :: cavity
      !> Nuclear positions in bohr, shape (3, nat)
      real(rk), intent(in) :: centers(:, :)
      !> Anchor offsets from their owner nuclei in bohr, shape (3, ngrid)
      real(rk), intent(in) :: anchors(:, :)
      !> SvdW blending sharpness
      real(rk), intent(in) :: blend
      !> Refined branch positions in bohr, shape (3, ngrid)
      real(rk), intent(out) :: points(:, :)
      !> Gaussian widths, shape (ngrid)
      real(rk), intent(out) :: widths(:)
      !> Anchor switching values, shape (ngrid)
      real(rk), intent(out) :: switches(:)
      !> Reference failure
      type(error_type), allocatable, intent(out) :: error
      real(rk) :: weight(cavity%ngrid), phi(cavity%ngrid), branch(cavity%ngrid)
      real(rk) :: anchor(3), width, rij, rplus, rminus, radius
      integer :: i, a, owner, lo, hi
      ! Full-sum oracle; production screens contributions at tolerance*0.1
      ! The registered fixtures use tolerance=PROJ_TOL; central raw values are
      ! compared separately before any finite-difference derivative check
      select type (lsf => cavity%lsf_model)
      type is (moist_cavity_drop_lsf_svdw_type)
         if (lsf%param%blend_1b /= 1.0_wp .or. lsf%param%blend_2b /= 0.0_wp .or. &
             lsf%param%blend_3b /= 1.0_wp .or. real(lsf%param%blend_k, rk) /= blend .or. &
             lsf%screening_threshold > PROJ_TOL*0.1_wp .or. cavity%param%wleb_prune_level /= 0) then
            call test_failed(error, "CPCM reference requires the tightly screened fixture SvdW parameters")
            return
         end if
      class default
         call test_failed(error, "CPCM reference requires the SvdW level set")
         return
      end select
      do i = 1, cavity%ngrid
         owner = cavity%owner(i)
         anchor = centers(:, owner) + anchors(:, i)
         weight(i) = ref_weight(cavity, i, anchor, centers, blend, points(:, i), error)
         if (allocated(error)) return
         phi(i) = 0.5_rk*real(cavity%param%phi_alpha, rk)*sum((points(:, i) - anchor)**2)
         width = real(cavity%anchor_xi0(i), rk)
         switches(i) = 1.0_rk
         do a = 1, size(centers, 2)
            if (a == owner) cycle
            rij = sqrt(sum((anchor - centers(:, a))**2))
            radius = real(cavity%iswig%radii(a), rk)
            rplus = width*(radius + rij); rminus = width*(radius - rij)
            switches(i) = switches(i)*0.5_rk*(erfc(rplus) + erfc(rminus))
         end do
      end do
      branch = 1.0_rk
      lo = 1
      do while (lo <= cavity%ngrid)
         hi = lo
         do while (hi < cavity%ngrid)
            if (cavity%anchor_id(hi + 1) /= cavity%anchor_id(lo)) exit
            hi = hi + 1
         end do
         if (hi > lo) then
            branch(lo:hi) = exp(-(phi(lo:hi) - minval(phi(lo:hi)))/real(cavity%param%branch_weight_s, rk))
            branch(lo:hi) = branch(lo:hi)/sum(branch(lo:hi))
         end if
         lo = hi + 1
      end do
      do i = 1, cavity%ngrid
         widths(i) = real(cavity%iswig%swx, rk)/(real(cavity%iswig%radii(cavity%owner(i)), rk)*sqrt(weight(i)*branch(i)))
         if (.not. cavity%converged(i)) then
            ! Unconverged points retain their seeds and remain outside the FD mask
            points(:, i) = real(cavity%xyz(:, i), rk)
            widths(i) = real(cavity%xi0(i), rk)
            switches(i) = real(cavity%f(i), rk)
         end if
      end do
      if (.not. all(ieee_is_finite(points)) .or. .not. all(ieee_is_finite(widths)) .or. &
          .not. all(ieee_is_finite(switches)) .or. any(widths <= 0.0_rk) .or. any(switches <= 0.0_rk)) then
         call test_failed(error, "CPCM reference surface contains non-finite or non-positive values")
      end if
   end subroutine cpcm_reference_surface

   !> Independent raw Gaussian PCM matrix values
   !>
   !> @param[in] cavity Production branch seeds, identities, and fixed cavity parameters
   !> @param[in] centers Nuclear positions in bohr, shape (3, nat)
   !> @param[in] anchors Anchor offsets from their owner nuclei in bohr, shape (3, ngrid)
   !> @param[in] blend SvdW blending sharpness
   !> @param[out] values Raw interaction matrix, shape (ngrid, ngrid)
   !> @param[out] error Reference failure
   subroutine cpcm_reference_matrix(cavity, centers, anchors, blend, values, error)
      !> Production branch seeds, identities, and fixed cavity parameters
      type(cavity_type_drop), intent(in) :: cavity
      !> Nuclear positions in bohr, shape (3, nat)
      real(rk), intent(in) :: centers(:, :)
      !> Anchor offsets from their owner nuclei in bohr, shape (3, ngrid)
      real(rk), intent(in) :: anchors(:, :)
      !> SvdW blending sharpness
      real(rk), intent(in) :: blend
      !> Raw interaction matrix, shape (ngrid, ngrid)
      real(rk), intent(out) :: values(:, :)
      !> Reference failure
      type(error_type), allocatable, intent(out) :: error
      real(rk) :: points(3, cavity%ngrid), widths(cavity%ngrid), switches(cavity%ngrid)
      real(rk) :: r2, pair_width
      integer :: i, j
      call cpcm_reference_surface(cavity, centers, anchors, blend, points, widths, switches, error)
      if (allocated(error)) return
      do i = 1, cavity%ngrid
         values(i, i) = sqrt(2.0_rk/acos(-1.0_rk))*widths(i)/switches(i)
      end do
      do j = 1, cavity%ngrid
         do i = 1, j - 1
            r2 = sum((points(:, i) - points(:, j))**2)
            pair_width = widths(i)*widths(j)/sqrt(widths(i)**2 + widths(j)**2)
            if (r2 == 0.0_rk) then
               values(i, j) = 2.0_rk/sqrt(acos(-1.0_rk))*pair_width
            else
               values(i, j) = erf(pair_width*sqrt(r2))/sqrt(r2)
            end if
            values(j, i) = values(i, j)
         end do
      end do
   end subroutine cpcm_reference_matrix

end module test_cavity_drop_cpcm
