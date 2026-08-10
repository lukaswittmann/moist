module test_cavity_drop_cpcm
   use mctc_env_accuracy, only: wp
   use mctc_env_error, only: mctc_error => error_type
   use mctc_io, only: structure_type, new
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use testdrive, only: to_string, test_failed
   use moist_cavity_drop, only: cavity_type_drop, new_cavity_drop
   use moist_cavity_drop_lsf_svdw, only: moist_cavity_drop_lsf_svdw_type
   use moist_radii, only: default_cpcm_radii
   use moist_data_radii_legacy, only: get_radius_func
   use mstore, only: get_structure
   use moist_math_lapack, only: getrf, getri
   use moist_model_component_pcm_amat, only: assemble_pcm_amat, &
      & assemble_pcm_amat_with_gradient, pcm_amat_surface_weights, &
      & pcm_amat_nuclear_gradient
   use moist_model_component_pcm_electrostatics, only: &
      & pcm_electrostatic_nuclear_gradient
   use moist_context, only: moist_context_type, new_context
   use, intrinsic :: iso_fortran_env, only: error_unit
   implicit none (type, external)
   private

   public :: collect_cavity_drop_cpcm

   integer, parameter :: ndim = 3

   real(wp), parameter :: k = 1.5
   real(wp), parameter :: gamma = 1.0
   integer, parameter :: NUM_LEB = 26

   real(wp), parameter :: STEP_SIZE = 1.0E-4_wp

   real(wp), parameter :: ATHR_OFFDIAG = 1.0E-9_wp
   real(wp), parameter :: RTHR_OFFDIAG = 1.0E-8_wp

   real(wp), parameter :: ATHR_DIAG = 1.0E-6_wp
   real(wp), parameter :: RTHR_DIAG = 2.0E-6_wp

   !> Fraction of the quantities largest derivative
   real(wp), parameter :: MTHR_OFFDIAG = 2.0E-10_wp
   real(wp), parameter :: MTHR_DIAG = 2.0E-10_wp

   real(wp), parameter :: PROJ_TOL = 1E-14_wp
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
                  new_unittest("heavy28_h2o", test_heavy28_h2o), &
                  new_unittest("mb16_43_01", test_mb16_43_01), &
                  new_unittest("mb16_43_19", test_mb16_43_19), &
                  new_unittest("but14diol_1", test_but14diol_1), &
                  new_unittest("but14diol_32", test_but14diol_32), &
                  new_unittest("il16_008", test_il16_008) &
                  ]
   end subroutine collect_cavity_drop_cpcm

   !> Test the contracted A-matrix gradient against the explicit tensor
   !> contraction of the dense derivative built by
   !> `assemble_pcm_amat_with_gradient`.
   subroutine test_contract_amat1_q1q2_rA(error)
      type(error_type), allocatable, intent(out) :: error

      type(structure_type) :: mol
      type(cavity_type_drop), allocatable :: cavity
      real(wp), allocatable :: radii(:)
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

      call fill_cpcm_radii(mol, radii, error)
      if (allocated(error)) return

      allocate (cavity)
      block
         type(moist_cavity_drop_lsf_svdw_type) :: svdw_template
         call svdw_template%new(blend_k=k, blend_3b=gamma)
         call new_cavity_drop(cavity, ctx, nleb=NUM_LEB, &
                              tolerance=PROJ_TOL, proj_maxiter=PROJ_MAXITER, proj_level=PROJ_LEVEL, &
                              radius_model=default_cpcm_radii(), &
                              lsf_model=svdw_template, error=cavity_error)
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

      call cavity%get_gradient()
      if (allocated(cavity%error)) then
         call test_failed(error, cavity%error%message)
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

   !> Test the fused nuclear/electronic contraction against finite differences.
   !> This checks the point-charge case with qefield = 0.
   subroutine test_contract_nuc_elec_pointcharge_fd(error)
      type(error_type), allocatable, intent(out) :: error

      type(structure_type) :: mol, mol_fd
      type(cavity_type_drop), allocatable :: cavity
      real(wp), allocatable :: radii(:)
      real(wp), allocatable :: surface_q(:), qefield(:, :), za(:)
      real(wp), allocatable :: grad_ctr(:, :), grad_num(:, :)
      integer, allocatable :: numbering_ref(:)
      integer :: iat, iaxis, igrid, ngrid
      real(wp) :: e_plus, e_minus
      real(wp), parameter :: step = 1.0e-5_wp
      type(mctc_error), allocatable :: cavity_error
      !> Local run context borrowed by the cavities built here
      type(moist_context_type), target :: ctx

      call new_context(ctx, verbosity=0)

      call get_structure(mol, "MB16-43", "15")

      call fill_cpcm_radii(mol, radii, error)
      if (allocated(error)) return

      allocate (cavity)
      block
         type(moist_cavity_drop_lsf_svdw_type) :: svdw_template
         call svdw_template%new(blend_k=k, blend_3b=gamma)
         call new_cavity_drop(cavity, ctx, nleb=NUM_LEB, &
                              tolerance=PROJ_TOL, proj_maxiter=PROJ_MAXITER, proj_level=PROJ_LEVEL, &
                              wleb_prune_level=3, &
                              radius_model=default_cpcm_radii(), &
                              lsf_model=svdw_template, error=cavity_error)
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

      call cavity%get_gradient()
      if (allocated(cavity%error)) then
         call test_failed(error, cavity%error%message)
         return
      end if

      allocate (surface_q(ngrid), qefield(3, ngrid), za(mol%nat))
      allocate (grad_ctr(3, mol%nat), grad_num(3, mol%nat))

      do igrid = 1, ngrid
         if (mod(igrid, 2) == 0) then
            surface_q(igrid) = -0.07_wp/(real(igrid, wp) + 0.5_wp)
         else
            surface_q(igrid) = 0.09_wp/(real(igrid, wp) + 0.25_wp)
         end if
      end do

      qefield = 0.0_wp
      do iat = 1, mol%nat
         za(iat) = real(mol%num(mol%id(iat)), wp)
      end do

      call pcm_electrostatic_nuclear_gradient(cavity%xyz, cavity%sphxyz, &
                                              cavity%xyz1_rA, surface_q, qefield, za, &
                                              grad_ctr, cavity_error)
      if (allocated(cavity_error)) then
         call test_failed(error, cavity_error%message)
         return
      end if

      grad_num = 0.0_wp
      do iat = 1, mol%nat
         do iaxis = 1, 3
            mol_fd = mol
            mol_fd%xyz(iaxis, iat) = mol_fd%xyz(iaxis, iat) + step
            call cavity%update(mol_fd, error=cavity_error)
            if (allocated(cavity_error)) then
               call test_failed(error, cavity_error%message)
               return
            end if
            if (cavity%ngrid /= ngrid) then
               call test_failed(error, "contract_nuc_elec FD: ngrid changed for +step")
               return
            end if
            if (any(cavity%numbering(1:ngrid) /= numbering_ref)) then
               call test_failed(error, "contract_nuc_elec FD: numbering changed for +step")
               return
            end if
            e_plus = weighted_nuclear_potential(cavity, mol_fd, surface_q, za)

            mol_fd = mol
            mol_fd%xyz(iaxis, iat) = mol_fd%xyz(iaxis, iat) - step
            call cavity%update(mol_fd, error=cavity_error)
            if (allocated(cavity_error)) then
               call test_failed(error, cavity_error%message)
               return
            end if
            if (cavity%ngrid /= ngrid) then
               call test_failed(error, "contract_nuc_elec FD: ngrid changed for -step")
               return
            end if
            if (any(cavity%numbering(1:ngrid) /= numbering_ref)) then
               call test_failed(error, "contract_nuc_elec FD: numbering changed for -step")
               return
            end if
            e_minus = weighted_nuclear_potential(cavity, mol_fd, surface_q, za)

            grad_num(iaxis, iat) = (e_plus - e_minus)/(2.0_wp*step)
         end do
      end do

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
                       thr_abs=2.0e-6_wp, thr_rel=2.0e-5_wp, &
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
      real(wp), allocatable :: radii(:)

      ! Create single oxygen atom
      call new(mol, [8], reshape([0.0_wp, 0.0_wp, 0.0_wp], [3, 1]))

      call fill_cpcm_radii(mol, radii, error)
      if (allocated(error)) return

      call do_test(error, mol, radii)
   end subroutine test_single_atom

   !> Test A matrix gradient for a dimer
   subroutine test_dimer(error)
      type(error_type), allocatable, intent(out) :: error
      type(structure_type) :: mol
      real(wp), allocatable :: radii(:)

      ! Create dimer (two oxygen atoms)
      call new(mol, [8, 8], reshape([0.0_wp, 0.0_wp, 0.0_wp, &
                                     3.0_wp, 0.0_wp, 0.0_wp], [3, 2]))

      call fill_cpcm_radii(mol, radii, error)
      if (allocated(error)) return

      call do_test(error, mol, radii)
   end subroutine test_dimer

   !> Test A matrix gradient for 5-argon geometry with custom blend-k
   subroutine test_ar5_blendk_09(error)
      type(error_type), allocatable, intent(out) :: error
      type(structure_type) :: mol
      real(wp), allocatable :: radii(:)

      call new(mol, [18, 18, 18, 18, 18], reshape([ &
                                                  0.2_wp, 0.0_wp, 5.1_wp, &
                                                  -2.2_wp, -2.2_wp, 0.0_wp, &
                                                  2.2_wp, -2.2_wp, 0.0_wp, &
                                                  -2.2_wp, 2.2_wp, 0.0_wp, &
                                                  2.2_wp, 2.2_wp, 0.0_wp], [3, 5]))

      call fill_cpcm_radii(mol, radii, error)
      if (allocated(error)) return

      call do_test(error, mol, radii, blend_k_override=0.8_wp)
   end subroutine test_ar5_blendk_09

   !> Test A matrix gradient for MB16-43 h2
   subroutine test_mb16_43_h2(error)
      type(error_type), allocatable, intent(out) :: error
      type(structure_type) :: mol
      real(wp), allocatable :: radii(:)

      call get_structure(mol, "MB16-43", "H2")

      call fill_cpcm_radii(mol, radii, error)
      if (allocated(error)) return

      call do_test(error, mol, radii)
   end subroutine test_mb16_43_h2

   !> Test A matrix gradient for bih3_h2o system
   subroutine test_bih3_h2o(error)
      type(error_type), allocatable, intent(out) :: error
      type(structure_type) :: mol
      real(wp), allocatable :: radii(:)

      call get_structure(mol, "Heavy28", "bih3_h2o")

      call fill_cpcm_radii(mol, radii, error)
      if (allocated(error)) return

      call do_test(error, mol, radii)
   end subroutine test_bih3_h2o

   !> Test A matrix gradient for Heavy28 h2o
   subroutine test_heavy28_h2o(error)
      type(error_type), allocatable, intent(out) :: error
      type(structure_type) :: mol
      real(wp), allocatable :: radii(:)

      call get_structure(mol, "Heavy28", "h2o")

      call fill_cpcm_radii(mol, radii, error)
      if (allocated(error)) return

      call do_test(error, mol, radii)
   end subroutine test_heavy28_h2o

   !> Test A matrix gradient for Heavy28 pbh4
   subroutine test_heavy28_pbh4(error)
      type(error_type), allocatable, intent(out) :: error
      type(structure_type) :: mol
      real(wp), allocatable :: radii(:)

      call get_structure(mol, "Heavy28", "pbh4")

      call fill_cpcm_radii(mol, radii, error)
      if (allocated(error)) return

      call do_test(error, mol, radii)
   end subroutine test_heavy28_pbh4

   !> Test A matrix gradient for MB16-43 01
   subroutine test_mb16_43_01(error)
      type(error_type), allocatable, intent(out) :: error
      type(structure_type) :: mol
      real(wp), allocatable :: radii(:)

      call get_structure(mol, "MB16-43", "01")

      call fill_cpcm_radii(mol, radii, error)
      if (allocated(error)) return

      call do_test(error, mol, radii)
   end subroutine test_mb16_43_01

   !> Test A matrix gradient for MB16-43 19
   subroutine test_mb16_43_19(error)
      type(error_type), allocatable, intent(out) :: error
      type(structure_type) :: mol
      real(wp), allocatable :: radii(:)

      call get_structure(mol, "MB16-43", "19")

      call fill_cpcm_radii(mol, radii, error)
      if (allocated(error)) return

      call do_test(error, mol, radii)
   end subroutine test_mb16_43_19

   !> Test A matrix gradient for But14diol 1
   subroutine test_but14diol_1(error)
      type(error_type), allocatable, intent(out) :: error
      type(structure_type) :: mol
      real(wp), allocatable :: radii(:)

      call get_structure(mol, "But14diol", "1")

      call fill_cpcm_radii(mol, radii, error)
      if (allocated(error)) return

      call do_test(error, mol, radii)
   end subroutine test_but14diol_1

   !> Test A matrix gradient for But14diol 32
   subroutine test_but14diol_32(error)
      type(error_type), allocatable, intent(out) :: error
      type(structure_type) :: mol
      real(wp), allocatable :: radii(:)

      call get_structure(mol, "But14diol", "32")

      call fill_cpcm_radii(mol, radii, error)
      if (allocated(error)) return

      call do_test(error, mol, radii)
   end subroutine test_but14diol_32

   !> Test A matrix gradient for IL16 008
   subroutine test_il16_008(error)
      type(error_type), allocatable, intent(out) :: error
      type(structure_type) :: mol
      real(wp), allocatable :: radii(:)

      call get_structure(mol, "IL16", "008")

      call fill_cpcm_radii(mol, radii, error)
      if (allocated(error)) return

      call do_test(error, mol, radii)
   end subroutine test_il16_008

   !> Test A matrix gradient w.r.t. atomic positions
   subroutine do_test(error, mol, radii, blend_k_override)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      !> Molecular structure
      type(structure_type), intent(in) :: mol
      !> Atomic radii
      real(wp), intent(in) :: radii(:)
      !> Optional override for blending steepness parameter
      real(wp), intent(in), optional :: blend_k_override

      type(structure_type) :: mol_fd
      type(cavity_type_drop), allocatable :: cavity
      real(wp), allocatable :: Amat0(:, :), Amat1_rA(:, :, :, :)
      real(wp), allocatable :: en_Amat1_rA(:, :, :, :)
      real(wp), allocatable :: num_Amat1_rA(:, :, :, :)
      real(wp), allocatable :: nn_Amat(:, :), n_Amat(:, :)
      real(wp), allocatable :: p_Amat(:, :), pp_Amat(:, :)
      integer :: iat, idir, igrid, jgrid, ngrid, num_idn, num_jdn
      integer :: max_numbering, idx_i, idx_j
      real(wp), allocatable :: ref_f(:)
      integer, allocatable :: ref_numbering(:), numbering_to_idx(:)
      logical, allocatable :: valid_gridpoint_ref(:)
      logical, allocatable :: valid_gridpoint(:, :, :)
      logical, allocatable :: nn_converged(:, :, :), n_converged(:, :, :), &
                              p_converged(:, :, :), pp_converged(:, :, :)
      integer, allocatable :: num_nn(:, :, :), num_n(:, :, :), num_p(:, :, :), num_pp(:, :, :)

      real(wp), allocatable :: Amat_copy(:, :), Amat_inv(:, :), identity_test(:, :)
      integer, allocatable :: ipiv(:)
      integer :: info
      real(wp) :: max_err
      real(wp) :: blend_k_local
      !> Analytic and numeric value of the entry under inspection
      real(wp) :: analytic, numeric
      !> Per-channel scale, max |numeric| over the channel; sets the MTHR term
      real(wp) :: off_scale, diag_scale
      !> Worst |deviation|/tolerance seen, and where it occurred
      real(wp) :: worst_ratio, worst_a, worst_n, worst_tol
      integer :: worst_iat, worst_idir, worst_i, worst_j
      character(len=16) :: worst_kind
      !> Tolerance applied to the entry under inspection
      real(wp) :: tol
      type(mctc_error), allocatable :: cavity_error
      type(moist_context_type), target :: ctx

      !> Initialize cavity with 26-point Lebedev grid
      blend_k_local = k
      off_scale = 0.0_wp; diag_scale = 0.0_wp
      worst_ratio = 0.0_wp; worst_a = 0.0_wp; worst_n = 0.0_wp; worst_tol = 0.0_wp
      worst_iat = 0; worst_idir = 0; worst_i = 0; worst_j = 0
      worst_kind = "none"
      if (present(blend_k_override)) blend_k_local = blend_k_override
      allocate (cavity)
      block
         type(moist_cavity_drop_lsf_svdw_type) :: svdw_template
         call svdw_template%new(blend_k=blend_k_local, blend_3b=gamma)
         call new_context(ctx, verbosity=0, debug=.false.)
         call new_cavity_drop(cavity, ctx, nleb=NUM_LEB, &
                              tolerance=PROJ_TOL, proj_maxiter=PROJ_MAXITER, proj_level=PROJ_LEVEL, &
                              radius_model=default_cpcm_radii(), &
                              lsf_model=svdw_template, error=cavity_error)
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

      allocate (ref_numbering(ngrid))
      if (ngrid > 0) then
         ref_numbering = cavity%numbering(1:ngrid)
         max_numbering = maxval(ref_numbering)
      else
         max_numbering = 0
      end if
      allocate (numbering_to_idx(max(1, max_numbering)), source=0)
      do igrid = 1, ngrid
         num_idn = ref_numbering(igrid)
         if (num_idn > 0 .and. num_idn <= size(numbering_to_idx)) then
            numbering_to_idx(num_idn) = igrid
         end if
      end do

      ! Store reference convergence status before finite differences
      allocate (valid_gridpoint_ref(ngrid), source=.false.)
      do jgrid = 1, cavity%ngrid
         num_idn = cavity%numbering(jgrid)
         if (num_idn <= 0 .or. num_idn > size(numbering_to_idx)) cycle
         idx_i = numbering_to_idx(num_idn)
         if (idx_i <= 0) cycle
         valid_gridpoint_ref(idx_i) = cavity%converged(jgrid)
      end do

      !> Get analytic geometry derivatives.
      call cavity%get_gradient()

      !> Assemble A matrix and its gradient
      allocate (Amat0(ngrid, ngrid), Amat1_rA(ndim, mol%nat, ngrid, ngrid))
      call assemble_pcm_amat_with_gradient(cavity%xi0, cavity%f, cavity%xyz, &
                                           cavity%xi1_rA, cavity%f1_rA, cavity%xyz1_rA, &
                                           Amat0, Amat1_rA, cavity_error)
      if (allocated(cavity_error)) then
         call test_failed(error, "assemble_pcm_amat_with_gradient failed: "//cavity_error%message)
         return
      end if

      !> Allocate arrays for numerical and analytical gradients (mapped to reference grid IDs)
      allocate (num_Amat1_rA(ndim, mol%nat, ngrid, ngrid), source=0.0_wp)
      allocate (en_Amat1_rA(ndim, mol%nat, ngrid, ngrid), source=0.0_wp)

      !> Store reference grid properties for diagnostics
      allocate (ref_f(ngrid), source=0.0_wp)
      do igrid = 1, cavity%ngrid
         num_idn = cavity%numbering(igrid)
         if (num_idn <= 0 .or. num_idn > size(numbering_to_idx)) cycle
         idx_i = numbering_to_idx(num_idn)
         if (idx_i <= 0) cycle
         ref_f(idx_i) = cavity%f(igrid)
      end do

      ! Remap analytical gradient using numbering (do this BEFORE FD loop while cavity is at reference)
      do iat = 1, mol%nat
         do idir = 1, ndim
            do igrid = 1, cavity%ngrid
               num_idn = cavity%numbering(igrid)
               if (num_idn <= 0 .or. num_idn > size(numbering_to_idx)) cycle
               idx_i = numbering_to_idx(num_idn)
               if (idx_i <= 0) cycle
               do jgrid = 1, cavity%ngrid
                  num_jdn = cavity%numbering(jgrid)
                  if (num_jdn <= 0 .or. num_jdn > size(numbering_to_idx)) cycle
                  idx_j = numbering_to_idx(num_jdn)
                  if (idx_j <= 0) cycle
                  ! Store analytical gradient (mapped by numbering)
                  en_Amat1_rA(idir, iat, idx_i, idx_j) = Amat1_rA(idir, iat, igrid, jgrid)
               end do
            end do
         end do
      end do

      !> Test inversion of A matrix using LAPACK

      ! Copy A matrix for inversion (getrf/getri overwrite the input)
      allocate (Amat_copy(ngrid, ngrid), source=Amat0(1:ngrid, 1:ngrid))
      allocate (Amat_inv(ngrid, ngrid))
      allocate (identity_test(ngrid, ngrid))
      allocate (ipiv(ngrid))

      ! LU factorization
      call getrf(Amat_copy, ipiv, info)
      if (info /= 0) then
         call test_failed(error, "LAPACK getrf failed with info = "//to_string(info))
         return
      end if

      ! Compute inverse from LU factorization
      call getri(Amat_copy, ipiv, info)
      if (info /= 0) then
         call test_failed(error, "LAPACK getri failed with info = "//to_string(info))
         return
      end if

      ! Amat_copy now contains A^{-1}
      Amat_inv = Amat_copy

      ! Test: A * A^{-1} should be identity
      ! Use explicit matrix multiplication
      identity_test = matmul(Amat0, Amat_inv)

      ! Check that result is close to identity matrix
      max_err = 0.0_wp
      do igrid = 1, ngrid
         do jgrid = 1, ngrid
            if (igrid == jgrid) then
               max_err = max(max_err, abs(identity_test(igrid, jgrid) - 1.0_wp))
            else
               max_err = max(max_err, abs(identity_test(igrid, jgrid)))
            end if
         end do
      end do

      ! Check that max error is small (A * A^{-1} = I)
      if (max_err >= 1.0e-10_wp) then
         write (error_unit, "(A)") "A * A^{-1} differs from identity. Problematic entries:"
         do igrid = 1, ngrid
            do jgrid = 1, ngrid
               if (igrid == jgrid) then
                  if (abs(identity_test(igrid, jgrid) - 1.0_wp) > 1.0e-10_wp) then
                     write (error_unit, "(A,I4,A,I4,A,ES15.6,A,ES15.6)") &
                        "  Diagonal (", igrid, ",", jgrid, "): got ", &
                        identity_test(igrid, jgrid), ", expected 1.0, err = ", &
                        abs(identity_test(igrid, jgrid) - 1.0_wp)
                  end if
               else
                  if (abs(identity_test(igrid, jgrid)) > 1.0e-10_wp) then
                     write (error_unit, "(A,I4,A,I4,A,ES15.6,A,ES15.6)") &
                        "  Off-diag (", igrid, ",", jgrid, "): got ", &
                        identity_test(igrid, jgrid), ", expected 0.0, err = ", &
                        abs(identity_test(igrid, jgrid))
                  end if
               end if
            end do
         end do
         call test_failed(error, "A * A^{-1} differs from identity by "//to_string(max_err))
         return
      end if

      !> Compute numerical gradients via finite differences

      ! Allocate temporary arrays for FD bookkeeping
      allocate (nn_Amat(ngrid, ngrid), source=0.0_wp)
      allocate (n_Amat(ngrid, ngrid), source=0.0_wp)
      allocate (p_Amat(ngrid, ngrid), source=0.0_wp)
      allocate (pp_Amat(ngrid, ngrid), source=0.0_wp)

      ! Allocate convergence and numbering tracking arrays
      allocate (nn_converged(ndim, mol%nat, ngrid), source=.false.)
      allocate (n_converged(ndim, mol%nat, ngrid), source=.false.)
      allocate (p_converged(ndim, mol%nat, ngrid), source=.false.)
      allocate (pp_converged(ndim, mol%nat, ngrid), source=.false.)
      allocate (num_nn(ndim, mol%nat, ngrid), source=0)
      allocate (num_n(ndim, mol%nat, ngrid), source=0)
      allocate (num_p(ndim, mol%nat, ngrid), source=0)
      allocate (num_pp(ndim, mol%nat, ngrid), source=0)

      !> Compute numerical gradients via finite differences with bookkeeping
      do iat = 1, mol%nat
         do idir = 1, ndim
            ! Initialize arrays to zero
            nn_Amat = 0.0_wp
            n_Amat = 0.0_wp
            p_Amat = 0.0_wp
            pp_Amat = 0.0_wp

            ! -2h step
            mol_fd = mol
            mol_fd%xyz(idir, iat) = mol_fd%xyz(idir, iat) - 2.0_wp*STEP_SIZE
            call cavity%update(mol_fd, error=cavity_error)
            if (allocated(cavity_error)) then
               call test_failed(error, cavity_error%message)
               return
            end if
            call assemble_pcm_amat(cavity%xi0, cavity%f, cavity%xyz, Amat0, cavity_error)
            if (allocated(cavity_error)) then
               call test_failed(error, "assemble_pcm_amat failed: "//cavity_error%message)
               return
            end if
            ! Store using numbering for bookkeeping
            do igrid = 1, cavity%ngrid
               num_idn = cavity%numbering(igrid)
               if (num_idn <= 0 .or. num_idn > size(numbering_to_idx)) cycle
               idx_i = numbering_to_idx(num_idn)
               if (idx_i <= 0) cycle
               do jgrid = 1, cavity%ngrid
                  num_jdn = cavity%numbering(jgrid)
                  if (num_jdn <= 0 .or. num_jdn > size(numbering_to_idx)) cycle
                  idx_j = numbering_to_idx(num_jdn)
                  if (idx_j <= 0) cycle
                  nn_Amat(idx_i, idx_j) = Amat0(igrid, jgrid)
               end do
               nn_converged(idir, iat, idx_i) = cavity%converged(igrid)
               num_nn(idir, iat, idx_i) = num_idn
            end do

            ! -1h step
            mol_fd = mol
            mol_fd%xyz(idir, iat) = mol_fd%xyz(idir, iat) - 1.0_wp*STEP_SIZE
            call cavity%update(mol_fd, error=cavity_error)
            if (allocated(cavity_error)) then
               call test_failed(error, cavity_error%message)
               return
            end if
            call assemble_pcm_amat(cavity%xi0, cavity%f, cavity%xyz, Amat0, cavity_error)
            if (allocated(cavity_error)) then
               call test_failed(error, "assemble_pcm_amat failed: "//cavity_error%message)
               return
            end if
            do igrid = 1, cavity%ngrid
               num_idn = cavity%numbering(igrid)
               if (num_idn <= 0 .or. num_idn > size(numbering_to_idx)) cycle
               idx_i = numbering_to_idx(num_idn)
               if (idx_i <= 0) cycle
               do jgrid = 1, cavity%ngrid
                  num_jdn = cavity%numbering(jgrid)
                  if (num_jdn <= 0 .or. num_jdn > size(numbering_to_idx)) cycle
                  idx_j = numbering_to_idx(num_jdn)
                  if (idx_j <= 0) cycle
                  n_Amat(idx_i, idx_j) = Amat0(igrid, jgrid)
               end do
               n_converged(idir, iat, idx_i) = cavity%converged(igrid)
               num_n(idir, iat, idx_i) = num_idn
            end do

            ! +1h step
            mol_fd = mol
            mol_fd%xyz(idir, iat) = mol_fd%xyz(idir, iat) + 1.0_wp*STEP_SIZE
            call cavity%update(mol_fd, error=cavity_error)
            if (allocated(cavity_error)) then
               call test_failed(error, cavity_error%message)
               return
            end if
            call assemble_pcm_amat(cavity%xi0, cavity%f, cavity%xyz, Amat0, cavity_error)
            if (allocated(cavity_error)) then
               call test_failed(error, "assemble_pcm_amat failed: "//cavity_error%message)
               return
            end if
            do igrid = 1, cavity%ngrid
               num_idn = cavity%numbering(igrid)
               if (num_idn <= 0 .or. num_idn > size(numbering_to_idx)) cycle
               idx_i = numbering_to_idx(num_idn)
               if (idx_i <= 0) cycle
               do jgrid = 1, cavity%ngrid
                  num_jdn = cavity%numbering(jgrid)
                  if (num_jdn <= 0 .or. num_jdn > size(numbering_to_idx)) cycle
                  idx_j = numbering_to_idx(num_jdn)
                  if (idx_j <= 0) cycle
                  p_Amat(idx_i, idx_j) = Amat0(igrid, jgrid)
               end do
               p_converged(idir, iat, idx_i) = cavity%converged(igrid)
               num_p(idir, iat, idx_i) = num_idn
            end do

            ! +2h step
            mol_fd = mol
            mol_fd%xyz(idir, iat) = mol_fd%xyz(idir, iat) + 2.0_wp*STEP_SIZE
            call cavity%update(mol_fd, error=cavity_error)
            if (allocated(cavity_error)) then
               call test_failed(error, cavity_error%message)
               return
            end if
            call assemble_pcm_amat(cavity%xi0, cavity%f, cavity%xyz, Amat0, cavity_error)
            if (allocated(cavity_error)) then
               call test_failed(error, "assemble_pcm_amat failed: "//cavity_error%message)
               return
            end if
            do igrid = 1, cavity%ngrid
               num_idn = cavity%numbering(igrid)
               if (num_idn <= 0 .or. num_idn > size(numbering_to_idx)) cycle
               idx_i = numbering_to_idx(num_idn)
               if (idx_i <= 0) cycle
               do jgrid = 1, cavity%ngrid
                  num_jdn = cavity%numbering(jgrid)
                  if (num_jdn <= 0 .or. num_jdn > size(numbering_to_idx)) cycle
                  idx_j = numbering_to_idx(num_jdn)
                  if (idx_j <= 0) cycle
                  pp_Amat(idx_i, idx_j) = Amat0(igrid, jgrid)
               end do
               pp_converged(idir, iat, idx_i) = cavity%converged(igrid)
               num_pp(idir, iat, idx_i) = num_idn
            end do

            ! Central difference formula: f'(x) ~ [-f(x+2h) + 8f(x+h) - 8f(x-h) + f(x-2h)] / (12h)
            do igrid = 1, ngrid
               do jgrid = 1, ngrid
                  ! Skip if any of the 4 FD steps are missing for either gridpoint
                  if (num_nn(idir, iat, igrid) == 0 .or. num_nn(idir, iat, jgrid) == 0) cycle
                  if (num_n(idir, iat, igrid) == 0 .or. num_n(idir, iat, jgrid) == 0) cycle
                  if (num_p(idir, iat, igrid) == 0 .or. num_p(idir, iat, jgrid) == 0) cycle
                  if (num_pp(idir, iat, igrid) == 0 .or. num_pp(idir, iat, jgrid) == 0) cycle

                  ! Compute numerical gradient
                  num_Amat1_rA(idir, iat, igrid, jgrid) = &
                     (-pp_Amat(igrid, jgrid) + 8.0_wp*p_Amat(igrid, jgrid) &
                      - 8.0_wp*n_Amat(igrid, jgrid) + nn_Amat(igrid, jgrid)) &
                     /(12.0_wp*STEP_SIZE)
               end do
            end do
         end do
      end do

      ! Build comprehensive validity tracking array
      ! A gridpoint pair (igrid, jgrid) is valid only if:
      ! 1. Both gridpoints converged in the reference configuration
      ! 2. Both exist (numbering > 0) and converged in all 4 FD steps
      allocate (valid_gridpoint(ndim, mol%nat, ngrid), source=.false.)

      do iat = 1, mol%nat
         do idir = 1, ndim
            do igrid = 1, ngrid
               ! Start with reference convergence
               valid_gridpoint(idir, iat, igrid) = valid_gridpoint_ref(igrid)

               ! Check all 4 FD steps: must exist (numbering > 0) AND converged
               if (valid_gridpoint(idir, iat, igrid)) then
                  valid_gridpoint(idir, iat, igrid) = &
                     (num_nn(idir, iat, igrid) > 0 .and. nn_converged(idir, iat, igrid)) .and. &
                     (num_n(idir, iat, igrid) > 0 .and. n_converged(idir, iat, igrid)) .and. &
                     (num_p(idir, iat, igrid) > 0 .and. p_converged(idir, iat, igrid)) .and. &
                     (num_pp(idir, iat, igrid) > 0 .and. pp_converged(idir, iat, igrid))
               end if
            end do
         end do
      end do

      !> Compare analytic vs numeric A matrix gradients (only valid gridpoint
      !> pairs). Two passes over each channel: the first fixes its scale, which
      !> the scale-relative MTHR term of the tolerance needs, the second judges
      !> every entry against that tolerance and keeps the worst violation.
      !> The passes are plain traversals of arrays that are already complete, so
      !> they cost nothing next to the finite differences that filled them.

      !> Off-diagonal elements, pass 1: channel scale.
      do iat = 1, mol%nat
         do idir = 1, ndim
            do igrid = 1, ngrid
               do jgrid = 1, ngrid
                  if (jgrid == igrid) cycle
                  if (.not. pair_is_valid(idir, iat, igrid, jgrid)) cycle
                  off_scale = max(off_scale, abs(num_Amat1_rA(idir, iat, igrid, jgrid)))
               end do
            end do
         end do
      end do

      !> Off-diagonal elements, pass 2: tolerance.
      do iat = 1, mol%nat
         do idir = 1, ndim
            do igrid = 1, ngrid
               do jgrid = 1, ngrid
                  if (jgrid == igrid) cycle
                  if (.not. pair_is_valid(idir, iat, igrid, jgrid)) cycle

                  analytic = en_Amat1_rA(idir, iat, igrid, jgrid)
                  numeric = num_Amat1_rA(idir, iat, igrid, jgrid)
                  tol = max(ATHR_OFFDIAG, RTHR_OFFDIAG*abs(numeric), &
                            MTHR_OFFDIAG*off_scale)
                  call record_worst("off-diagonal", analytic, numeric, tol, &
                                    iat, idir, igrid, jgrid)
               end do
            end do
         end do
      end do

      !> Diagonal elements, pass 1: channel scale.
      do iat = 1, mol%nat
         do idir = 1, ndim
            do igrid = 1, ngrid
               if (.not. pair_is_valid(idir, iat, igrid, igrid)) cycle
               ! Cycle if switching function (otherwise num. noise dominates as A_ii ~ 1/f)
               if (ref_f(igrid) < 1.0e-4_wp) cycle
               diag_scale = max(diag_scale, abs(num_Amat1_rA(idir, iat, igrid, igrid)))
            end do
         end do
      end do

      !> Diagonal elements, pass 2: tolerance.
      do iat = 1, mol%nat
         do idir = 1, ndim
            do igrid = 1, ngrid
               if (.not. pair_is_valid(idir, iat, igrid, igrid)) cycle
               if (ref_f(igrid) < 1.0e-4_wp) cycle

               analytic = en_Amat1_rA(idir, iat, igrid, igrid)
               numeric = num_Amat1_rA(idir, iat, igrid, igrid)
               tol = max(ATHR_DIAG, RTHR_DIAG*abs(numeric), MTHR_DIAG*diag_scale)
               call record_worst("diagonal", analytic, numeric, tol, &
                                 iat, idir, igrid, igrid)
            end do
         end do
      end do

      ! One assertion at the end rather than a check per entry, so the message
      ! can name the entry that actually decided the outcome.
      if (worst_ratio > 1.0_wp) then
         write (error_unit, "(a)") "A-matrix gradient exceeds its tolerance:"
         write (error_unit, "(2x,a,a)") "channel   : ", trim(worst_kind)
         write (error_unit, "(2x,a,i0,a,i0)") "atom/axis : ", worst_iat, " / ", worst_idir
         write (error_unit, "(2x,a,i0,a,i0)") "entry     : ", worst_i, " , ", worst_j
         write (error_unit, "(2x,a,es15.6)") "analytic  : ", worst_a
         write (error_unit, "(2x,a,es15.6)") "numeric   : ", worst_n
         write (error_unit, "(2x,a,es15.6)") "deviation : ", abs(worst_a - worst_n)
         write (error_unit, "(2x,a,es15.6)") "tolerance : ", worst_tol
         call test_failed(error, "A-matrix gradient off by "// &
                          to_string(worst_ratio)//" times its tolerance")
         return
      end if

   contains

      !> True iff a gridpoint pair carries a usable finite difference.
      !>
      !> Both points must have converged in the reference configuration and
      !> must have existed (nonzero numbering) and converged in all four
      !> displaced ones; otherwise `num_Amat1_rA` was left at zero and there is
      !> nothing to compare against.
      !> @param[in] idir   Displacement axis
      !> @param[in] iat    Displaced atom
      !> @param[in] ig     Row gridpoint
      !> @param[in] jg     Column gridpoint
      !> @return    valid  Whether the pair may be compared
      pure logical function pair_is_valid(idir, iat, ig, jg) result(valid)
         !> Displacement axis and displaced atom
         integer, intent(in) :: idir, iat
         !> Row and column gridpoints
         integer, intent(in) :: ig, jg

         valid = valid_gridpoint(idir, iat, ig) .and. valid_gridpoint(idir, iat, jg)
         if (.not. valid) return
         valid = num_nn(idir, iat, ig) /= 0 .and. num_nn(idir, iat, jg) /= 0 .and. &
                 num_n(idir, iat, ig) /= 0 .and. num_n(idir, iat, jg) /= 0 .and. &
                 num_p(idir, iat, ig) /= 0 .and. num_p(idir, iat, jg) /= 0 .and. &
                 num_pp(idir, iat, ig) /= 0 .and. num_pp(idir, iat, jg) /= 0
      end function pair_is_valid

      !> Keep the entry with the largest deviation measured in its own tolerance.
      !> @param[in] kind  Channel label used in the failure report
      !> @param[in] a     Analytic derivative
      !> @param[in] n     Numeric derivative
      !> @param[in] t     Tolerance applied to this entry
      !> @param[in] iat   Displaced atom
      !> @param[in] idir  Displacement axis
      !> @param[in] ig    Row gridpoint
      !> @param[in] jg    Column gridpoint
      subroutine record_worst(kind, a, n, t, iat, idir, ig, jg)
         !> Channel label
         character(len=*), intent(in) :: kind
         !> Analytic derivative, numeric derivative, and tolerance
         real(wp), intent(in) :: a, n, t
         !> Displaced atom, displacement axis, and the two gridpoints
         integer, intent(in) :: iat, idir, ig, jg

         !> Deviation in units of the tolerance
         real(wp) :: ratio

         ratio = abs(a - n)/t
         if (ratio <= worst_ratio) return
         worst_ratio = ratio
         worst_kind = kind
         worst_a = a
         worst_n = n
         worst_tol = t
         worst_iat = iat
         worst_idir = idir
         worst_i = ig
         worst_j = jg
      end subroutine record_worst

   end subroutine do_test

   !> Fill per-atom CPCM radii, turning a failed lookup into a test failure.
   subroutine fill_cpcm_radii(mol, radii, error)
      !> Structure whose per-atom radii are filled
      type(structure_type), intent(in) :: mol
      !> Allocated on exit to mol%nat
      real(wp), allocatable, intent(out) :: radii(:)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      type(mctc_error), allocatable :: err
      integer :: iat

      allocate (radii(mol%nat))
      do iat = 1, mol%nat
         radii(iat) = get_radius_func(mol%num(mol%id(iat)), err)
         if (allocated(err)) then
            call test_failed(error, "radius lookup failed: "//trim(err%message))
            return
         end if
      end do
   end subroutine fill_cpcm_radii

end module test_cavity_drop_cpcm
