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
   use, intrinsic :: iso_fortran_env, only: error_unit
   implicit none
   private

   public :: collect_cavity_drop_cpcm

   integer, parameter :: ndim = 3

   real(wp), parameter :: k = 1.5
   real(wp), parameter :: gamma = 1.0
   integer, parameter :: NUM_LEB = 26

   real(wp), parameter :: STEP_SIZE = 3.0E-4_wp

   real(wp), parameter :: ATHR = 1.0E-11_wp
   real(wp), parameter :: RTHR = 1.0E-10_wp

   real(wp), parameter :: ATHR_OFFDIAG = 5.0E-9_wp
   real(wp), parameter :: RTHR_OFFDIAG = 5.0E-8_wp

   real(wp), parameter :: ATHR_DIAG = 5.0E-6_wp
   real(wp), parameter :: RTHR_DIAG = 1.0E-5_wp

   real(wp), parameter :: PROJ_TOL = 1E-14_wp
   integer, parameter :: PROJ_MAXITER = 150
   integer, parameter :: PROJ_LEVEL = 2

contains

   subroutine collect_cavity_drop_cpcm(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("contract_amat1_q1q2", test_contract_amat1_q1q2_rA), &
                  new_unittest("contract_nuc_elec_pointcharge_fd", test_contract_nuc_elec_pointcharge_fd), &
                  ! new_unittest("single_atom", test_single_atom), &
                  ! new_unittest("dimer", test_dimer), &
                  ! new_unittest("ar5_blendk_09", test_ar5_blendk_09), &
                  ! new_unittest("mb16-43_h2", test_mb16_43_h2), &
                  ! new_unittest("heavy28_h2o", test_heavy28_h2o), &
                  ! new_unittest("heavy28_pbh4", test_heavy28_pbh4), &
                  ! new_unittest("mb16-43_01", test_mb16_43_01) &
                  ! new_unittest("mb16-43_19", test_mb16_43_19), &
                  ! new_unittest("but14diol_1", test_but14diol_1), &
                  new_unittest("but14diol_32", test_but14diol_32), &
                  new_unittest("il16_008", test_il16_008) &
                  ]
   end subroutine collect_cavity_drop_cpcm

   !> Test A matrix gradient for a single atom
   subroutine test_single_atom(error)
      type(error_type), allocatable, intent(out) :: error
      type(structure_type) :: mol
      real(wp), allocatable :: radii(:)
      integer :: iat

      ! Create single oxygen atom
      call new(mol, [8], reshape([0.0_wp, 0.0_wp, 0.0_wp], [3, 1]))

      allocate (radii(mol%nat))
      do iat = 1, mol%nat
         radii(iat) = get_radius_func(mol%num(mol%id(iat)))
      end do

      call do_test(error, mol, radii)
   end subroutine test_single_atom

   !> Test A matrix gradient for a dimer
   subroutine test_dimer(error)
      type(error_type), allocatable, intent(out) :: error
      type(structure_type) :: mol
      real(wp), allocatable :: radii(:)
      integer :: iat

      ! Create dimer (two oxygen atoms)
      call new(mol, [8, 8], reshape([0.0_wp, 0.0_wp, 0.0_wp, &
                                     3.0_wp, 0.0_wp, 0.0_wp], [3, 2]))

      allocate (radii(mol%nat))
      do iat = 1, mol%nat
         radii(iat) = get_radius_func(mol%num(mol%id(iat)))
      end do

      call do_test(error, mol, radii)
   end subroutine test_dimer

   !> Test A matrix gradient for 5-argon geometry with custom blend-k
   subroutine test_ar5_blendk_09(error)
      type(error_type), allocatable, intent(out) :: error
      type(structure_type) :: mol
      real(wp), allocatable :: radii(:)
      integer :: iat

      call new(mol, [18, 18, 18, 18, 18], reshape([ &
                                                  0.2_wp, 0.0_wp, 5.1_wp, &
                                                  -2.2_wp, -2.2_wp, 0.0_wp, &
                                                  2.2_wp, -2.2_wp, 0.0_wp, &
                                                  -2.2_wp, 2.2_wp, 0.0_wp, &
                                                  2.2_wp, 2.2_wp, 0.0_wp], [3, 5]))

      allocate (radii(mol%nat))
      do iat = 1, mol%nat
         radii(iat) = get_radius_func(mol%num(mol%id(iat)))
      end do

      call do_test(error, mol, radii, blend_k_override=0.8_wp)
   end subroutine test_ar5_blendk_09

   !> Test A matrix gradient for MB16-43 h2
   subroutine test_mb16_43_h2(error)
      type(error_type), allocatable, intent(out) :: error
      type(structure_type) :: mol
      real(wp), allocatable :: radii(:)
      integer :: iat

      call get_structure(mol, "MB16-43", "H2")

      allocate (radii(mol%nat))
      do iat = 1, mol%nat
         radii(iat) = get_radius_func(mol%num(mol%id(iat)))
      end do

      call do_test(error, mol, radii)
   end subroutine test_mb16_43_h2

   !> Test A matrix gradient for bih3_h2o system
   subroutine test_bih3_h2o(error)
      type(error_type), allocatable, intent(out) :: error
      type(structure_type) :: mol
      real(wp), allocatable :: radii(:)
      integer :: iat

      call get_structure(mol, "Heavy28", "bih3_h2o")

      allocate (radii(mol%nat))
      do iat = 1, mol%nat
         radii(iat) = get_radius_func(mol%num(mol%id(iat)))
      end do

      call do_test(error, mol, radii)
   end subroutine test_bih3_h2o

   !> Test A matrix gradient for Heavy28 h2o
   subroutine test_heavy28_h2o(error)
      type(error_type), allocatable, intent(out) :: error
      type(structure_type) :: mol
      real(wp), allocatable :: radii(:)
      integer :: iat

      call get_structure(mol, "Heavy28", "h2o")

      allocate (radii(mol%nat))
      do iat = 1, mol%nat
         radii(iat) = get_radius_func(mol%num(mol%id(iat)))
      end do

      call do_test(error, mol, radii)
   end subroutine test_heavy28_h2o

   !> Test A matrix gradient for Heavy28 pbh4
   subroutine test_heavy28_pbh4(error)
      type(error_type), allocatable, intent(out) :: error
      type(structure_type) :: mol
      real(wp), allocatable :: radii(:)
      integer :: iat

      call get_structure(mol, "Heavy28", "pbh4")

      allocate (radii(mol%nat))
      do iat = 1, mol%nat
         radii(iat) = get_radius_func(mol%num(mol%id(iat)))
      end do

      call do_test(error, mol, radii)
   end subroutine test_heavy28_pbh4

   !> Test A matrix gradient for MB16-43 01
   subroutine test_mb16_43_01(error)
      type(error_type), allocatable, intent(out) :: error
      type(structure_type) :: mol
      real(wp), allocatable :: radii(:)
      integer :: iat

      call get_structure(mol, "MB16-43", "01")

      allocate (radii(mol%nat))
      do iat = 1, mol%nat
         radii(iat) = get_radius_func(mol%num(mol%id(iat)), 'cpcm')
      end do

      call do_test(error, mol, radii)
   end subroutine test_mb16_43_01

   !> Test A matrix gradient for MB16-43 19
   subroutine test_mb16_43_19(error)
      type(error_type), allocatable, intent(out) :: error
      type(structure_type) :: mol
      real(wp), allocatable :: radii(:)
      integer :: iat

      call get_structure(mol, "MB16-43", "19")

      allocate (radii(mol%nat))
      do iat = 1, mol%nat
         radii(iat) = get_radius_func(mol%num(mol%id(iat)))
      end do

      call do_test(error, mol, radii)
   end subroutine test_mb16_43_19

   !> Test A matrix gradient for But14diol 1
   subroutine test_but14diol_1(error)
      type(error_type), allocatable, intent(out) :: error
      type(structure_type) :: mol
      real(wp), allocatable :: radii(:)
      integer :: iat

      call get_structure(mol, "But14diol", "1")

      allocate (radii(mol%nat))
      do iat = 1, mol%nat
         radii(iat) = get_radius_func(mol%num(mol%id(iat)))
      end do

      call do_test(error, mol, radii)
   end subroutine test_but14diol_1

   !> Test A matrix gradient for But14diol 32
   subroutine test_but14diol_32(error)
      type(error_type), allocatable, intent(out) :: error
      type(structure_type) :: mol
      real(wp), allocatable :: radii(:)
      integer :: iat

      call get_structure(mol, "But14diol", "32")

      allocate (radii(mol%nat))
      do iat = 1, mol%nat
         radii(iat) = get_radius_func(mol%num(mol%id(iat)))
      end do

      call do_test(error, mol, radii)
   end subroutine test_but14diol_32

   !> Test A matrix gradient for IL16 008
   subroutine test_il16_008(error)
      type(error_type), allocatable, intent(out) :: error
      type(structure_type) :: mol
      real(wp), allocatable :: radii(:)
      integer :: iat

      call get_structure(mol, "IL16", "008")

      allocate (radii(mol%nat))
      do iat = 1, mol%nat
         radii(iat) = get_radius_func(mol%num(mol%id(iat)))
      end do

      call do_test(error, mol, radii)
   end subroutine test_il16_008

   !> Test contracted A-matrix gradient against explicit tensor contraction.
   subroutine test_contract_amat1_q1q2_rA(error)
      type(error_type), allocatable, intent(out) :: error

      type(structure_type) :: mol
      type(cavity_type_drop), allocatable :: cavity
      real(wp), allocatable :: radii(:)
      real(wp), allocatable :: Amat0(:, :), Amat1_rA(:, :, :, :)
      real(wp), allocatable :: q1(:), q2(:)
      real(wp), allocatable :: grad_ref(:, :), grad_ctr(:, :)
      integer :: iat, iaxis, igrid, jgrid, ngrid
      type(mctc_error), allocatable :: cavity_error

      call get_structure(mol, "MB16-43", "04")

      allocate (radii(mol%nat))
      do iat = 1, mol%nat
         radii(iat) = get_radius_func(mol%num(mol%id(iat)))
      end do

      allocate (cavity)
      block
         type(moist_cavity_drop_lsf_svdw_type) :: svdw_template
         call svdw_template%new(blend_k=k, blend_3b=gamma)
         call new_cavity_drop(cavity, nleb=NUM_LEB, &
                              tolerance=PROJ_TOL, proj_maxiter=PROJ_MAXITER, proj_level=PROJ_LEVEL, &
                              debug=.false., verbose=0, radius_model=default_cpcm_radii(), &
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
      call cavity%Amat012_rA(Amat0, Amat1_rA, error=cavity_error)
      if (allocated(cavity_error)) then
         call test_failed(error, cavity_error%message)
         return
      end if

      ngrid = cavity%ngrid

      allocate (q1(ngrid), q2(ngrid))
      allocate (grad_ref(3, mol%nat), grad_ctr(3, mol%nat))

      do igrid = 1, ngrid
         q1(igrid) = real(igrid, wp)/(real(ngrid, wp) + 1.0_wp)
         if (mod(igrid, 2) == 0) then
            q2(igrid) = -1.0_wp/(real(igrid, wp) + 0.5_wp)
         else
            q2(igrid) = 1.0_wp/(real(igrid, wp) + 0.25_wp)
         end if
      end do

      grad_ref = 0.0_wp
      do iat = 1, mol%nat
         do iaxis = 1, 3
            do igrid = 1, ngrid
               do jgrid = 1, ngrid
                  grad_ref(iaxis, iat) = grad_ref(iaxis, iat) &
                                         + q1(igrid)*Amat1_rA(iaxis, iat, igrid, jgrid)*q2(jgrid)
               end do
            end do
         end do
      end do

      call cavity%contract_amat1_q1q2_rA(q1, q2, grad_ctr, error=cavity_error)
      if (allocated(cavity_error)) then
         call test_failed(error, cavity_error%message)
         return
      end if

      do iat = 1, mol%nat
         do iaxis = 1, 3
            call check(error, &
                       grad_ctr(iaxis, iat), &
                       grad_ref(iaxis, iat), &
                       thr_abs=5.0e-11_wp, thr_rel=5.0e-11_wp, &
                       more="contract_amat1_q1q2_rA mismatch")
            if (allocated(error)) return
         end do
      end do
   end subroutine test_contract_amat1_q1q2_rA

   !> Test fused nuclear/electronic contraction against finite differences.
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

      call get_structure(mol, "MB16-43", "15")

      allocate (radii(mol%nat))
      do iat = 1, mol%nat
         radii(iat) = get_radius_func(mol%num(mol%id(iat)))
      end do

      allocate (cavity)
      block
         type(moist_cavity_drop_lsf_svdw_type) :: svdw_template
         call svdw_template%new(blend_k=k, blend_3b=gamma)
         call new_cavity_drop(cavity, nleb=NUM_LEB, &
                              tolerance=PROJ_TOL, proj_maxiter=PROJ_MAXITER, proj_level=PROJ_LEVEL, &
                              wleb_prune_level=3, &
                              debug=.false., verbose=0, radius_model=default_cpcm_radii(), &
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

      call cavity%contract_nuc_elec_qefield_rA(surface_q, qefield, za, grad_ctr, error=cavity_error)
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
                       more="contract_nuc_elec_qefield_rA point-charge FD mismatch")
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
      integer :: iat, idir, igrid, jgrid, ngrid, num_idn, num_jdn, idim
      integer :: max_numbering, idx_i, idx_j
      real(wp), allocatable :: ref_f(:), ref_xi(:)
      integer, allocatable :: ref_owner(:)
      integer, allocatable :: ref_numbering(:), numbering_to_idx(:)
      real(wp) :: diff, max_diff
      logical, allocatable :: valid_gridpoint_ref(:)
      logical, allocatable :: valid_gridpoint(:, :, :)
      logical, allocatable :: nn_converged(:, :, :), n_converged(:, :, :), &
                              p_converged(:, :, :), pp_converged(:, :, :)
      integer, allocatable :: num_nn(:, :, :), num_n(:, :, :), num_p(:, :, :), num_pp(:, :, :)

      real(wp), allocatable :: Amat_copy(:, :), Amat_inv(:, :), identity_test(:, :)
      integer, allocatable :: ipiv(:)
      integer :: info, ii
      real(wp) :: max_err
      real(wp) :: blend_k_local
      type(mctc_error), allocatable :: cavity_error

      ! CPCM energy gradient FD arrays
      real(wp), allocatable :: en_cpcm_gradient(:, :)
      real(wp), allocatable :: num_cpcm_gradient(:, :)
      real(wp) :: nn_energy, n_energy, p_energy, pp_energy

      !> Initialize cavity with 26-point Lebedev grid
      blend_k_local = k
      if (present(blend_k_override)) blend_k_local = blend_k_override
      allocate (cavity)
      block
         type(moist_cavity_drop_lsf_svdw_type) :: svdw_template
         call svdw_template%new(blend_k=blend_k_local, blend_3b=gamma)
         call new_cavity_drop(cavity, nleb=NUM_LEB, &
                              tolerance=PROJ_TOL, proj_maxiter=PROJ_MAXITER, proj_level=PROJ_LEVEL, &
                              debug=.false., verbose=0, do_cpcm=.true., &
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

      !> Get analytic gradients (also computes CPCM energy gradient)
      call cavity%get_gradient()

      !> Store analytical CPCM energy gradient
      allocate (en_cpcm_gradient(ndim, mol%nat), source=cavity%cpcm_gradient)
      allocate (num_cpcm_gradient(ndim, mol%nat), source=0.0_wp)

      !> Assemble A matrix and its gradient
      call cavity%Amat012_rA(Amat0, Amat1_rA, error=cavity_error)
      if (allocated(cavity_error)) then
         call test_failed(error, "Amat012_rA failed: "//cavity_error%message)
         return
      end if

      if (.not. allocated(Amat1_rA)) then
         error stop "Amat1_rA not allocated after Amat012_rA"
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
         write (error_unit, '(A)') "A * A^{-1} differs from identity. Problematic entries:"
         do igrid = 1, ngrid
            do jgrid = 1, ngrid
               if (igrid == jgrid) then
                  if (abs(identity_test(igrid, jgrid) - 1.0_wp) > 1.0e-10_wp) then
                     write (error_unit, '(A,I4,A,I4,A,ES15.6,A,ES15.6)') &
                        "  Diagonal (", igrid, ",", jgrid, "): got ", &
                        identity_test(igrid, jgrid), ", expected 1.0, err = ", &
                        abs(identity_test(igrid, jgrid) - 1.0_wp)
                  end if
               else
                  if (abs(identity_test(igrid, jgrid)) > 1.0e-10_wp) then
                     write (error_unit, '(A,I4,A,I4,A,ES15.6,A,ES15.6)') &
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
            call cavity%get_gradient()
            call cavity%Amat012_rA(Amat0, error=cavity_error)
            if (allocated(cavity_error)) then
               call test_failed(error, "Amat012_rA failed: "//cavity_error%message)
               return
            end if
            nn_energy = cavity%cpcm_energy
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
            call cavity%get_gradient()
            call cavity%Amat012_rA(Amat0, error=cavity_error)
            if (allocated(cavity_error)) then
               call test_failed(error, "Amat012_rA failed: "//cavity_error%message)
               return
            end if
            n_energy = cavity%cpcm_energy
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
            call cavity%get_gradient()
            call cavity%Amat012_rA(Amat0, error=cavity_error)
            if (allocated(cavity_error)) then
               call test_failed(error, "Amat012_rA failed: "//cavity_error%message)
               return
            end if
            p_energy = cavity%cpcm_energy
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
            call cavity%get_gradient()
            call cavity%Amat012_rA(Amat0, error=cavity_error)
            if (allocated(cavity_error)) then
               call test_failed(error, "Amat012_rA failed: "//cavity_error%message)
               return
            end if
            pp_energy = cavity%cpcm_energy
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

            ! CPCM energy gradient via FD
            num_cpcm_gradient(idir, iat) = &
               (-pp_energy + 8.0_wp*p_energy &
                - 8.0_wp*n_energy + nn_energy) &
               /(12.0_wp*STEP_SIZE)

            ! Central difference formula: f'(x) ~= [-f(x+2h) + 8f(x+h) - 8f(x-h) + f(x-2h)] / (12h)
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

      !> Compare analytic vs numeric A matrix gradients (only valid gridpoint pairs)
      !> Off-diagonal elements
      do iat = 1, mol%nat
         do idir = 1, ndim
            do igrid = 1, ngrid
               ! Skip invalid gridpoints for igrid
               if (.not. valid_gridpoint(idir, iat, igrid)) cycle

               do jgrid = 1, ngrid
                  ! Skip invalid gridpoints for jgrid
                  if (.not. valid_gridpoint(idir, iat, jgrid)) cycle

                  if (jgrid == igrid) cycle  ! Skip diagonal elements here

                  ! Additional check: we need both gridpoints to have existed in all 4 FD steps
                  ! This was already checked during FD computation - if any were missing,
                  ! num_Amat1_rA was left as zero. So check if we actually computed a FD value:
                  if (num_nn(idir, iat, igrid) == 0 .or. num_nn(idir, iat, jgrid) == 0) cycle
                  if (num_n(idir, iat, igrid) == 0 .or. num_n(idir, iat, jgrid) == 0) cycle
                  if (num_p(idir, iat, igrid) == 0 .or. num_p(idir, iat, jgrid) == 0) cycle
                  if (num_pp(idir, iat, igrid) == 0 .or. num_pp(idir, iat, jgrid) == 0) cycle

                  call check(error, &
                             en_Amat1_rA(idir, iat, igrid, jgrid), &
                             num_Amat1_rA(idir, iat, igrid, jgrid), &
                             thr_abs=ATHR_OFFDIAG, thr_rel=RTHR_OFFDIAG, &
                             more="Off-diagonal A matrix gradient error")
                  if (allocated(error)) return
               end do
            end do
         end do
      end do

      !> Compare analytic vs numeric A matrix gradients (only valid gridpoint pairs)
      !> Diagonal elements
      do iat = 1, mol%nat
         do idir = 1, ndim
            do igrid = 1, ngrid
               ! Skip invalid gridpoints for igrid
               if (.not. valid_gridpoint(idir, iat, igrid)) cycle

               jgrid = igrid  ! Test diagonal elements only

               ! Additional check: we need both gridpoints to have existed in all 4 FD steps
               ! This was already checked during FD computation - if any were missing,
               ! num_Amat1_rA was left as zero. So check if we actually computed a FD value:
               if (num_nn(idir, iat, igrid) == 0 .or. num_nn(idir, iat, jgrid) == 0) cycle
               if (num_n(idir, iat, igrid) == 0 .or. num_n(idir, iat, jgrid) == 0) cycle
               if (num_p(idir, iat, igrid) == 0 .or. num_p(idir, iat, jgrid) == 0) cycle
               if (num_pp(idir, iat, igrid) == 0 .or. num_pp(idir, iat, jgrid) == 0) cycle

               ! Cycle if switching function (otherwise num. noise dominates as A_ii ~ 1/f)
               if (ref_f(igrid) < 1.0e-4_wp) cycle

               call check(error, &
                          en_Amat1_rA(idir, iat, igrid, jgrid), &
                          num_Amat1_rA(idir, iat, igrid, jgrid), &
                          thr_abs=ATHR_DIAG, thr_rel=RTHR_DIAG, &
                          more="Diagonal A matrix gradient error")
               if (allocated(error)) return
            end do
         end do
      end do

      !> Compare CPCM energy gradient: analytical vs numerical
      do iat = 1, mol%nat
         do idir = 1, ndim
            call check(error, en_cpcm_gradient(idir, iat), &
                       num_cpcm_gradient(idir, iat), &
                       thr_abs=ATHR, thr_rel=RTHR, &
                       more="CPCM energy gradient mismatch")
            if (allocated(error)) return
         end do
      end do

   end subroutine do_test

end module test_cavity_drop_cpcm
