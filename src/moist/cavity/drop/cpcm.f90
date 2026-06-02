!> DROP CPCM assembly, energy, and gradients.
submodule(moist_cavity_drop) moist_cavity_drop_cpcm
   use moist_model_component_pcm_solvers, only: solve_pcm_iterative, solve_pcm_cholesky
   use moist_utils_prettyprint, only: new_prettyprinter
   use omp_lib, only: omp_get_max_threads, omp_get_thread_num
   implicit none

contains

   !* ================================================================================= *!
   !*                                   CPCM A Matrix                                   *!
   !* ================================================================================= *!

   !> Asemble the CPCM A matrix
   !>
   !> @note: This is a naive (reference) implementation that computes all entries and derivatives; it should
   !>        not be used as it requires O(ngrid^2) memory for Amat0 and O(3*nsph*ngrid^2) for Amat1_rA.
   !>        Instead, use contract_amat1_q1q2_rA to compute gradients without materializing Amat1_rA.
   module subroutine assemble_Amat012_rA(self, Amat0, Amat1_rA, error)
      class(cavity_type_drop), intent(in) :: self
      real(wp), allocatable, intent(out) :: Amat0(:, :)
      real(wp), allocatable, optional, intent(out) :: Amat1_rA(:, :, :, :)
      type(error_type), allocatable, intent(out) :: error
      ! real(wp), allocatable, optional, intent(out) :: Amat2_rArB(:, :, :, :)

      integer :: igrid, jgrid, iatom, iaxis
      real(wp) :: xi_ij0, r_ij, xi_i0, xi_j, xi2_sum
      real(wp) :: xi_i1_rA(3, self%nsph), xi_j1_rA(3, self%nsph)
      real(wp) :: xi_ij1_rA(3, self%nsph)
      real(wp) :: r_ij_vec(3), r_ij_unit(3), r_ij1_rA(3, self%nsph)
      real(wp) :: exp_term, T, boys_F_01(0:1), F0, F1, pref, pref2

      logical :: do_grad, do_hess

      do_grad = .false.
      if (present(Amat1_rA)) do_grad = .true.
      ! do_hess = .false.
      ! if (present(Amat2_rArB)) do_hess = .true. ; do_grad = .true.

      ! Allocate A matrix
      allocate (Amat0(self%ngrid, self%ngrid), source=0.0_wp)

      if (do_grad) then
         allocate (Amat1_rA(3, self%nsph, self%ngrid, self%ngrid), source=0.0_wp)
      end if

      pref = 2.0_wp/sqrt(pi)
      pref2 = sqrt(2.0_wp/pi)

      !$omp parallel do default(none) &
      !$omp shared(self, Amat0, Amat1_rA, do_grad, pref, pref2) &
      !$omp private(igrid, jgrid, iatom, iaxis, xi_i0, xi_j, xi2_sum, xi_ij0, r_ij, &
      !$omp         xi_i1_rA, xi_j1_rA, xi_ij1_rA, r_ij_vec, r_ij_unit, r_ij1_rA, &
      !$omp         exp_term, T, boys_F_01, F0, F1) &
      !$omp schedule(static)
      do igrid = 1, self%ngrid

         xi_i0 = self%xi0(igrid)
         if (do_grad) then
            xi_i1_rA = self%xi1_rA(:, :, igrid)
         end if

         Amat0(igrid, igrid) = xi_i0*pref2/self%f(igrid)

         if (do_grad) then
            ! \gradrA A_{ii} =
            !   \sqrt{\dfrac{2}{\pi}}
            !   \left[
            !         + \dfrac{1}{f_i}\left(\gradrA\xi_i\right)
            !         -\xi_i\left(\dfrac{1}{f_i^2}\right)\left(\gradrA f_i\right)
            !   \right]
            do iatom = 1, self%nsph
               do iaxis = 1, 3
                  Amat1_rA(iaxis, iatom, igrid, igrid) = pref2*( &
                                                         (1.0_wp/self%f(igrid))*self%xi1_rA(iaxis, iatom, igrid) &
                                                         - xi_i0*(1.0_wp/self%f(igrid)**2)*self%f1_rA(iaxis, iatom, igrid) &
                                                         )
               end do
            end do
         end if ! do_grad

         do jgrid = 1, igrid
            ! Skip diagonal (already handled above)
            if (jgrid == igrid) cycle

            xi_j = self%xi0(jgrid)
            xi2_sum = xi_i0**2 + xi_j**2

            ! Off-diagonal \xi_{ij}=\frac{\xi_i\xi_j}{\sqrt{\xi_i^2+\xi_j^2}}
            xi_ij0 = xi_i0*xi_j/sqrt(xi2_sum)

            if (do_grad) then
               ! \gradrA \xi_{ij} =
               ! + \frac{
               !         + \xi_j\left(\gradrA \xi_i\right)
               !         + \xi_i\left(\gradrA \xi_j\right)
               !         }{\sqrt{\xi_i^2 + \xi_j^2}}
               ! - \xi_{ij}\frac{
               !                 + \xi_i \left(\gradrA \xi_i\right)
               !                 + \xi_j \left(\gradrA \xi_j\right)
               !                 }{\xi_i^2 + \xi_j^2}
               xi_j1_rA = self%xi1_rA(:, :, jgrid)

               do iatom = 1, self%nsph
                  do iaxis = 1, 3
                     xi_ij1_rA(iaxis, iatom) = &
                        (xi_j*xi_i1_rA(iaxis, iatom) + xi_i0*xi_j1_rA(iaxis, iatom)) &
                        /sqrt(xi2_sum) - xi_ij0*(xi_i0*xi_i1_rA(iaxis, iatom) &
                                                 + xi_j*xi_j1_rA(iaxis, iatom))/xi2_sum
                  end do
               end do
            end if ! do_grad

            ! Compute r_ij vector (needed for both value and gradient)
            r_ij_vec(:) = self%xyz(:, igrid) - self%xyz(:, jgrid)

            ! Distance r_ij (inter-gridpoint) - use shared r_ij_vec for consistency
            r_ij = norm2(r_ij_vec)

            ! A_{ij}=\dfrac{\mathrm{erf}\left(\xi_{ij}r_{ij}\right)}{r_{ij}}
            T = (xi_ij0*r_ij)**2
            call dboysfun1(T, boys_F_01)
            F0 = boys_F_01(0)
            F1 = boys_F_01(1)
            Amat0(igrid, jgrid) = pref*xi_ij0*F0
            ! Amat0(igrid, jgrid) = erf(xi_ij * r_ij) / r_ij
            ! Symmetric matrix entry
            Amat0(jgrid, igrid) = Amat0(igrid, jgrid)

            if (do_grad) then
               ! \gradrA \A_{ij} =
               ! \frac{2}{\sqrt{\pi}}\exp\left(-\xi_{ij}^2 r_{ij}^2\right)\left(\gradrA \xi_{ij}\right)
               !  +\left[
               !        + \frac{2}{\sqrt{\pi}}\frac{\xi_{ij}}{r_{ij}}
               !           \exp\left(-\xi_{ij}^2 r_{ij}^2\right)
               !        - \frac{\operatorname{erf}\!\left(\xi_{ij} r_{ij}\right)}{r_{ij}^2}
               !   \right]\left(\gradrA r_{ij}\right)

               ! Compute r_ij derivative using reformulated expression:
               ! Since second term has r_ij x rhat_ij = r_ij x (r_vec/r_ij) = r_vec,
               ! we can work directly with r_vec to avoid division:
               ! \gradrA r_{ij} x r_{ij} = \mathbf{r}_{ij} \cdot (\gradrA \mathbf{r}_i - \gradrA \mathbf{r}_j)
               do iatom = 1, self%nsph
                  do iaxis = 1, 3
                     r_ij1_rA(iaxis, iatom) = dot_product(r_ij_vec, &
                                                          self%xyz1_rA(:, iaxis, iatom, igrid) &
                                                          - self%xyz1_rA(:, iaxis, iatom, jgrid))
                  end do
               end do

               ! Precompute exponential term (uses same T as value for consistency)
               exp_term = exp(-T)
               ! exp_term = exp(-xi_ij0**2 * r_ij**2)

               ! Compute A_ij derivative (uses exp_term and Boys F1)
               ! Note: r_ij1_rA contains r_ij x (dr_ij/dR_A), so no r_ij factor needed
               do iatom = 1, self%nsph
                  do iaxis = 1, 3
                     Amat1_rA(iaxis, iatom, igrid, jgrid) = pref*exp_term &
                                                            *xi_ij1_rA(iaxis, iatom) - 2.0_wp*pref*(xi_ij0**3) &
                                                            *F1*r_ij1_rA(iaxis, iatom)

                     ! Symmetric entry
                     Amat1_rA(iaxis, iatom, jgrid, igrid) = Amat1_rA(iaxis, iatom, igrid, jgrid)
                  end do
               end do
            end if ! do_grad

         end do ! jgrid
      end do ! igrid
      !$omp end parallel do

   end subroutine assemble_Amat012_rA

   !> Wrapper around assemble_Amat012_rA for the cavity get_amat interface
   !> Assembles only the A-matrix (no derivatives) and copies into the
   !> caller-provided array.
   module subroutine get_amat_drop(self, amat, error)
      class(cavity_type_drop), intent(in) :: self
      real(wp), intent(out) :: amat(:, :)
      type(error_type), allocatable, intent(out) :: error

      real(wp), allocatable :: amat_local(:, :)
      type(error_type), allocatable :: local_error

      ! Check dimensions
      if (size(amat, 1) /= self%ngrid .or. size(amat, 2) /= self%ngrid) then
         call fatal_error(error, &
            & "[get_amat_drop] Matrix dimension mismatch")
         return
      end if

      ! Assemble via the DROP Gaussian-based routine (value only, no gradients)
      call self%Amat012_rA(amat_local, error=local_error)
      if (allocated(local_error)) then
         error = local_error
         return
      end if

      ! Copy into caller-provided array
      amat(:, :) = amat_local(:, :)

   end subroutine get_amat_drop

   !> Contract first derivatives of CPCM A with two grid vectors
   !>
   !> Computes:
   !> - `grad_rA = \sum_{ij} q1_i (\partial A_{ij}/\partial R_A) q2_j`
   !>
   !> This routine avoids materializing `Amat1_rA(3, nsph, ngrid, ngrid)` by
   !> contracting on-the-fly.
   module subroutine contract_amat1_q1q2_rA(self, q1, q2, grad_rA, error)
      class(cavity_type_drop), intent(in) :: self
      real(wp), intent(in) :: q1(:)
      real(wp), intent(in) :: q2(:)
      real(wp), intent(out) :: grad_rA(3, self%nsph)
      type(error_type), allocatable, intent(out) :: error

      integer :: igrid, iatom, iaxis
      real(wp), allocatable :: w_xi(:), w_f(:), b_dr(:, :)
      ! Explicit thread-local buffers: gfortran's OMP reduction on array dummy
      ! arguments (grad_rA(3, self%nsph)) was observed to race intermittently.
      ! Each thread writes its own buffer slice; serial final sum is deterministic.
      real(wp), allocatable :: thread_grad(:, :, :)
      integer :: tid, nthreads_max

      real(wp), parameter :: q_weight_tol = 1.0e-30_wp

      if (.not. allocated(self%xi1_rA) .or. .not. allocated(self%f1_rA) &
          .or. .not. allocated(self%xyz1_rA)) then
         call fatal_error(error, &
                          "contract_amat1_q1q2_rA: missing xi1_rA/f1_rA/xyz1_rA. Did you call get_gradient?")
         return
      end if

      grad_rA = 0.0_wp

      ! TODO: move allocation into omp region for memory optimization
      nthreads_max = omp_get_max_threads()
      allocate (thread_grad(3, self%nsph, nthreads_max), source=0.0_wp)
      allocate (w_xi(self%ngrid), source=0.0_wp)
      allocate (w_f(self%ngrid), source=0.0_wp)
      allocate (b_dr(3, self%ngrid), source=0.0_wp)

      call self%contract_amat1_q1q2_surface_weights(q1, q2, w_xi, w_f, b_dr, error)
      if (allocated(error)) return

      !$omp parallel do default(none) &
      !$omp shared(self, w_xi, w_f, b_dr, thread_grad) &
      !$omp private(igrid, iatom, iaxis, tid) schedule(static)
      do igrid = 1, self%ngrid
         if (abs(w_xi(igrid)) <= q_weight_tol .and. abs(w_f(igrid)) <= q_weight_tol &
             .and. abs(b_dr(1, igrid)) <= q_weight_tol &
             .and. abs(b_dr(2, igrid)) <= q_weight_tol &
             .and. abs(b_dr(3, igrid)) <= q_weight_tol) cycle

         tid = omp_get_thread_num() + 1
         do iatom = 1, self%nsph
            do iaxis = 1, 3
               thread_grad(iaxis, iatom, tid) = thread_grad(iaxis, iatom, tid) &
                                                + w_xi(igrid)*self%xi1_rA(iaxis, iatom, igrid) &
                                                + w_f(igrid)*self%f1_rA(iaxis, iatom, igrid) &
                                                + dot_product(self%xyz1_rA(:, iaxis, iatom, igrid), b_dr(:, igrid))
            end do
         end do
      end do
      !$omp end parallel do

      ! FIXME: Deterministic serial reduction over thread-local buffers
      do tid = 1, nthreads_max
         grad_rA = grad_rA + thread_grad(:, :, tid)
      end do

   end subroutine contract_amat1_q1q2_rA

   !> Contract first derivatives of CPCM A to per-grid surface weights
   module subroutine contract_amat1_q1q2_surface_weights(self, q1, q2, w_xi, w_f, w_xyz, error)
      class(cavity_type_drop), intent(in) :: self
      real(wp), intent(in) :: q1(:)
      real(wp), intent(in) :: q2(:)
      real(wp), intent(out) :: w_xi(:)
      real(wp), intent(out) :: w_f(:)
      real(wp), intent(out) :: w_xyz(:, :)
      type(error_type), allocatable, intent(out) :: error

      integer :: igrid, jgrid
      real(wp) :: q1_i, q2_i, qdiag, qsym
      real(wp) :: xi_i0, xi_j0, xi_ij0, xi2_sum, inv_xi2_sum, inv_sqrt_xi2_sum
      real(wp) :: dxi_coeff_i, dxi_coeff_j
      real(wp) :: inv_f_i, inv_f_i2
      real(wp) :: r_ij_vec(3), r_ij2
      real(wp) :: T, boys_F_01(0:1), F1, exp_term
      real(wp) :: pref, pref2, term_dxi, term_dr

      real(wp), parameter :: q_weight_tol = 1.0e-30_wp

      if (size(q1) /= self%ngrid .or. size(q2) /= self%ngrid) then
         call fatal_error(error, &
                          "contract_amat1_q1q2_surface_weights: q1/q2 size mismatch with ngrid")
         return
      end if

      if (size(w_xi) /= self%ngrid .or. size(w_f) /= self%ngrid &
          .or. size(w_xyz, 1) /= 3 .or. size(w_xyz, 2) /= self%ngrid) then
         call fatal_error(error, &
                          "contract_amat1_q1q2_surface_weights: output size mismatch with cavity")
         return
      end if

      if (.not. allocated(self%xi0) .or. .not. allocated(self%f) .or. .not. allocated(self%xyz)) then
         call fatal_error(error, &
                          "contract_amat1_q1q2_surface_weights: cavity is missing xyz/xi/f data")
         return
      end if

      w_xi = 0.0_wp
      w_f = 0.0_wp
      w_xyz = 0.0_wp

      pref = 2.0_wp/sqrt(pi)
      pref2 = sqrt(2.0_wp/pi)

      ! Reformulate into scalar/vector weights per grid point:
      ! - w_xi(i), w_f(i) multiply dxi_i and df_i
      ! - w_xyz(:, i) multiplies dxyz_i
      do igrid = 1, self%ngrid
         q1_i = q1(igrid)
         q2_i = q2(igrid)

         ! If both i-side weights are tiny, all terms involving gridpoint i vanish.
         if (abs(q1_i) <= q_weight_tol .and. abs(q2_i) <= q_weight_tol) cycle

         xi_i0 = self%xi0(igrid)

         ! Diagonal contraction: q1_i * dA_ii * q2_i
         qdiag = q1_i*q2_i
         if (abs(qdiag) > q_weight_tol) then
            inv_f_i = 1.0_wp/self%f(igrid)
            inv_f_i2 = inv_f_i*inv_f_i
            w_xi(igrid) = w_xi(igrid) + qdiag*pref2*inv_f_i
            w_f(igrid) = w_f(igrid) - qdiag*pref2*xi_i0*inv_f_i2
         end if

         ! Off-diagonal contraction over unique pairs with symmetric q-weight.
         do jgrid = 1, igrid - 1
            qsym = q1_i*q2(jgrid) + q1(jgrid)*q2_i
            if (abs(qsym) <= q_weight_tol) cycle

            xi_j0 = self%xi0(jgrid)
            xi2_sum = xi_i0*xi_i0 + xi_j0*xi_j0
            inv_xi2_sum = 1.0_wp/xi2_sum
            inv_sqrt_xi2_sum = sqrt(inv_xi2_sum)

            xi_ij0 = xi_i0*xi_j0*inv_sqrt_xi2_sum
            dxi_coeff_i = xi_j0*inv_sqrt_xi2_sum - xi_ij0*xi_i0*inv_xi2_sum
            dxi_coeff_j = xi_i0*inv_sqrt_xi2_sum - xi_ij0*xi_j0*inv_xi2_sum

            r_ij_vec(:) = self%xyz(:, igrid) - self%xyz(:, jgrid)
            r_ij2 = r_ij_vec(1)*r_ij_vec(1) &
                    + r_ij_vec(2)*r_ij_vec(2) &
                    + r_ij_vec(3)*r_ij_vec(3)

            T = (xi_ij0*xi_ij0)*r_ij2
            call dboysfun1(T, boys_F_01)
            F1 = boys_F_01(1)
            exp_term = exp(-T)

            term_dxi = qsym*pref*exp_term
            term_dr = -2.0_wp*qsym*pref*(xi_ij0**3)*F1

            w_xi(igrid) = w_xi(igrid) + term_dxi*dxi_coeff_i
            w_xi(jgrid) = w_xi(jgrid) + term_dxi*dxi_coeff_j

            w_xyz(:, igrid) = w_xyz(:, igrid) + term_dr*r_ij_vec(:)
            w_xyz(:, jgrid) = w_xyz(:, jgrid) - term_dr*r_ij_vec(:)
         end do
      end do

   end subroutine contract_amat1_q1q2_surface_weights

   !> Contract combined nuclear and electronic contributions for CPCM gradients
   !>
   !> This routine computes both terms in one grid-point loop:
   !> - direct nuclear contribution from `\delta_{AK}` terms
   !> - chain-rule contribution through `\partial r_i / \partial R_A`
   module subroutine contract_nuc_elec_qefield_rA(self, surface_q, qefield, za, grad_rA, error)
      class(cavity_type_drop), intent(in) :: self
      real(wp), intent(in) :: surface_q(:)
      real(wp), intent(in) :: qefield(:, :)
      real(wp), intent(in) :: za(:)
      real(wp), intent(out) :: grad_rA(3, self%nsph)
      type(error_type), allocatable, intent(out) :: error

      integer :: igrid, iatom, katom, nsource
      real(wp) :: qi, zk, r2, inv_r3
      real(wp) :: r_iK_vec(3), e_nuc(3), chain_vec(3), direct_scale
      real(wp), parameter :: r2_tol = 1.0e-30_wp
      real(wp), parameter :: q_weight_tol = 1.0e-30_wp
      real(wp), allocatable :: thread_grad(:, :, :)
      integer :: tid, nthreads_max

      if (size(surface_q) /= self%ngrid) then
         call fatal_error(error, &
                          "contract_nuc_elec_qefield_rA: surface_q size mismatch with ngrid")
         return
      end if

      if (size(qefield, 1) /= 3 .or. size(qefield, 2) /= self%ngrid) then
         call fatal_error(error, &
                          "contract_nuc_elec_qefield_rA: qefield must have shape (3, ngrid)")
         return
      end if

      if (size(za) /= self%nsph) then
         call fatal_error(error, &
                          "contract_nuc_elec_qefield_rA: za size mismatch with nsph")
         return
      end if

      if (.not. allocated(self%xyz1_rA)) then
         call fatal_error(error, &
                          "contract_nuc_elec_qefield_rA: xyz1_rA is not allocated. Did you call get_gradient?")
         return
      end if

      grad_rA = 0.0_wp
      nsource = min(self%nsph, self%mol%nat)
      ! TODO: move allocation into omp region for memory optimization
      nthreads_max = omp_get_max_threads()
      allocate (thread_grad(3, self%nsph, nthreads_max), source=0.0_wp)

      !$omp parallel do default(none) &
      !$omp shared(self, surface_q, qefield, za, nsource, thread_grad) &
      !$omp private(igrid, iatom, katom, qi, zk, r2, inv_r3, r_iK_vec, e_nuc, chain_vec, &
      !$omp         direct_scale, tid) &
      !$omp schedule(static)
      do igrid = 1, self%ngrid
         qi = surface_q(igrid)
         if (abs(qi) <= q_weight_tol &
             .and. abs(qefield(1, igrid)) <= q_weight_tol &
             .and. abs(qefield(2, igrid)) <= q_weight_tol &
             .and. abs(qefield(3, igrid)) <= q_weight_tol) cycle

         tid = omp_get_thread_num() + 1
         e_nuc = 0.0_wp

         ! Build nuclear electric field at grid point i and direct nuclear term
         do katom = 1, nsource
            r_iK_vec(:) = self%xyz(:, igrid) - self%mol%xyz(:, katom)
            r2 = r_iK_vec(1)*r_iK_vec(1) &
                 + r_iK_vec(2)*r_iK_vec(2) &
                 + r_iK_vec(3)*r_iK_vec(3)
            if (r2 <= r2_tol) cycle

            inv_r3 = 1.0_wp/(sqrt(r2)*r2)
            zk = za(katom)

            e_nuc(1) = e_nuc(1) + zk*inv_r3*r_iK_vec(1)
            e_nuc(2) = e_nuc(2) + zk*inv_r3*r_iK_vec(2)
            e_nuc(3) = e_nuc(3) + zk*inv_r3*r_iK_vec(3)

            direct_scale = qi*zk*inv_r3
            thread_grad(1, katom, tid) = thread_grad(1, katom, tid) + direct_scale*r_iK_vec(1)
            thread_grad(2, katom, tid) = thread_grad(2, katom, tid) + direct_scale*r_iK_vec(2)
            thread_grad(3, katom, tid) = thread_grad(3, katom, tid) + direct_scale*r_iK_vec(3)
         end do

         ! Chain-rule vector shared by all perturbed atoms A
         chain_vec(1) = qefield(1, igrid) - qi*e_nuc(1)
         chain_vec(2) = qefield(2, igrid) - qi*e_nuc(2)
         chain_vec(3) = qefield(3, igrid) - qi*e_nuc(3)

         do iatom = 1, self%nsph
            thread_grad(1, iatom, tid) = thread_grad(1, iatom, tid) &
                                         + dot_product(self%xyz1_rA(:, 1, iatom, igrid), chain_vec)

            thread_grad(2, iatom, tid) = thread_grad(2, iatom, tid) &
                                         + dot_product(self%xyz1_rA(:, 2, iatom, igrid), chain_vec)

            thread_grad(3, iatom, tid) = thread_grad(3, iatom, tid) &
                                         + dot_product(self%xyz1_rA(:, 3, iatom, igrid), chain_vec)
         end do
      end do
      !$omp end parallel do

      ! FIXME: Deterministic serial reduction over thread-local buffers
      do tid = 1, nthreads_max
         grad_rA = grad_rA + thread_grad(:, :, tid)
      end do

   end subroutine contract_nuc_elec_qefield_rA

   !* ================================================================================= *!
   !*                         Fixed charge routines for testing                         *!
   !* ================================================================================= *!

   !> Compute CPCM solvation energy using fixed partial charges
   module subroutine compute_cpcm_energy(self, error)
      class(cavity_type_drop), intent(inout) :: self
      type(error_type), allocatable, intent(out) :: error

      real(wp), allocatable :: Amat(:, :), V(:), q(:)
      type(prettyprinter) :: pp
      integer :: iat, igrid
      real(wp) :: eps_r, r_ij
      real(wp), parameter :: au_to_angstrom = 0.529177_wp

      ! Dielectric constant (eps_r = 100 for now)
      eps_r = 100.0_wp

      ! Allocate arrays
      allocate (Amat(self%ngrid, self%ngrid), source=0.0_wp)
      allocate (V(self%ngrid), source=0.0_wp)
      allocate (q(self%ngrid), source=0.0_wp)

      ! Allocate debug storage
      if (allocated(self%cpcm_q)) deallocate (self%cpcm_q)
      if (allocated(self%cpcm_pot)) deallocate (self%cpcm_pot)
      allocate (self%cpcm_q(self%ngrid))
      allocate (self%cpcm_pot(self%ngrid))
      if (allocated(self%cpcm_source_charges)) deallocate (self%cpcm_source_charges)
      allocate (self%cpcm_source_charges(self%nsph))

      ! Assign simple alternating charges that sum to +1.0 e for testing
      do iat = 1, self%nsph
         if (mod(iat, 2) == 1) then
            self%cpcm_source_charges(iat) = 0.7_wp
         else
            self%cpcm_source_charges(iat) = -0.7_wp
         end if
      end do
      self%cpcm_source_charges(self%nsph) = self%cpcm_source_charges(self%nsph) &
                                            - sum(self%cpcm_source_charges) + 1.0_wp
      ! Assemble A matrix (CPCM A-matrix)
      call self%Amat012_rA(Amat, error=error)
      if (allocated(error)) return

      ! Compute electrostatic potential V at each grid point from atomic charges
      ! V_i = \sum_{atoms} q_j / r_ij
      !$omp parallel do default(shared) private(igrid, iat, r_ij)
      do igrid = 1, self%ngrid
         V(igrid) = 0.0_wp
         do iat = 1, self%mol%nat
            r_ij = norm2(self%xyz(:, igrid) - self%mol%xyz(:, iat))
            V(igrid) = V(igrid) + self%cpcm_source_charges(iat)/r_ij
         end do
      end do
      !$omp end parallel do

      ! Store potential for debugging
      self%cpcm_pot = V

      ! Solve: A*q = -f(eps)*V for surface charges q using Cholesky solver
      call solve_pcm_cholesky(Amat, -((eps_r - 1.0_wp)/eps_r)*V, q, error=error)
      if (allocated(error)) return

      ! Store surface charges for debugging
      self%cpcm_q = q

      ! Compute CPCM solvation energy: E = 1/2 q^T . V
      self%cpcm_energy = 0.5_wp*dot_product(q, V)

      ! Print debug output
      if (self%verbosity > 1) then
         pp = new_prettyprinter(unit=output_unit, fmt_len=20)
         call pp%blank()
         call pp%push('=== CPCM Debug Output [TEMPORARY] ===')
         call pp%kv('Dielectric constant', eps_r)
         call pp%kv('Sum of atomic charge', sum(self%cpcm_source_charges), 'e')
         call pp%kv('Mean atomic charge', sum(self%cpcm_source_charges)/real(self%mol%nat, wp), 'e')
         call pp%kv('Abs. mean atomic charge', sum(abs(self%cpcm_source_charges))/real(self%mol%nat, wp), 'e')
         call pp%kv('CPCM solvation energy', self%cpcm_energy, 'Hartree')
         call pp%kv('CPCM solvation energy', self%cpcm_energy*627.509_wp, 'kcal/mol')
         call pp%pop()
         call pp%blank()
      end if

      deallocate (Amat, V, q)

   end subroutine compute_cpcm_energy

   !> Compute CPCM gradient using fixed partial charges
   module subroutine compute_cpcm_energy_gradient(self, error)
      class(cavity_type_drop), intent(inout) :: self
      type(error_type), allocatable, intent(out) :: error

      real(wp), allocatable :: grad_amat(:, :), grad_nuc_elec(:, :)
      real(wp), allocatable :: qefield(:, :)
      integer :: iatom
      real(wp) :: eps_r, eps_factor

      ! Dielectric constant (must match compute_cpcm_energy)
      eps_r = 100.0_wp
      eps_factor = (eps_r - 1.0_wp)/eps_r

      ! Allocate gradient storage
      if (allocated(self%cpcm_gradient)) deallocate (self%cpcm_gradient)
      allocate (self%cpcm_gradient(3, self%nsph), source=0.0_wp)

      ! Need surface charges from energy calculation
      if (.not. allocated(self%cpcm_q)) then
         call fatal_error(error, &
                          "compute_cpcm_energy_gradient: cpcm_q not allocated. Call compute_cpcm_energy first.")
         return
      end if
      if (.not. allocated(self%cpcm_source_charges)) then
         call fatal_error(error, &
                          "compute_cpcm_energy_gradient: cpcm_source_charges not allocated. Call compute_cpcm_energy first.")
         return
      end if
      if (size(self%cpcm_source_charges) /= self%nsph) then
         call fatal_error(error, &
                          "compute_cpcm_energy_gradient: cpcm_source_charges size mismatch with nsph.")
         return
      end if

      ! Contract A-matrix derivative term without materializing Amat1_rA:
      ! grad_A = q^T * (dA/dR) * q
      allocate (grad_amat(3, self%nsph), source=0.0_wp)
      call self%contract_amat1_q1q2_rA(self%cpcm_q, self%cpcm_q, grad_amat, error=error)
      if (allocated(error)) return

      ! A-matrix gradient term: +1/(2*f(eps)) * q^T * (dA/dR) * q
      self%cpcm_gradient(:, :) = 0.5_wp/eps_factor*grad_amat(:, :)

      ! Add potential-derivative term via fused nuclear/electronic contraction.
      ! In the fixed-charge test model there is no electronic field term.
      allocate (qefield(3, self%ngrid), source=0.0_wp)
      allocate (grad_nuc_elec(3, self%nsph), source=0.0_wp)

      call self%contract_nuc_elec_qefield_rA(self%cpcm_q, qefield, self%cpcm_source_charges, &
                                             grad_nuc_elec, error=error)
      if (allocated(error)) return

      self%cpcm_gradient(:, :) = self%cpcm_gradient(:, :) + grad_nuc_elec(:, :)

      ! Debug output
      if (self%verbosity > 1) then
         write (output_unit, '(/,A)') "=== CPCM Gradient Debug ==="
         write (output_unit, '(A,F20.12,A)') "Gradient norm: ", &
            norm2(self%cpcm_gradient), " Hartree/Bohr"
         write (output_unit, '(A,F20.12,A)') "Max gradient:  ", &
            maxval(abs(self%cpcm_gradient)), " Hartree/Bohr"
         write (output_unit, '(A,/)') "==========================="
      end if

      ! Machine-readable gradient output for MD code
      if (self%verbosity > 1) then
         write (output_unit, '(/,A)') "CPCM_GRADIENT_START"
         do iatom = 1, self%nsph
            write (output_unit, '(3(ES22.15,1X))') self%cpcm_gradient(1:3, iatom)
         end do
         write (output_unit, '(A,/)') "CPCM_GRADIENT_END"
      end if

   end subroutine compute_cpcm_energy_gradient

end submodule moist_cavity_drop_cpcm
