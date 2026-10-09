!> Tables of a MOZ potential: sites resolved through the coupling, Ng-split
!> kernels, and the adjoints handed back to the terms
!>
!> - Sides resolve into plain sites; the kernels see only those and the fields
!> - Adjoint recipients: every term supplying a field or owning a charge, host-fed or not
!> - Charge adjoint masked to the atoms a term owns; zero-sized when it owns none
submodule(moist_model_moz_potential) moist_model_moz_potential_compute
   use mctc_env, only: fatal_error
   use moist_channels_coupling, only: coupling_view_type, coupling_make_view, coupling_close_view
   use moist_math_packed_sym, only: packed_index, npair_from_ns
   use moist_model_moz_potential_kernel_table_1d, only: evaluate_table_1d, evaluate_table_1d_adjoint
   use moist_model_moz_potential_kernel_table_3d, only: evaluate_table_3d, evaluate_table_3d_adjoint
   implicit none(type, external)

contains

   module subroutine potential_site_pairs(self, pairs, error, solvent)
      class(moz_potential_type), intent(in) :: self
      integer, allocatable, intent(out) :: pairs(:, :)
      type(error_type), allocatable, intent(out) :: error
      class(moz_potential_type), intent(in), optional :: solvent
      integer :: i, j, nu, nv
      call self%require_updated(error)
      if (allocated(error)) return
      if (present(solvent)) then
         call solvent%require_updated(error)
         if (allocated(error)) return
         nu = self%nat
         nv = solvent%nat
         allocate (pairs(2, nu*nv))
         do i = 1, nu
            do j = 1, nv
               pairs(:, (i - 1)*nv + j) = [i, j]
            end do
         end do
      else
         allocate (pairs(2, npair_from_ns(self%nat)))
         do j = 1, self%nat
            do i = 1, j
               pairs(:, packed_index(i, j)) = [i, j]
            end do
         end do
      end if
   end subroutine potential_site_pairs

   module subroutine potential_compute_1d(self, rgrid, kgrid, u_sr, ur_lr, uk_lr, error, solvent, coupling)
      class(moz_potential_type), intent(in) :: self
      type(moist_math_grid_radial_type), intent(in) :: rgrid
      type(moist_math_grid_radial_type), intent(in) :: kgrid
      real(wp), allocatable, intent(out) :: u_sr(:, :)
      real(wp), allocatable, intent(out) :: ur_lr(:, :)
      real(wp), allocatable, intent(out) :: uk_lr(:, :)
      type(error_type), allocatable, intent(out) :: error
      class(moz_potential_type), intent(in), optional :: solvent
      class(coupling_type), intent(inout), target, optional :: coupling
      type(potential_sites_type) :: site_u, site_v
      real(wp), allocatable :: phi_u(:, :)
      integer, allocatable :: pairs(:, :)
      call resolve_sites(self, site_u, site_v, error, solvent, coupling)
      if (allocated(error)) return
      call self%site_pairs(pairs, error, solvent)
      if (allocated(error)) return
      if (site_u%has_field) then
         call field_1d(self, rgrid, phi_u, error, coupling)
         if (allocated(error)) return
         call evaluate_table_1d(rgrid, kgrid, site_u, site_v, pairs, self%mixing_rule, self%split, self%split_alpha, &
            & u_sr, ur_lr, uk_lr, error, phi_u)
      else
         call evaluate_table_1d(rgrid, kgrid, site_u, site_v, pairs, self%mixing_rule, self%split, self%split_alpha, &
            & u_sr, ur_lr, uk_lr, error)
      end if
   end subroutine potential_compute_1d

   module subroutine potential_compute_3d(self, grid, u_sr, ur_lr, uk_lr, error, solvent, coupling)
      class(moz_potential_type), intent(in) :: self
      class(moist_math_grid_3d_type), intent(in) :: grid
      real(wp), allocatable, intent(out) :: u_sr(:, :)
      real(wp), allocatable, intent(out) :: ur_lr(:, :)
      complex(wp), allocatable, intent(out) :: uk_lr(:, :)
      type(error_type), allocatable, intent(out) :: error
      class(moz_potential_type), intent(in), optional :: solvent
      class(coupling_type), intent(inout), target, optional :: coupling
      type(potential_sites_type) :: site_u, site_v
      real(wp), allocatable :: phi_u(:)
      if (.not. present(solvent)) then
         call fatal_error(error, "3D tables are solute-solvent tables; pass the solvent potential")
         return
      end if
      call resolve_sites(self, site_u, site_v, error, solvent, coupling)
      if (allocated(error)) return
      if (site_u%has_field) then
         call field_3d(self, grid, phi_u, error, coupling)
         if (allocated(error)) return
         call evaluate_table_3d(grid, self%mol%xyz, site_u, site_v, self%mixing_rule, self%split, self%split_alpha, &
            & u_sr, ur_lr, uk_lr, error, phi_u)
      else
         call evaluate_table_3d(grid, self%mol%xyz, site_u, site_v, self%mixing_rule, self%split, self%split_alpha, &
            & u_sr, ur_lr, uk_lr, error)
      end if
   end subroutine potential_compute_3d

   module subroutine potential_compute_adjoint_1d(self, rgrid, kgrid, solvent, coupling, w_sr, w_lr, w_k, &
         & response, error)
      class(moz_potential_type), intent(inout) :: self
      type(moist_math_grid_radial_type), intent(in) :: rgrid
      type(moist_math_grid_radial_type), intent(in) :: kgrid
      class(moz_potential_type), intent(in) :: solvent
      class(coupling_type), intent(inout), target :: coupling
      real(wp), intent(in) :: w_sr(:, :)
      real(wp), intent(in) :: w_lr(:, :)
      real(wp), intent(in) :: w_k(:, :)
      type(response_type), intent(inout) :: response
      type(error_type), allocatable, intent(out) :: error
      type(potential_sites_type) :: site_u, site_v
      real(wp), allocatable :: w_phi(:, :), w_q(:)
      integer, allocatable :: pairs(:, :)
      call resolve_sites(self, site_u, site_v, error, solvent, coupling)
      if (allocated(error)) return
      call self%site_pairs(pairs, error, solvent)
      if (allocated(error)) return
      call evaluate_table_1d_adjoint(rgrid, kgrid, site_u, site_v, pairs, self%split, self%split_alpha, &
         & w_sr, w_lr, w_k, w_phi, w_q, error)
      if (allocated(error)) return
      if (.not. allocated(w_phi)) return
      call accumulate_adjoint_1d(self, rgrid, coupling, w_phi, w_q, response, error)
   end subroutine potential_compute_adjoint_1d

   module subroutine potential_compute_adjoint_3d(self, grid, solvent, coupling, w_sr, w_lr, w_k, response, &
         & error, gradient, grid_adjoint)
      class(moz_potential_type), intent(inout) :: self
      class(moist_math_grid_3d_type), intent(in) :: grid
      class(moz_potential_type), intent(in) :: solvent
      class(coupling_type), intent(inout), target :: coupling
      real(wp), intent(in) :: w_sr(:, :)
      real(wp), intent(in) :: w_lr(:, :)
      complex(wp), intent(in) :: w_k(:, :)
      type(response_type), intent(inout) :: response
      type(error_type), allocatable, intent(out) :: error
      real(wp), intent(inout), optional :: gradient(:, :)
      type(volume_adjoint_type), intent(inout), optional :: grid_adjoint
      type(potential_sites_type) :: site_u, site_v
      real(wp), allocatable :: w_phi(:), w_q(:)
      call resolve_sites(self, site_u, site_v, error, solvent, coupling)
      if (allocated(error)) return
      call evaluate_table_3d_adjoint(grid, self%mol%xyz, site_u, site_v, self%mixing_rule, self%split, &
         & self%split_alpha, w_sr, w_lr, w_k, w_phi, w_q, error, gradient, grid_adjoint)
      if (allocated(error)) return
      if (.not. allocated(w_phi)) return
      call accumulate_adjoint_3d(self, grid, coupling, w_phi, w_q, response, error, grid_adjoint, gradient)
   end subroutine potential_compute_adjoint_3d

   !> Resolve the sites of both sides
   !>
   !> - UV: first side through `coupling`; the solvent must read no host data
   !> - VV (no `solvent`): this potential on both sides, which must read no host data
   !>
   !> @param[in]     self      updated potential, first side
   !> @param[out]    site_u    resolved sites of the first side
   !> @param[out]    site_v    resolved sites of the second side
   !> @param[out]    error     not updated, host data on the wrong side or a charge read failure
   !> @param[in]     solvent   optional updated solvent potential, second side
   !> @param[in,out] coupling  optional coupling of the model; host-data term i reads scope i
   subroutine resolve_sites(self, site_u, site_v, error, solvent, coupling)
      !> Updated potential, first side
      class(moz_potential_type), intent(in) :: self
      !> Resolved sites of the first side
      type(potential_sites_type), intent(out) :: site_u
      !> Resolved sites of the second side
      type(potential_sites_type), intent(out) :: site_v
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      !> Optional solvent potential, second side
      class(moz_potential_type), intent(in), optional :: solvent
      !> Optional coupling of the model
      class(coupling_type), intent(inout), target, optional :: coupling
      call self%require_updated(error)
      if (allocated(error)) return
      if (.not. present(solvent)) then
         if (self%reads_host_data()) then
            call fatal_error(error, "VV tables take no host-data terms; pass the solvent potential for UV tables")
            return
         end if
         call self%sites(site_u, error)
         if (allocated(error)) return
         site_v = site_u
         return
      end if
      call solvent%require_updated(error)
      if (allocated(error)) return
      if (solvent%reads_host_data()) then
         call fatal_error(error, "Solvent potential must not read host data")
         return
      end if
      call self%sites(site_u, error, coupling)
      if (allocated(error)) return
      call solvent%sites(site_v, error)
   end subroutine resolve_sites

   !> Sum the radial fields of the field-supplying terms
   !>
   !> - Host terms through scoped view of `coupling`; refusal without coupling
   !>
   !> @param[in]     self      updated potential
   !> @param[in]     grid      radial grid
   !> @param[out]    phi       summed field per atom, Hartree/e (npts, natom); unallocated when no term supplies one
   !> @param[out]    error     field read failure, host-data term without coupling or a field not matching the grid
   !> @param[in,out] coupling  optional coupling of the model; host-data term i reads scope i
   subroutine field_1d(self, grid, phi, error, coupling)
      !> Updated potential
      class(moz_potential_type), intent(in) :: self
      !> Radial grid
      type(moist_math_grid_radial_type), intent(in) :: grid
      !> Summed field per atom (npts, natom)
      real(wp), allocatable, intent(out) :: phi(:, :)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      !> Optional coupling of the model
      class(coupling_type), intent(inout), target, optional :: coupling
      !> Scoped view of a host-data term and the unbound view handed to parameter terms
      type(coupling_view_type) :: view, unbound
      real(wp), allocatable :: phi_term(:, :)
      integer :: i
      do i = 1, self%n_terms()
         associate (term => self%slots(i)%term)
            if (.not. term%has_field()) cycle
            if (term%is_host_fed()) then
               if (.not. present(coupling)) then
                  call fatal_error(error, "Potential term '"//trim(term%name())// &
                     & "' reads host data; its field needs the coupling")
                  return
               end if
               call coupling_make_view(coupling, i, view)
               call term%field(grid, view, phi_term, error)
               call coupling_close_view(coupling)
            else
               call term%field(grid, unbound, phi_term, error)
            end if
            if (allocated(error)) return
            if (any(shape(phi_term) /= [grid%npts, self%nat])) then
               call fatal_error(error, "Potential term '"//trim(term%name())// &
                  & "' radial field is not (npts, natom) of its side")
               return
            end if
            if (allocated(phi)) then
               phi = phi + phi_term
            else
               call move_alloc(phi_term, phi)
            end if
         end associate
      end do
   end subroutine field_1d

   !> Sum the volume-grid fields of the field-supplying terms
   !>
   !> - Host terms through scoped view of `coupling`; refusal without coupling
   !>
   !> @param[in]     self      updated potential
   !> @param[in]     grid      volume grid
   !> @param[out]    phi       summed field, Hartree/e (ngrid); unallocated when no term supplies one
   !> @param[out]    error     field read failure, host-data term without coupling or a field not matching the grid
   !> @param[in,out] coupling  optional coupling of the model; host-data term i reads scope i
   subroutine field_3d(self, grid, phi, error, coupling)
      !> Updated potential
      class(moz_potential_type), intent(in) :: self
      !> Volume grid
      class(moist_math_grid_3d_type), intent(in) :: grid
      !> Summed field (ngrid)
      real(wp), allocatable, intent(out) :: phi(:)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      !> Optional coupling of the model
      class(coupling_type), intent(inout), target, optional :: coupling
      !> Scoped view of a host-data term and the unbound view handed to parameter terms
      type(coupling_view_type) :: view, unbound
      real(wp), allocatable :: phi_term(:)
      integer :: i
      do i = 1, self%n_terms()
         associate (term => self%slots(i)%term)
            if (.not. term%has_field()) cycle
            if (term%is_host_fed()) then
               if (.not. present(coupling)) then
                  call fatal_error(error, "Potential term '"//trim(term%name())// &
                     & "' reads host data; its field needs the coupling")
                  return
               end if
               call coupling_make_view(coupling, i, view)
               call term%field(grid, view, phi_term, error)
               call coupling_close_view(coupling)
            else
               call term%field(grid, unbound, phi_term, error)
            end if
            if (allocated(error)) return
            if (size(phi_term) /= grid%ngrid) then
               call fatal_error(error, "Potential term '"//trim(term%name())// &
                  & "' field does not match the grid")
               return
            end if
            if (allocated(phi)) then
               phi = phi + phi_term
            else
               call move_alloc(phi_term, phi)
            end if
         end associate
      end do
   end subroutine field_3d

   !> Hand the field and charge adjoints to every contributing term (1D)
   !>
   !> @param[in,out] self      updated potential
   !> @param[in]     grid      radial grid
   !> @param[in,out] coupling  coupling of the model; host-data term i reads scope i
   !> @param[in]     w_phi     adjoint of the summed field per atom (npts, natom)
   !> @param[in]     w_q       adjoint of the merged tail charges (natom); zero-sized when the side has none
   !> @param[in,out] response  host response items
   !> @param[out]    error     term adjoint failure
   subroutine accumulate_adjoint_1d(self, grid, coupling, w_phi, w_q, response, error)
      !> Updated potential
      class(moz_potential_type), intent(inout) :: self
      !> Radial grid
      type(moist_math_grid_radial_type), intent(in) :: grid
      !> Coupling of the model
      class(coupling_type), intent(inout), target :: coupling
      !> Adjoint of the summed field per atom (npts, natom)
      real(wp), intent(in) :: w_phi(:, :)
      !> Adjoint of the merged tail charges (natom)
      real(wp), intent(in) :: w_q(:)
      !> Host response items
      type(response_type), intent(inout) :: response
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      !> Scoped view of a host-data term and the unbound view handed to parameter terms
      type(coupling_view_type) :: view, unbound
      real(wp), allocatable :: w_qinf(:)
      logical :: contributes
      integer :: i
      call require_charge_adjoint(self, w_q, error)
      if (allocated(error)) return
      do i = 1, self%n_terms()
         associate (term => self%slots(i)%term)
            call term_charge_adjoint(self, i, w_q, w_qinf, contributes)
            if (.not. contributes) cycle
            if (term%is_host_fed()) then
               call coupling_make_view(coupling, i, view)
               call term%accumulate_adjoint(grid, view, w_phi, w_qinf, response, error)
               call coupling_close_view(coupling)
            else
               call term%accumulate_adjoint(grid, unbound, w_phi, w_qinf, response, error)
            end if
            deallocate (w_qinf)
            if (allocated(error)) return
         end associate
      end do
   end subroutine accumulate_adjoint_1d

   !> Hand the field and charge adjoints to every contributing term (3D)
   !>
   !> @param[in,out] self          updated potential
   !> @param[in]     grid          volume grid
   !> @param[in,out] coupling      coupling of the model; host-data term i reads scope i
   !> @param[in]     w_phi         adjoint of the summed field (ngrid)
   !> @param[in]     w_q           adjoint of the merged tail charges (natom); zero-sized when the side has none
   !> @param[in,out] response      host response items
   !> @param[out]    error         term adjoint failure
   !> @param[in,out] grid_adjoint  optional volume adjoint accumulator; present when the geometry derivative is requested
   !> @param[in,out] gradient      optional explicit nuclear gradient, Hartree/bohr (3, natom)
   subroutine accumulate_adjoint_3d(self, grid, coupling, w_phi, w_q, response, error, grid_adjoint, gradient)
      !> Updated potential
      class(moz_potential_type), intent(inout) :: self
      !> Volume grid
      class(moist_math_grid_3d_type), intent(in) :: grid
      !> Coupling of the model
      class(coupling_type), intent(inout), target :: coupling
      !> Adjoint of the summed field (ngrid)
      real(wp), intent(in) :: w_phi(:)
      !> Adjoint of the merged tail charges (natom)
      real(wp), intent(in) :: w_q(:)
      !> Host response items
      type(response_type), intent(inout) :: response
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      !> Optional volume adjoint accumulator
      type(volume_adjoint_type), intent(inout), optional :: grid_adjoint
      !> Optional explicit nuclear gradient, Hartree/bohr (3, natom)
      real(wp), intent(inout), optional :: gradient(:, :)
      !> Scoped view of a host-data term and the unbound view handed to parameter terms
      type(coupling_view_type) :: view, unbound
      real(wp), allocatable :: w_qinf(:)
      logical :: contributes
      integer :: i
      call require_charge_adjoint(self, w_q, error)
      if (allocated(error)) return
      do i = 1, self%n_terms()
         associate (term => self%slots(i)%term)
            call term_charge_adjoint(self, i, w_q, w_qinf, contributes)
            if (.not. contributes) cycle
            if (term%is_host_fed()) then
               call coupling_make_view(coupling, i, view)
               call term%accumulate_adjoint(grid, view, w_phi, w_qinf, response, error, grid_adjoint, gradient)
               call coupling_close_view(coupling)
            else
               call term%accumulate_adjoint(grid, unbound, w_phi, w_qinf, response, error, grid_adjoint, gradient)
            end if
            deallocate (w_qinf)
            if (allocated(error)) return
         end associate
      end do
   end subroutine accumulate_adjoint_3d

   !> Require a charge adjoint per atom of the side, or none
   !>
   !> @param[in]  self   updated potential
   !> @param[in]  w_q    adjoint of the merged tail charges
   !> @param[out] error  adjoint of another length
   subroutine require_charge_adjoint(self, w_q, error)
      !> Updated potential
      class(moz_potential_type), intent(in) :: self
      !> Adjoint of the merged tail charges
      real(wp), intent(in) :: w_q(:)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      if (size(w_q) /= 0 .and. size(w_q) /= self%nat) then
         call fatal_error(error, "Charge adjoint does not match the atoms of the side")
      end if
   end subroutine require_charge_adjoint

   !> Whether term i receives the adjoint, and its masked charge adjoint
   !>
   !> - Adjoint recipients: field suppliers or charge owners
   !> - Charge adjoint: w_q masked to owned atoms
   !> - No owned atoms or no side charges: zero-sized charge adjoint
   !>
   !> @param[in]  self         updated potential
   !> @param[in]  i            term index
   !> @param[in]  w_q          adjoint of the merged tail charges (natom); zero-sized when the side has none
   !> @param[out] w_qinf       charge adjoint of the term (natom or 0); allocated only when it contributes
   !> @param[out] contributes  whether the term receives the adjoint
   subroutine term_charge_adjoint(self, i, w_q, w_qinf, contributes)
      !> Updated potential
      class(moz_potential_type), intent(in) :: self
      !> Term index
      integer, intent(in) :: i
      !> Adjoint of the merged tail charges (natom)
      real(wp), intent(in) :: w_q(:)
      !> Charge adjoint of the term
      real(wp), allocatable, intent(out) :: w_qinf(:)
      !> Whether the term receives the adjoint
      logical, intent(out) :: contributes
      logical :: owns
      owns = any(self%charge_owner == i)
      contributes = owns .or. self%slots(i)%term%has_field()
      if (.not. contributes) return
      if (owns .and. size(w_q) > 0) then
         w_qinf = merge(w_q, 0.0_wp, self%charge_owner == i)
      else
         allocate (w_qinf(0))
      end if
   end subroutine term_charge_adjoint

end submodule moist_model_moz_potential_compute
