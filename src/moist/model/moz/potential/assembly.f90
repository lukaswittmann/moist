!> Assembly of a MOZ potential: terms, update, coverage merge, settings and site table
submodule(moist_model_moz_potential) moist_model_moz_potential_assembly
   use mctc_env, only: fatal_error
   use mctc_io_convert, only: autokcal
   use moist_channels_coupling, only: coupling_view_type, coupling_make_view, coupling_close_view, coupling_set_scope
   use moist_model_moz_potential_term, only: charge_label_len, atom_label, atom_label_len
   use moist_model_moz_potential_lj_base, only: lj_label_len, check_lj_mixing
   use moist_model_moz_potential_kernel_ng, only: require_split_alpha
   implicit none(type, external)

   !> Last change stamp handed out to any potential
   integer(int64) :: last_stamp = 0_int64

   !> Largest accepted difference between the charge sum and the total charge, e
   real(wp), parameter :: total_charge_tolerance = 1.0e-4_wp

   !> Largest difference between the charge sum and the total charge that `neutralize` shifts away, e
   real(wp), parameter :: neutralize_limit = 1.0e-3_wp

contains

   !> Next change stamp, unique across all potentials and threads
   function next_stamp() result(stamp)
      !> New stamp
      integer(int64) :: stamp
      !$omp atomic capture
      last_stamp = last_stamp + 1_int64
      stamp = last_stamp
      !$omp end atomic
   end function next_stamp

   module subroutine potential_add(self, term, error)
      class(moz_potential_type), intent(inout) :: self
      class(potential_term_type), intent(in) :: term
      type(error_type), allocatable, intent(out) :: error
      type(potential_slot_type), allocatable :: slots(:)
      integer :: n, i, stat
      if (self%refuse_host_fed .and. term%is_host_fed()) then
         call fatal_error(error, "Potential term '"//trim(term%name())// &
            & "' reads host data and cannot serve a solvent")
         return
      end if
      n = self%n_terms()
      allocate (slots(n + 1))
      do i = 1, n
         call move_alloc(self%slots(i)%term, slots(i)%term)
      end do
      allocate (slots(n + 1)%term, source=term, stat=stat)
      if (stat /= 0) then
         ! Hand the existing terms back before failing
         do i = 1, n
            call move_alloc(slots(i)%term, self%slots(i)%term)
         end do
         call fatal_error(error, "Failed to copy potential term '"//trim(term%name())//"'")
         return
      end if
      call move_alloc(slots, self%slots)
      call drop_update(self)
   end subroutine potential_add

   module subroutine potential_update(self, mol, error, solvent_id)
      class(moz_potential_type), intent(inout) :: self
      class(structure_type), intent(in) :: mol
      type(error_type), allocatable, intent(out) :: error
      integer, intent(in), optional :: solvent_id
      integer :: i
      call drop_update(self)
      if (mol%nat < 1) then
         call fatal_error(error, "Potential update needs at least one atom")
         return
      end if
      if (self%n_terms() < 1) then
         call fatal_error(error, "Potential has no terms; add one before update")
         return
      end if
      do i = 1, self%n_terms()
         call self%slots(i)%term%update(mol, error, solvent_id)
         if (allocated(error)) return
         call require_contribution(self%slots(i)%term, mol%nat, error)
         if (allocated(error)) return
      end do
      call merge_pair(self, mol, error)
      if (.not. allocated(error)) call merge_charges(self, mol, error)
      if (allocated(error)) then
         call drop_update(self)
         return
      end if
      self%mol = mol
      self%nat = mol%nat
   end subroutine potential_update

   module subroutine potential_update_1d(self, mol, grid, error)
      class(moz_potential_type), intent(inout) :: self
      class(structure_type), intent(in) :: mol
      type(moist_math_grid_radial_type), intent(in) :: grid
      type(error_type), allocatable, intent(out) :: error
      integer :: i
      call potential_update(self, mol, error)
      if (allocated(error)) return
      do i = 1, self%n_terms()
         call self%slots(i)%term%update_grid(grid, error)
         if (allocated(error)) then
            call drop_update(self)
            return
         end if
      end do
   end subroutine potential_update_1d

   module subroutine potential_update_3d(self, mol, grid, error)
      class(moz_potential_type), intent(inout) :: self
      class(structure_type), intent(in) :: mol
      class(moist_math_grid_3d_type), intent(in) :: grid
      type(error_type), allocatable, intent(out) :: error
      integer :: i
      call potential_update(self, mol, error)
      if (allocated(error)) return
      do i = 1, self%n_terms()
         call self%slots(i)%term%update_grid(grid, error)
         if (allocated(error)) then
            call drop_update(self)
            return
         end if
      end do
   end subroutine potential_update_3d

   module subroutine potential_declare_1d(self, grid, coupling, error)
      class(moz_potential_type), intent(inout) :: self
      type(moist_math_grid_radial_type), intent(in) :: grid
      type(coupling_type), intent(inout) :: coupling
      type(error_type), allocatable, intent(out) :: error
      integer :: i
      ! TODO: skip terms whose charges are all shadowed by earlier owners; the host answers unread requests
      do i = 1, self%n_terms()
         call coupling_set_scope(coupling, i)
         call self%slots(i)%term%declare_pass(grid, coupling, error)
         if (allocated(error)) return
      end do
   end subroutine potential_declare_1d

   module subroutine potential_declare_3d(self, grid, coupling, error)
      class(moz_potential_type), intent(inout) :: self
      class(moist_math_grid_3d_type), intent(in) :: grid
      type(coupling_type), intent(inout) :: coupling
      type(error_type), allocatable, intent(out) :: error
      integer :: i
      ! TODO: skip terms whose charges are all shadowed by earlier owners; the host answers unread requests
      do i = 1, self%n_terms()
         call coupling_set_scope(coupling, i)
         call self%slots(i)%term%declare_pass(grid, coupling, error)
         if (allocated(error)) return
      end do
   end subroutine potential_declare_3d

   module subroutine potential_refuse_host_fed_terms(self)
      class(moz_potential_type), intent(inout) :: self
      self%refuse_host_fed = .true.
   end subroutine potential_refuse_host_fed_terms

   module subroutine potential_set_ng_split(self, ng_split, alpha, error)
      class(moz_potential_type), intent(inout) :: self
      logical, intent(in) :: ng_split
      real(wp), intent(in) :: alpha
      type(error_type), allocatable, intent(out) :: error
      call require_split_alpha(ng_split, alpha, error)
      if (allocated(error)) return
      self%split = ng_split
      self%split_alpha = alpha
      self%generation_stamp = next_stamp()
   end subroutine potential_set_ng_split

   module subroutine potential_set_mixing(self, mixing, error)
      class(moz_potential_type), intent(inout) :: self
      integer, intent(in) :: mixing
      type(error_type), allocatable, intent(out) :: error
      call check_lj_mixing(mixing, error)
      if (allocated(error)) return
      self%mixing_rule = mixing
      self%generation_stamp = next_stamp()
   end subroutine potential_set_mixing

   pure module function potential_ng_split(self) result(ng_split)
      class(moz_potential_type), intent(in) :: self
      logical :: ng_split
      ng_split = self%split
   end function potential_ng_split

   pure module function potential_alpha(self) result(alpha)
      class(moz_potential_type), intent(in) :: self
      real(wp) :: alpha
      alpha = self%split_alpha
   end function potential_alpha

   pure module function potential_mixing(self) result(mixing)
      class(moz_potential_type), intent(in) :: self
      integer :: mixing
      mixing = self%mixing_rule
   end function potential_mixing

   pure module function potential_n_terms(self) result(n)
      class(moz_potential_type), intent(in) :: self
      integer :: n
      n = 0
      if (allocated(self%slots)) n = size(self%slots)
   end function potential_n_terms

   pure module function potential_natom(self) result(natom)
      class(moz_potential_type), intent(in) :: self
      integer :: natom
      natom = self%nat
   end function potential_natom

   pure module function potential_is_updated(self) result(updated)
      class(moz_potential_type), intent(in) :: self
      logical :: updated
      updated = self%nat > 0
   end function potential_is_updated

   pure module function potential_generation(self) result(generation)
      class(moz_potential_type), intent(in) :: self
      integer(int64) :: generation
      generation = self%generation_stamp
   end function potential_generation

   module subroutine potential_sites(self, sites, error, coupling)
      class(moz_potential_type), intent(in) :: self
      type(potential_sites_type), intent(out) :: sites
      type(error_type), allocatable, intent(out) :: error
      class(coupling_type), intent(inout), target, optional :: coupling
      !> Scoped view of a host-data term and the unbound view handed to parameter terms
      type(coupling_view_type) :: view, unbound
      real(wp), allocatable :: q_term(:)
      integer :: i
      call self%require_updated(error)
      if (allocated(error)) return
      sites%natom = self%nat
      sites%has_field = self%has_field()
      if (allocated(self%pair)) sites%pair = self%pair
      if (.not. any(self%charge_owner > 0)) return
      allocate (sites%q(self%nat), source=0.0_wp)
      do i = 1, self%n_terms()
         if (.not. any(self%charge_owner == i)) cycle
         associate (term => self%slots(i)%term)
            if (term%is_host_fed()) then
               if (.not. present(coupling)) then
                  call fatal_error(error, "Potential term '"//trim(term%name())// &
                     & "' reads host data; its tail charges need the coupling")
                  return
               end if
               call coupling_make_view(coupling, i, view)
               call term%tail_charges(view, q_term, error)
               call coupling_close_view(coupling)
            else
               call term%tail_charges(unbound, q_term, error)
            end if
            if (allocated(error)) return
            if (.not. allocated(q_term)) then
               call fatal_error(error, "Potential term '"//trim(term%name())//"' supplied no tail charges")
               return
            end if
            if (size(q_term) /= self%nat) then
               call fatal_error(error, "Tail charges of potential term '"//trim(term%name())// &
                  & "' do not match the atoms of its side")
               return
            end if
            where (self%charge_owner == i) sites%q = q_term
            deallocate (q_term)
         end associate
      end do
      sites%q = sites%q - self%charge_shift
      if (.not. host_charges(self)) return
      if (abs(sum(sites%q) - self%mol%charge) > total_charge_tolerance) then
         call charge_sum_error(sum(sites%q), self%mol%charge, error)
      end if
   end subroutine potential_sites

   module subroutine potential_check_total_charge(self, error)
      class(moz_potential_type), intent(inout) :: self
      type(error_type), allocatable, intent(out) :: error
      type(potential_sites_type) :: sites
      real(wp) :: total
      call self%require_updated(error)
      if (allocated(error)) return
      if (host_charges(self)) return
      call self%sites(sites, error)
      if (allocated(error)) return
      total = 0.0_wp
      if (allocated(sites%q)) total = sum(sites%q)
      if (abs(total - self%mol%charge) > total_charge_tolerance) then
         call charge_sum_error(total, self%mol%charge, error)
         call drop_update(self)
      end if
   end subroutine potential_check_total_charge

   module subroutine potential_neutralize(self, shift, error)
      class(moz_potential_type), intent(inout) :: self
      real(wp), intent(out) :: shift
      type(error_type), allocatable, intent(out) :: error
      type(potential_sites_type) :: sites
      real(wp) :: mismatch
      shift = 0.0_wp
      call self%require_updated(error)
      if (allocated(error)) return
      if (host_charges(self)) then
         call fatal_error(error, "Charges from host data cannot be neutralized")
         return
      end if
      self%charge_shift = 0.0_wp
      call self%sites(sites, error)
      if (allocated(error)) return
      if (.not. allocated(sites%q)) return
      mismatch = sum(sites%q) - self%mol%charge
      ! Written so that a NaN sum fails too
      if (.not. abs(mismatch) <= neutralize_limit) then
         call charge_sum_error(sum(sites%q), self%mol%charge, error)
         call drop_update(self)
         return
      end if
      shift = mismatch/real(self%nat, wp)
      self%charge_shift = shift
      self%generation_stamp = next_stamp()
   end subroutine potential_neutralize

   module subroutine potential_print_table(self, unit, error)
      class(moz_potential_type), intent(in) :: self
      integer, intent(in) :: unit
      type(error_type), allocatable, intent(out) :: error
      type(coupling_view_type) :: unbound
      real(wp), allocatable :: q(:), q_term(:)
      character(len=12) :: q_col
      character(len=charge_label_len) :: q_src
      character(len=16) :: sigma_col, eps_col
      character(len=lj_label_len) :: label
      integer :: iat, i, owner
      call self%require_updated(error)
      if (allocated(error)) return
      allocate (q(self%nat), source=0.0_wp)
      do i = 1, self%n_terms()
         if (.not. any(self%charge_owner == i)) cycle
         if (self%slots(i)%term%is_host_fed()) cycle
         call self%slots(i)%term%tail_charges(unbound, q_term, error)
         if (allocated(error)) return
         where (self%charge_owner == i) q = q_term
      end do
      q = q - self%charge_shift
      write (unit, "(2a7, a12, 2x, a14, 2a16, 2x, a)") "N", "Sym", "qat", "q/Src         ", "sigma (bohr)", &
         & "eps (kcal/mol)", "Type/Src"
      do iat = 1, self%nat
         owner = self%charge_owner(iat)
         q_src = "-"
         if (owner == 0) then
            q_col = adjustr("-           ")
         else
            q_src = self%slots(owner)%term%charge_label()
            if (self%slots(owner)%term%is_host_fed()) then
               q_col = adjustr("host        ")
            else
               write (q_col, "(f12.4)") q(iat)
            end if
         end if
         if (allocated(self%pair)) then
            write (sigma_col, "(f16.4)") self%pair%sigma(iat)
            write (eps_col, "(f16.4)") self%pair%epsilon(iat)*autokcal
            label = self%pair%label(iat)
         else
            sigma_col = adjustr("-               ")
            eps_col = sigma_col
            label = "-"
         end if
         write (unit, "(i7, a7, a12, 2x, a14, 2a16, 2x, a)") iat, trim(self%mol%sym(self%mol%id(iat))), q_col, &
            & q_src, sigma_col, eps_col, trim(label)
      end do
   end subroutine potential_print_table

   module subroutine potential_require_updated(self, error)
      class(moz_potential_type), intent(in) :: self
      type(error_type), allocatable, intent(out) :: error
      if (self%nat < 1) call fatal_error(error, "Potential is not updated; update it after adding terms")
   end subroutine potential_require_updated

   pure module function potential_reads_host_data(self) result(host_fed)
      class(moz_potential_type), intent(in) :: self
      logical :: host_fed
      integer :: i
      host_fed = .false.
      do i = 1, self%n_terms()
         host_fed = host_fed .or. self%slots(i)%term%is_host_fed()
      end do
   end function potential_reads_host_data

   pure module function potential_has_field(self) result(has)
      class(moz_potential_type), intent(in) :: self
      logical :: has
      integer :: i
      has = .false.
      do i = 1, self%n_terms()
         has = has .or. self%slots(i)%term%has_field()
      end do
   end function potential_has_field

   !> Drop the updated state, keeping the terms and settings; takes a fresh change stamp
   !>
   !> @param[in,out] self  potential
   subroutine drop_update(self)
      !> Potential
      class(moz_potential_type), intent(inout) :: self
      self%nat = 0
      if (allocated(self%pair)) deallocate (self%pair)
      if (allocated(self%charge_owner)) deallocate (self%charge_owner)
      self%charge_shift = 0.0_wp
      self%generation_stamp = next_stamp()
   end subroutine drop_update

   !> Whether a host-data term owns the charge of any atom
   !>
   !> @param[in] self  updated potential
   pure function host_charges(self) result(host)
      !> Updated potential
      class(moz_potential_type), intent(in) :: self
      !> Whether host data supplies charges
      logical :: host
      integer :: i
      host = .false.
      do i = 1, self%n_terms()
         if (any(self%charge_owner == i)) host = host .or. self%slots(i)%term%is_host_fed()
      end do
   end function host_charges

   !> Error for charges not summing to the total charge of the structure
   !>
   !> @param[in]  total   charge sum, e
   !> @param[in]  charge  total charge of the structure, e
   !> @param[out] error   mismatch message
   subroutine charge_sum_error(total, charge, error)
      !> Charge sum, e
      real(wp), intent(in) :: total
      !> Total charge of the structure, e
      real(wp), intent(in) :: charge
      !> Mismatch message
      type(error_type), allocatable, intent(out) :: error
      character(len=24) :: sum_text, charge_text
      write (sum_text, "(f24.6)") total
      write (charge_text, "(f24.6)") charge
      call fatal_error(error, "Charges sum to "//trim(adjustl(sum_text))// &
         & " e, but the total charge of the structure is "//trim(adjustl(charge_text))//" e")
   end subroutine charge_sum_error

   !> Require nonempty updated term data, with `natom` LJ entries when present
   !>
   !> - Host terms: data through coupling
   !> - Parameter terms: LJ data, field or tail charges required
   !>
   !> @param[in]  term   updated term
   !> @param[in]  natom  number of atoms of the side
   !> @param[out] error  term supplying nothing or LJ data of the wrong length
   subroutine require_contribution(term, natom, error)
      !> Updated term
      class(potential_term_type), intent(in) :: term
      !> Number of atoms of the side
      integer, intent(in) :: natom
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      type(coupling_view_type) :: unbound
      real(wp), allocatable :: q(:)
      if (allocated(term%pair)) then
         if (size(term%pair%sigma) /= natom .or. size(term%pair%epsilon) /= natom &
            & .or. size(term%pair%covered) /= natom) then
            call fatal_error(error, "Lennard-Jones data of potential term '"//trim(term%name())// &
               & "' does not cover the atoms of its side")
         end if
         return
      end if
      if (term%is_host_fed() .or. term%has_field()) return
      call term%tail_charges(unbound, q, error)
      if (allocated(error)) return
      if (.not. allocated(q)) then
         call fatal_error(error, "Potential term '"//trim(term%name())//"' supplies no contribution")
      else if (size(q) /= natom) then
         call fatal_error(error, "Tail charges of potential term '"//trim(term%name())// &
            & "' do not match the atoms of its side")
      end if
   end subroutine require_contribution

   !> Merge the Lennard-Jones data of the updated terms in insertion order
   !>
   !> - First LJ term seeds `pair`; later terms fill uncovered atoms
   !> - Atoms still uncovered: named error
   !>
   !> @param[in,out] self   potential with updated terms
   !> @param[in]     mol    structure of this side, for the atom labels of the error
   !> @param[out]    error  uncovered atoms
   subroutine merge_pair(self, mol, error)
      !> Potential with updated terms
      class(moz_potential_type), intent(inout) :: self
      !> Structure of this side
      class(structure_type), intent(in) :: mol
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      character(len=:), allocatable :: missing
      character(len=atom_label_len) :: label
      logical, allocatable :: fill(:)
      integer :: i, iat
      do i = 1, self%n_terms()
         associate (term => self%slots(i)%term)
            if (.not. allocated(term%pair)) cycle
            if (.not. allocated(self%pair)) then
               self%pair = term%pair
               cycle
            end if
            fill = term%pair%covered .and. .not. self%pair%covered
            where (fill)
               self%pair%sigma = term%pair%sigma
               self%pair%epsilon = term%pair%epsilon
               self%pair%label = term%pair%label
               self%pair%covered = .true.
            end where
         end associate
      end do
      if (.not. allocated(self%pair)) return
      if (all(self%pair%covered)) return
      missing = ""
      do iat = 1, mol%nat
         if (self%pair%covered(iat)) cycle
         call atom_label(mol, iat, label)
         if (len(missing) > 0) missing = missing//", "
         missing = missing//trim(label)
      end do
      call fatal_error(error, "No pair potential parameters for atoms: "//missing)
   end subroutine merge_pair

   !> Assign each atom the first term, in insertion order, with a charge for it
   !>
   !> - Any charge term present, even one covering no atom: full atom coverage required,
   !>   named error for uncovered atoms
   !> - Whole-side grid field: no other charge supplier
   !>
   !> @param[in,out] self   potential with updated terms
   !> @param[in]     mol    structure of this side, for the atom names of the error
   !> @param[out]    error  atoms without a charge or a field term beside another charge owner
   subroutine merge_charges(self, mol, error)
      !> Potential with updated terms
      class(moz_potential_type), intent(inout) :: self
      !> Structure of this side
      class(structure_type), intent(in) :: mol
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      character(len=:), allocatable :: missing
      character(len=atom_label_len) :: label
      integer :: i, iat
      allocate (self%charge_owner(mol%nat), source=0)
      do i = 1, self%n_terms()
         where (self%charge_owner == 0 .and. self%slots(i)%term%charges_covered(mol%nat)) self%charge_owner = i
      end do
      do i = 1, self%n_terms()
         if (.not. self%slots(i)%term%has_field()) cycle
         if (.not. any(self%charge_owner > 0 .and. self%charge_owner /= i)) cycle
         call fatal_error(error, "Potential term '"//trim(self%slots(i)%term%name())// &
            & "' supplies the field of the whole side and must supply every charge; "// &
            & "add it before any other charge term")
         return
      end do
      if (all(self%charge_owner > 0)) return
      if (all(self%charge_owner == 0)) then
         ! Without a charge term the side has no charges; a charge term must cover every atom
         do i = 1, self%n_terms()
            if (len_trim(self%slots(i)%term%charge_label()) > 0) exit
         end do
         if (i > self%n_terms()) return
      end if
      missing = ""
      do iat = 1, mol%nat
         if (self%charge_owner(iat) > 0) cycle
         call atom_label(mol, iat, label)
         if (len(missing) > 0) missing = missing//", "
         missing = missing//trim(label)
      end do
      call fatal_error(error, "No charges for atoms: "//missing)
   end subroutine merge_charges

end submodule moist_model_moz_potential_assembly
