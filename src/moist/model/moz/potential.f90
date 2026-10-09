!> MOZ potential: the ordered terms of one side and the interactions they enter
!>
!> Owned by the solvent (`solvent%potential`) and by the 1D and 3D models
!> (`model%potential`); one type serves both sides:
!>
!>    call solvent%potential%add(lj_spce, error)
!>    call solvent%potential%add(coulomb_spce, error)
!>    call model%potential%add(lj_gaff, error)
!>    call model%potential%add(lj_uff, error)
!>
!> - `add`: copy of a term; insertion order is the order of preference
!> - `update(mol)`: every term from the structure, coverage merge, structure kept
!>   - LJ data and charges: first term covering an atom supplies it, never summed
!>   - Uncovered atoms: named error
!>   - Whole-side field term: no other charge supplier
!> - `add` drops the updated state; `compute` refuses until the next `update`
!> - Interaction settings: LJ mixing rule, Ng split flag and alpha
!>   - Defaults: Lorentz-Berthelot, split on, alpha = 1/bohr
!>   - The potential `compute` is called on decides
!> - `compute`: Ng-split tables; VV without `solvent`, UV with it
!>   - 1D: radial grids, column layout from `site_pairs`
!>   - 3D: volume grid, UV only, one column per solvent site
!>   - Host-data terms read their coupling scope; scope of term i is i
!> - `compute_adjoint` (UV): field and charge adjoints back to every contributing term
!> - Owners do not watch the potential: update (and solve) them again after a change
!> - Re-exports every term type and preset
!>   - Typed Lennard-Jones: `lj_gaff`, `lj_oplsaa`
!>   - Element Lennard-Jones: `lj_uff`, `lj_dreiding`, `lj_tm`
!>   - Water-site Lennard-Jones: `lj_spce`
!>   - Solvent charges: `coulomb_<model>_<environment>`, `coulomb_spce`
!>   - Table charge models: hirshfeld, resp, mbis, chelpg
!>   - Table environments: gas, solvent, conductor
!>   - Solvent multipoles: `multipoles_mbis_<environment>`, stub for 6D MOZ
!>   - User-data constructors: `new_fixed_charges`, `new_custom_lj`
!>   - Solute electrostatics: `monopole_type`, `multipole_type`, `ec_charges_type`
module moist_model_moz_potential
   use mctc_env, only: wp, error_type
   use mctc_io, only: structure_type
   use moist_channels_coupling, only: coupling_type
   use moist_channels_response, only: response_type
   use moist_math_grid_radial_grid, only: moist_math_grid_radial_type
   use moist_math_grid_3d_base, only: moist_math_grid_3d_type
   use moist_math_grid_3d_adjoint, only: volume_adjoint_type
   use moist_model_moz_potential_term, only: potential_term_type
   use moist_model_moz_potential_sites, only: potential_sites_type
   use moist_model_moz_potential_lj_base, only: lj_12_6_type, lj_mixing_lorentz_berthelot, lj_mixing_geometric
   use moist_model_moz_potential_lj_typed, only: lj_typed_type, lj_gaff, lj_oplsaa
   use moist_model_moz_potential_lj_element, only: lj_element_type, lj_uff, lj_dreiding, lj_tm
   use moist_model_moz_potential_lj_solvent, only: lj_solvent_type, lj_spce
   use moist_model_moz_potential_lj_custom, only: custom_lj_type, new_custom_lj
   use moist_model_moz_potential_electrostatic_fixed, only: fixed_charges_type, new_fixed_charges
   use moist_model_moz_potential_electrostatic_monopole, only: monopole_type
   use moist_model_moz_potential_electrostatic_multipole, only: multipole_type
   use moist_model_moz_potential_electrostatic_ec, only: ec_charges_type
   use moist_model_moz_potential_electrostatic_solvent, only: solvent_charges_type, solvent_model_charges_type, &
      & solvent_multipoles_type, coulomb_hirshfeld_gas, coulomb_hirshfeld_solvent, coulomb_hirshfeld_conductor, &
      & coulomb_resp_gas, coulomb_resp_solvent, coulomb_resp_conductor, coulomb_mbis_gas, coulomb_mbis_solvent, &
      & coulomb_mbis_conductor, coulomb_chelpg_gas, coulomb_chelpg_solvent, coulomb_chelpg_conductor, coulomb_spce, &
      & multipoles_mbis_gas, multipoles_mbis_solvent, multipoles_mbis_conductor
   implicit none(type, external)
   private

   public :: moz_potential_type
   public :: potential_term_type, potential_sites_type
   public :: lj_mixing_lorentz_berthelot, lj_mixing_geometric
   public :: lj_typed_type, lj_element_type, lj_solvent_type, custom_lj_type, new_custom_lj
   public :: lj_gaff, lj_oplsaa, lj_uff, lj_dreiding, lj_tm, lj_spce
   public :: fixed_charges_type, new_fixed_charges
   public :: monopole_type, multipole_type, ec_charges_type
   public :: solvent_charges_type, solvent_model_charges_type, solvent_multipoles_type
   public :: coulomb_hirshfeld_gas, coulomb_hirshfeld_solvent, coulomb_hirshfeld_conductor
   public :: coulomb_resp_gas, coulomb_resp_solvent, coulomb_resp_conductor
   public :: coulomb_mbis_gas, coulomb_mbis_solvent, coulomb_mbis_conductor
   public :: coulomb_chelpg_gas, coulomb_chelpg_solvent, coulomb_chelpg_conductor
   public :: coulomb_spce
   public :: multipoles_mbis_gas, multipoles_mbis_solvent, multipoles_mbis_conductor

   !> One owned term
   type :: potential_slot_type
      !> Term
      class(potential_term_type), allocatable :: term
   end type potential_slot_type

   !> Ordered terms of one side and the settings of its interactions
   type :: moz_potential_type
      private
      !> Owned terms in insertion order; the coupling scope of term i is i
      type(potential_slot_type), allocatable :: slots(:)
      !> Whether host-data terms are refused (solvent side)
      logical :: refuse_host_fed = .false.
      !> Atoms of the latest update; 0 until updated
      integer :: nat = 0
      !> Structure of the latest update, bohr
      type(structure_type) :: mol
      !> Merged Lennard-Jones data, complete after an update; absent when no term has any
      type(lj_12_6_type), allocatable :: pair
      !> Term supplying the charge of each atom, 0 when no term has charges (natom)
      integer, allocatable :: charge_owner(:)
      !> Whether the Coulomb tail is split off
      logical :: split = .true.
      !> Ng split parameter, 1/bohr
      real(wp) :: split_alpha = 1.0_wp
      !> Lennard-Jones mixing rule
      integer :: mixing_rule = lj_mixing_lorentz_berthelot
   contains
      !> Append a copy of a term; drops the updated state
      procedure :: add => potential_add
      !> Update every term from the structure and merge their data
      procedure :: update => potential_update
      !> Declare the coupling requests of every term in its own scope (radial grid)
      procedure :: declare_1d => potential_declare_1d
      !> Declare the coupling requests of every term in its own scope (volume grid)
      procedure :: declare_3d => potential_declare_3d
      generic :: declare => declare_1d, declare_3d
      !> Refuse host-data terms from now on (solvent side)
      procedure :: refuse_host_fed_terms => potential_refuse_host_fed_terms
      !> Set the Ng split of the Coulomb tail
      procedure :: set_ng_split => potential_set_ng_split
      !> Set the Lennard-Jones mixing rule
      procedure :: set_mixing => potential_set_mixing
      !> Whether the Coulomb tail is split off
      procedure :: ng_split => potential_ng_split
      !> Ng split parameter, 1/bohr
      procedure :: alpha => potential_alpha
      !> Lennard-Jones mixing rule
      procedure :: mixing => potential_mixing
      !> Number of terms
      procedure :: n_terms => potential_n_terms
      !> Number of atoms of the latest update; 0 until updated
      procedure :: natom => potential_natom
      !> Whether the potential is updated
      procedure :: is_updated => potential_is_updated
      !> Resolved per-atom data of the latest update
      procedure :: sites => potential_sites
      !> Print the per-atom site parameters as a table
      procedure :: print_table => potential_print_table
      !> Column layout of the 1D tables
      procedure :: site_pairs => potential_site_pairs
      !> Ng-split tables on radial grids, VV or UV
      procedure :: compute_1d => potential_compute_1d
      !> Ng-split UV tables on a volume grid
      procedure :: compute_3d => potential_compute_3d
      generic :: compute => compute_1d, compute_3d
      !> Reverse mode of the 1D UV tables
      procedure :: compute_adjoint_1d => potential_compute_adjoint_1d
      !> Reverse mode of the 3D UV tables
      procedure :: compute_adjoint_3d => potential_compute_adjoint_3d
      generic :: compute_adjoint => compute_adjoint_1d, compute_adjoint_3d
      !> Require an updated potential
      procedure, private :: require_updated => potential_require_updated
      !> Whether any term reads host data
      procedure, private :: reads_host_data => potential_reads_host_data
      !> Whether any term supplies a grid field
      procedure, private :: has_field => potential_has_field
   end type moz_potential_type

   interface

      !* ============================================================================== *!
      !*                              Assembly (assembly.f90)                            *!
      !* ============================================================================== *!

      !> Append a copy of a term
      !>
      !> - Solvent side: host-data terms refused
      !> - Updated state dropped; `compute` refuses until the next `update`
      !>
      !> @param[in,out] self   potential
      !> @param[in]     term   term to copy
      !> @param[out]    error  host-data term on a solvent side or copy failure
      module subroutine potential_add(self, term, error)
         implicit none(type, external)
         !> Potential
         class(moz_potential_type), intent(inout) :: self
         !> Term to copy
         class(potential_term_type), intent(in) :: term
         !> Error handling
         type(error_type), allocatable, intent(out) :: error
      end subroutine potential_add

      !> Update every term from the structure and merge their data
      !>
      !> - Previous update dropped first; a failure leaves the potential not updated
      !> - Structure kept for the positions (3D) and the site table
      !>
      !> @param[in,out] self        potential
      !> @param[in]     mol         structure of this side, bohr
      !> @param[out]    error       empty structure, no terms, term failure, a term supplying nothing or uncovered atoms
      !> @param[in]     solvent_id  optional solvent id of this side
      module subroutine potential_update(self, mol, error, solvent_id)
         implicit none(type, external)
         !> Potential
         class(moz_potential_type), intent(inout) :: self
         !> Structure of this side
         class(structure_type), intent(in) :: mol
         !> Error handling
         type(error_type), allocatable, intent(out) :: error
         !> Optional solvent id of this side
         integer, intent(in), optional :: solvent_id
      end subroutine potential_update

      !> Declare the coupling requests of every term in its own scope (radial grid)
      !>
      !> - Coupling scope of term i: i; no requests from parameter terms
      !> - Registration opened by the owning model
      !>
      !> @param[in,out] self      potential
      !> @param[in]     grid      radial grid of the owning model
      !> @param[in,out] coupling  coupling in registration
      !> @param[out]    error     failed declaration
      module subroutine potential_declare_1d(self, grid, coupling, error)
         implicit none(type, external)
         !> Potential
         class(moz_potential_type), intent(inout) :: self
         !> Radial grid of the owning model
         type(moist_math_grid_radial_type), intent(in) :: grid
         !> Coupling in registration
         type(coupling_type), intent(inout) :: coupling
         !> Error handling
         type(error_type), allocatable, intent(out) :: error
      end subroutine potential_declare_1d

      !> Declare the coupling requests of every term in its own scope (volume grid)
      !>
      !> - Coupling scope of term i: i; no requests from parameter terms
      !> - Registration opened by the owning model
      !>
      !> @param[in,out] self      potential
      !> @param[in]     grid      volume grid of the owning model
      !> @param[in,out] coupling  coupling in registration
      !> @param[out]    error     failed declaration
      module subroutine potential_declare_3d(self, grid, coupling, error)
         implicit none(type, external)
         !> Potential
         class(moz_potential_type), intent(inout) :: self
         !> Volume grid of the owning model
         class(moist_math_grid_3d_type), intent(in) :: grid
         !> Coupling in registration
         type(coupling_type), intent(inout) :: coupling
         !> Error handling
         type(error_type), allocatable, intent(out) :: error
      end subroutine potential_declare_3d

      !> Refuse host-data terms from now on (solvent side)
      !>
      !> @param[in,out] self  potential
      module subroutine potential_refuse_host_fed_terms(self)
         implicit none(type, external)
         !> Potential
         class(moz_potential_type), intent(inout) :: self
      end subroutine potential_refuse_host_fed_terms

      !> Set the Ng split of the Coulomb tail
      !>
      !> @param[in,out] self      potential
      !> @param[in]     ng_split  whether the Coulomb tail is split off
      !> @param[in]     alpha     Ng split parameter, 1/bohr; positive and finite when the split is on
      !> @param[out]    error     invalid split parameter; the previous split is kept
      module subroutine potential_set_ng_split(self, ng_split, alpha, error)
         implicit none(type, external)
         !> Potential
         class(moz_potential_type), intent(inout) :: self
         !> Whether the Coulomb tail is split off
         logical, intent(in) :: ng_split
         !> Ng split parameter, 1/bohr
         real(wp), intent(in) :: alpha
         !> Error handling
         type(error_type), allocatable, intent(out) :: error
      end subroutine potential_set_ng_split

      !> Set the Lennard-Jones mixing rule
      !>
      !> @param[in,out] self    potential
      !> @param[in]     mixing  `lj_mixing_lorentz_berthelot` or `lj_mixing_geometric`
      !> @param[out]    error   unknown mixing rule; the previous rule is kept
      module subroutine potential_set_mixing(self, mixing, error)
         implicit none(type, external)
         !> Potential
         class(moz_potential_type), intent(inout) :: self
         !> Mixing rule
         integer, intent(in) :: mixing
         !> Error handling
         type(error_type), allocatable, intent(out) :: error
      end subroutine potential_set_mixing

      !> Whether the Coulomb tail is split off
      !>
      !> @param[in] self  potential
      pure module function potential_ng_split(self) result(ng_split)
         implicit none(type, external)
         !> Potential
         class(moz_potential_type), intent(in) :: self
         !> Split flag
         logical :: ng_split
      end function potential_ng_split

      !> Ng split parameter, 1/bohr
      !>
      !> @param[in] self  potential
      pure module function potential_alpha(self) result(alpha)
         implicit none(type, external)
         !> Potential
         class(moz_potential_type), intent(in) :: self
         !> Split parameter, 1/bohr
         real(wp) :: alpha
      end function potential_alpha

      !> Lennard-Jones mixing rule
      !>
      !> @param[in] self  potential
      pure module function potential_mixing(self) result(mixing)
         implicit none(type, external)
         !> Potential
         class(moz_potential_type), intent(in) :: self
         !> Mixing rule, `lj_mixing_lorentz_berthelot` or `lj_mixing_geometric`
         integer :: mixing
      end function potential_mixing

      !> Number of terms
      !>
      !> @param[in] self  potential
      pure module function potential_n_terms(self) result(n)
         implicit none(type, external)
         !> Potential
         class(moz_potential_type), intent(in) :: self
         !> Count
         integer :: n
      end function potential_n_terms

      !> Number of atoms of the latest update; 0 until updated
      !>
      !> @param[in] self  potential
      pure module function potential_natom(self) result(natom)
         implicit none(type, external)
         !> Potential
         class(moz_potential_type), intent(in) :: self
         !> Atom count
         integer :: natom
      end function potential_natom

      !> Whether the potential is updated
      !>
      !> @param[in] self  potential
      pure module function potential_is_updated(self) result(updated)
         implicit none(type, external)
         !> Potential
         class(moz_potential_type), intent(in) :: self
         !> Update flag
         logical :: updated
      end function potential_is_updated

      !> Resolve the merged LJ data and tail charges of the latest update
      !>
      !> - Charge per atom from its owning term
      !> - Host owners through scoped view of `coupling`; refusal without coupling
      !>
      !> @param[in]     self      updated potential
      !> @param[out]    sites     resolved data of the side
      !> @param[out]    error     not updated, host-data owner without coupling, charge read failure
      !> @param[in,out] coupling  optional coupling of the model; host-data term i reads scope i
      module subroutine potential_sites(self, sites, error, coupling)
         implicit none(type, external)
         !> Updated potential
         class(moz_potential_type), intent(in) :: self
         !> Resolved data of the side
         type(potential_sites_type), intent(out) :: sites
         !> Error handling
         type(error_type), allocatable, intent(out) :: error
         !> Optional coupling of the model
         class(coupling_type), intent(inout), target, optional :: coupling
      end subroutine potential_sites

      !> Print the per-atom site parameters of the latest update as a table
      !>
      !> - One row per atom: index, symbol, charge (e), sigma (bohr), epsilon (kcal/mol)
      !> - Charge-source label and LJ "Type/Src" label from owning terms
      !> - Host charge: "host"; missing data: "-"
      !>
      !> @param[in]  self   updated potential
      !> @param[in]  unit   output unit
      !> @param[out] error  not updated
      module subroutine potential_print_table(self, unit, error)
         implicit none(type, external)
         !> Updated potential
         class(moz_potential_type), intent(in) :: self
         !> Output unit
         integer, intent(in) :: unit
         !> Error handling
         type(error_type), allocatable, intent(out) :: error
      end subroutine potential_print_table

      !> Require an updated potential
      !>
      !> @param[in]  self   potential
      !> @param[out] error  not updated
      module subroutine potential_require_updated(self, error)
         implicit none(type, external)
         !> Potential
         class(moz_potential_type), intent(in) :: self
         !> Error handling
         type(error_type), allocatable, intent(out) :: error
      end subroutine potential_require_updated

      !> Whether any term reads host data
      !>
      !> @param[in] self  potential
      pure module function potential_reads_host_data(self) result(host_fed)
         implicit none(type, external)
         !> Potential
         class(moz_potential_type), intent(in) :: self
         !> Host-data flag
         logical :: host_fed
      end function potential_reads_host_data

      !> Whether any term supplies a grid field
      !>
      !> @param[in] self  potential
      pure module function potential_has_field(self) result(has)
         implicit none(type, external)
         !> Potential
         class(moz_potential_type), intent(in) :: self
         !> Presence
         logical :: has
      end function potential_has_field

      !* ============================================================================== *!
      !*                               Compute (compute.f90)                             *!
      !* ============================================================================== *!

      !> Column layout of the 1D tables
      !>
      !> - VV (no `solvent`): packed upper triangle, column `packed_index(i, j)` holds `[i, j]`, i <= j
      !> - UV: every (solute atom, solvent atom) pair, solute major; column (iu - 1) nv + iv holds `[iu, iv]`
      !>
      !> @param[in]  self     updated potential
      !> @param[out] pairs    site pairs (2, npair)
      !> @param[out] error    not updated
      !> @param[in]  solvent  optional updated solvent potential; UV layout when present
      module subroutine potential_site_pairs(self, pairs, error, solvent)
         implicit none(type, external)
         !> Updated potential
         class(moz_potential_type), intent(in) :: self
         !> Site pairs (2, npair)
         integer, allocatable, intent(out) :: pairs(:, :)
         !> Error handling
         type(error_type), allocatable, intent(out) :: error
         !> Optional solvent potential
         class(moz_potential_type), intent(in), optional :: solvent
      end subroutine potential_site_pairs

      !> Ng-split tables on radial grids, one column per site pair of `site_pairs`
      !>
      !> - VV without `solvent`: this potential on both sides; no host-data terms
      !> - UV with `solvent`: this potential is the first side; the solvent must read no host data
      !> - Mixing rule and Ng split of this potential
      !> - Radial field of the first side summed from its terms (host terms through `coupling`)
      !>
      !> @param[in]     self      updated potential, first side
      !> @param[in]     rgrid     real-space radial grid
      !> @param[in]     kgrid     reciprocal radial grid
      !> @param[out]    u_sr      short-range potential, Hartree (npts, npair)
      !> @param[out]    ur_lr     long-range real-space potential, Hartree (npts, npair)
      !> @param[out]    uk_lr     long-range reciprocal potential, Hartree bohr^3 (nk, npair)
      !> @param[out]    error     not updated, host data on the wrong side, coupling, category or input failure
      !> @param[in]     solvent   optional updated solvent potential, second side; UV when present
      !> @param[in,out] coupling  optional coupling of the model; host-data term i reads scope i
      module subroutine potential_compute_1d(self, rgrid, kgrid, u_sr, ur_lr, uk_lr, error, solvent, coupling)
         implicit none(type, external)
         !> Updated potential, first side
         class(moz_potential_type), intent(in) :: self
         !> Real-space radial grid
         type(moist_math_grid_radial_type), intent(in) :: rgrid
         !> Reciprocal radial grid
         type(moist_math_grid_radial_type), intent(in) :: kgrid
         !> Short-range potential (npts, npair)
         real(wp), allocatable, intent(out) :: u_sr(:, :)
         !> Long-range real-space potential (npts, npair)
         real(wp), allocatable, intent(out) :: ur_lr(:, :)
         !> Long-range reciprocal potential (nk, npair)
         real(wp), allocatable, intent(out) :: uk_lr(:, :)
         !> Error handling
         type(error_type), allocatable, intent(out) :: error
         !> Optional solvent potential, second side
         class(moz_potential_type), intent(in), optional :: solvent
         !> Optional coupling of the model
         class(coupling_type), intent(inout), target, optional :: coupling
      end subroutine potential_compute_1d

      !> Ng-split UV tables on a volume grid, one column per solvent site
      !>
      !> - Solute positions from the latest update of this potential
      !> - `solvent` required: there are no 3D VV tables
      !> - Mixing rule and Ng split of this potential
      !>
      !> @param[in]     self      updated solute potential
      !> @param[in]     grid      volume grid with reciprocal points
      !> @param[out]    u_sr      short-range potential, Hartree (ngrid, nv)
      !> @param[out]    ur_lr     long-range real-space potential, Hartree (ngrid, nv)
      !> @param[out]    uk_lr     long-range reciprocal potential, Hartree bohr^3 (npts_k, nv)
      !> @param[out]    error     no solvent, not updated, host data on the solvent, coupling, category or input failure
      !> @param[in]     solvent   updated solvent potential; optional only for the generic, required here
      !> @param[in,out] coupling  optional coupling of the model; host-data term i reads scope i
      module subroutine potential_compute_3d(self, grid, u_sr, ur_lr, uk_lr, error, solvent, coupling)
         implicit none(type, external)
         !> Updated solute potential
         class(moz_potential_type), intent(in) :: self
         !> Volume grid with reciprocal points
         class(moist_math_grid_3d_type), intent(in) :: grid
         !> Short-range potential (ngrid, nv)
         real(wp), allocatable, intent(out) :: u_sr(:, :)
         !> Long-range real-space potential (ngrid, nv)
         real(wp), allocatable, intent(out) :: ur_lr(:, :)
         !> Long-range reciprocal potential (npts_k, nv)
         complex(wp), allocatable, intent(out) :: uk_lr(:, :)
         !> Error handling
         type(error_type), allocatable, intent(out) :: error
         !> Solvent potential
         class(moz_potential_type), intent(in), optional :: solvent
         !> Optional coupling of the model
         class(coupling_type), intent(inout), target, optional :: coupling
      end subroutine potential_compute_3d

      !> Reverse mode of the 1D UV tables
      !>
      !> - Convention: L = sum w_sr u_sr + sum w_lr ur_lr + sum w_k uk_lr
      !> - Kernel adjoints of the solute radial field and tail charges handed to every contributing term
      !> - No geometry on radial grids, so no gradient
      !>
      !> @param[in,out] self      updated solute potential
      !> @param[in]     rgrid     real-space radial grid
      !> @param[in]     kgrid     reciprocal radial grid
      !> @param[in]     solvent   updated solvent potential
      !> @param[in,out] coupling  coupling of the model; host-data term i reads scope i
      !> @param[in]     w_sr      adjoint of u_sr (npts, npair)
      !> @param[in]     w_lr      adjoint of ur_lr (npts, npair)
      !> @param[in]     w_k       adjoint of uk_lr (nk, npair)
      !> @param[in,out] response  host response items
      !> @param[out]    error     category, coupling, input or term adjoint failure
      module subroutine potential_compute_adjoint_1d(self, rgrid, kgrid, solvent, coupling, w_sr, w_lr, w_k, &
            & response, error)
         implicit none(type, external)
         !> Updated solute potential
         class(moz_potential_type), intent(inout) :: self
         !> Real-space radial grid
         type(moist_math_grid_radial_type), intent(in) :: rgrid
         !> Reciprocal radial grid
         type(moist_math_grid_radial_type), intent(in) :: kgrid
         !> Updated solvent potential
         class(moz_potential_type), intent(in) :: solvent
         !> Coupling of the model
         class(coupling_type), intent(inout), target :: coupling
         !> Adjoint of u_sr (npts, npair)
         real(wp), intent(in) :: w_sr(:, :)
         !> Adjoint of ur_lr (npts, npair)
         real(wp), intent(in) :: w_lr(:, :)
         !> Adjoint of uk_lr (nk, npair)
         real(wp), intent(in) :: w_k(:, :)
         !> Host response items
         type(response_type), intent(inout) :: response
         !> Error handling
         type(error_type), allocatable, intent(out) :: error
      end subroutine potential_compute_adjoint_1d

      !> Reverse mode of the 3D UV tables
      !>
      !> - Convention: L = sum w_sr u_sr + sum w_lr ur_lr + Re sum conj(w_k) uk_lr
      !> - Kernel adjoints of the solute field and tail charges handed to every contributing term
      !> - With `gradient` and `grid_adjoint`: geometric derivative of every piece built from the sites
      !>
      !> @param[in,out] self          updated solute potential
      !> @param[in]     grid          volume grid with reciprocal points
      !> @param[in]     solvent       updated solvent potential
      !> @param[in,out] coupling      coupling of the model; host-data term i reads scope i
      !> @param[in]     w_sr          adjoint of u_sr (ngrid, nv)
      !> @param[in]     w_lr          adjoint of ur_lr (ngrid, nv)
      !> @param[in]     w_k           adjoint of uk_lr (npts_k, nv)
      !> @param[in,out] response      host response items
      !> @param[out]    error         category, coupling, input or term adjoint failure
      !> @param[in,out] gradient      optional solute nuclear gradient to accumulate, Hartree/bohr (3, nu)
      !> @param[in,out] grid_adjoint  optional volume adjoint accumulator of the model; passed with `gradient`
      module subroutine potential_compute_adjoint_3d(self, grid, solvent, coupling, w_sr, w_lr, w_k, response, &
            & error, gradient, grid_adjoint)
         implicit none(type, external)
         !> Updated solute potential
         class(moz_potential_type), intent(inout) :: self
         !> Volume grid with reciprocal points
         class(moist_math_grid_3d_type), intent(in) :: grid
         !> Updated solvent potential
         class(moz_potential_type), intent(in) :: solvent
         !> Coupling of the model
         class(coupling_type), intent(inout), target :: coupling
         !> Adjoint of u_sr (ngrid, nv)
         real(wp), intent(in) :: w_sr(:, :)
         !> Adjoint of ur_lr (ngrid, nv)
         real(wp), intent(in) :: w_lr(:, :)
         !> Adjoint of uk_lr (npts_k, nv)
         complex(wp), intent(in) :: w_k(:, :)
         !> Host response items
         type(response_type), intent(inout) :: response
         !> Error handling
         type(error_type), allocatable, intent(out) :: error
         !> Optional solute nuclear gradient to accumulate (3, nu)
         real(wp), intent(inout), optional :: gradient(:, :)
         !> Optional volume adjoint accumulator of the model
         type(volume_adjoint_type), intent(inout), optional :: grid_adjoint
      end subroutine potential_compute_adjoint_3d

   end interface

end module moist_model_moz_potential
