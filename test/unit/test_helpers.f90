!> Shared support routines for the unit-test suite; not a test itself
!>
!> Provides:
!>   * `center_at_origin(mol)` - centroid-shift a structure_type
!>   * `get_test_structures(mols, n)` - sample `n` MB16-43 + Heavy28 +
!>                                         Amino20x4 + But14diol + UPU23
!>                                         records (20% each, no replacement)
!>   * `get_test_radii(mol, radii)` - CPCM-table radii
!>   * `get_test_points(mol, points, n)` - `n` deterministic random sampling
!>                                         points inside/near mol's box
!>   * `get_test_cavity_iswig(mol, cavity, error, ...)` - ready-to-use iSwiG surface
!>   * `build_test_cavity(mol, nleb, ctx, ...)` - iSwiG surface on a caller-owned
!>                                         context and COSMO radii
!>   * `fill_point_charge_potential(cavity, coupling, qat, mol)` - analytic point-charge
!>                                         potential trace answered on the coupling
!>   * `set_host_potential(coupling, phi)` - answer the potential request directly
!>   * `fill_missing_with_zeros(cavity, coupling)` - zero answers to the
!>                                         derivative outputs
!>   * `submit(coupling, name, values)` - answer the current request or stop
!>   * `copy_potential_adjoint(response, item)`, `copy_density(response, item)`,
!>     `copy_gostshyp_amplitude(response, item)` - copy of one response item,
!>                                         unallocated when absent; walks a
!>                                         full pass, so the cursor is rewound
!>   * `stage_point_charge_energy(error, component, cavity, qat, mol, coupling)`
!>                                       - bare-component coupling, energy phase
!>                                         answered with the point-charge trace
!>   * `stage_model_point_charge_energy(error, model, qat, mol, coupling)`
!>                                       - the same for a general model
!>   * `get_test_cross(mol)` - five-carbon cross with concave seams
!>   * `check_moist_error(error, err, context)` - moist error -> testdrive failure
!>   * `fd4_scalar(fpp, fp, fm, fmm, h)` - 4-point central FD formula
!>   * `fd4_offsets` - the matching stencil offsets, in units of h
!>   * `rel_deviation(a, b)` - |a - b| / (1 + |b|)
!>   * `fill_legacy_radii(mol, radii, error)` - legacy per-element radius table
!>   * `build_numbering_map(numbering, map)` - persistent grid numbering ->
!>                                             current array index
!>
!> No global Fortran RNG state is touched (self-contained LCG), so the
!> point and structure samplers are safe under parallel test execution
module test_helpers
   use moist_cavity_iswig, only: moist_cavity_iswig_parameters_type
   use, intrinsic :: iso_fortran_env, only: int64
   use mctc_env, only: wp
   use mctc_io, only: structure_type, new
   use mctc_io_convert, only: aatoau
   use mctc_env_error, only: moist_error_type => error_type
   use mstore, only: get_structure
   use mstore_data_record, only: record_type
   use mstore_mb16_43, only: get_mb16_43_records
   use mstore_heavy28, only: get_heavy28_records
   use mstore_amino20x4, only: get_amino20x4_records
   use mstore_but14diol, only: get_but14diol_records
   use mstore_upu23, only: get_upu23_records
   use moist_cavity_type, only: cavity_type
   use moist_model_type, only: solvation_model_component_type
   use moist_model_general, only: solvation_model_general
   use moist_cavity_iswig, only: cavity_type_iswig, new_cavity_iswig
   use moist_context, only: moist_context_type, new_context
   use moist_radii, only: default_cpcm_radii, radius_type, new_radii_custom_atoms, &
                          radius_type_static, new_cosmo_radii
   use moist_channels_coupling, only: coupling_type, coupling_view_type, &
      & point_potential_request_type, gaussian_potential_request_type, moist_phase_energy, &
      & moist_phase_gradient, coupling_arm, coupling_make_view, coupling_close_view
   use moist_channels_response, only: response_type, potential_adjoint_response_type, &
      density_response_type, gostshyp_amplitude_response_type
   use moist_data_radii_legacy, only: get_radius_func
   use testdrive, only: error_type, test_failed
   implicit none(type, external)
   private

   public :: component_view, submit, read_fixture_moments
   public :: copy_potential_adjoint, copy_density, copy_gostshyp_amplitude
   public :: center_at_origin
   public :: get_test_structures
   public :: get_test_radii
   public :: get_test_points
   public :: get_test_cavity_iswig
   public :: build_test_cavity
   public :: fill_point_charge_potential
   public :: fill_point_charge_field
   public :: set_host_potential
   public :: fill_missing_with_zeros
   public :: stage_point_charge_energy
   public :: stage_model_point_charge_energy
   public :: get_test_cross
   public :: fd4_scalar
   public :: fd4_offsets
   public :: rel_deviation
   public :: check_moist_error
   public :: fill_legacy_radii
   public :: build_numbering_map

   !> Default n for get_test_structures (must be a multiple of 5)
   integer, parameter :: default_n_structures = 5
   !> Default n for get_test_points
   integer, parameter :: default_n_points = 7
   !> Default Lebedev order for get_test_cavity_iswig
   integer, parameter :: default_nleb = 26

   !> Stencil offsets, in units of h, matching `fd4_scalar`'s argument order
   real(wp), parameter :: fd4_offsets(4) = [2.0_wp, 1.0_wp, -1.0_wp, -2.0_wp]

   !> The 5 mstore collections that get_test_structures samples from
   integer, parameter :: n_datasets = 5
   character(len=*), parameter :: datasets(n_datasets) = [character(len=10):: &
                                                          "MB16-43", "Heavy28", "Amino20x4", "But14diol", "UPU23"]

   !> Answer one output of the current request, failing immediately on a native error
   interface submit
      module procedure submit_1, submit_2, submit_3
   end interface submit

contains

   !> Translate `mol` so its arithmetic centroid sits at the origin
   !> Pure positional shift; atomic identities and ordering preserved
   !>
   !> @param[inout] mol  structure whose %xyz is shifted in place
   subroutine center_at_origin(mol)
      !> Structure whose Cartesian coordinates are shifted to put the
      !> arithmetic centroid at the origin
      type(structure_type), intent(inout) :: mol
      real(wp) :: centroid(3)
      integer :: iat

      centroid = 0.0_wp
      do iat = 1, mol%nat
         centroid = centroid + mol%xyz(:, iat)
      end do
      centroid = centroid/real(mol%nat, wp)
      do iat = 1, mol%nat
         mol%xyz(:, iat) = mol%xyz(:, iat) - centroid
      end do
   end subroutine center_at_origin

   !> Populate `structures` with `n` mstore records sampled evenly across
   !> 5 collections (MB16-43, Heavy28, Amino20x4, But14diol, UPU23). Each
   !> collection contributes 20% of the total via Fisher-Yates shuffle on
   !> its record list (no duplicates within a collection). Different `n`
   !> values produce different but reproducible samples; the LCG state is
   !> self-contained
   !>
   !> `n` must be a multiple of 5 and >= 5. Default = 10 (2 per set)
   !>
   !> @param[out] structures  allocated array of mstore-sourced structures
   !> @param[in]  n           optional total count; default 10
   subroutine get_test_structures(structures, n)
      !> Output array of mstore-sourced structures
      type(structure_type), allocatable, intent(out) :: structures(:)
      !> Optional total count (must be multiple of 5)
      integer, optional, intent(in) :: n

      integer :: total, per_set, ids, k
      integer(int64) :: rng
      type(record_type), allocatable :: records(:)
      integer, allocatable :: order(:)

      total = default_n_structures
      if (present(n)) total = n
      if (total < n_datasets .or. mod(total, n_datasets) /= 0) then
         error stop "get_test_structures: n must be a positive multiple of 5"
      end if
      per_set = total/n_datasets

      allocate (structures(total))
      !* LCG seed mixes the total count so different n values yield
      !* different draws while each n is fully reproducible
      rng = int(total, int64)*1009_int64 + 12345_int64

      k = 0
      do ids = 1, n_datasets
         call load_dataset(trim(datasets(ids)), records)
         call shuffle_indices(size(records), rng, order)
         do k = 1, per_set
            call get_structure(structures((ids - 1)*per_set + k), &
                               trim(datasets(ids)), &
                               trim(records(order(k))%id))
         end do
         deallocate (records, order)
      end do
   end subroutine get_test_structures

   !> Return CPCM-table radii for `mol` using the project's standard
   !> `default_cpcm_radii()` model
   !>
   !> @param[in]  mol    structure to look up radii for
   !> @param[out] radii  allocated to size mol%nat; filled with CPCM radii
   subroutine get_test_radii(mol, radii)
      !> Structure whose per-atom radii are to be filled
      type(structure_type), intent(in) :: mol
      !> Allocated on exit to mol%nat; filled with CPCM-table radii
      real(wp), allocatable, intent(out) :: radii(:)
      class(radius_type), allocatable :: model
      type(moist_error_type), allocatable :: err

      model = default_cpcm_radii()
      call model%update(mol, err)
      if (allocated(err)) error stop "get_test_radii: "//trim(err%message)
      allocate (radii, source=model%f0)
   end subroutine get_test_radii

   !> Generate `n` deterministic random sampling points inside and near
   !> `mol`. Points are drawn uniformly from the atom bounding box padded
   !> by +/-2 bohr in each axis; candidates within 0.2 bohr of any nucleus
   !> are rejected. A self-contained 64-bit LCG seeded from mol%nat keeps
   !> results reproducible without touching the global Fortran RNG state,
   !> so the routine is safe under parallel test execution
   !>
   !> @param[in]  mol     structure used to anchor the bounding box
   !> @param[out] points  (3, n) allocatable array of sampling points
   !> @param[in]  n       optional point count; default 12
   subroutine get_test_points(mol, points, n)
      !> Structure used to seed the bounding-box random sampler
      type(structure_type), intent(in) :: mol
      !> Allocated on exit to (3, n); filled with sampling points
      real(wp), allocatable, intent(out) :: points(:, :)
      !> Optional point count; default `default_n_points`
      integer, optional, intent(in) :: n

      real(wp), parameter :: pad = 2.0_wp
      real(wp), parameter :: min_dist = 2.0e-1_wp
      integer, parameter :: max_attempts = 10000

      integer :: total
      integer(int64) :: rng
      real(wp) :: box_min(3), box_max(3), span(3), cand(3)
      integer :: ax, attempt, accepted

      total = default_n_points
      if (present(n)) total = n
      if (total < 1) error stop "get_test_points: n must be >= 1"
      if (mol%nat < 1) error stop "get_test_points: empty mol"

      do ax = 1, 3
         box_min(ax) = minval(mol%xyz(ax, :)) - pad
         box_max(ax) = maxval(mol%xyz(ax, :)) + pad
      end do
      span = box_max - box_min

      allocate (points(3, total))
      !* Seed mixes mol%nat with a fixed offset so different-sized
      !* molecules get distinct draws while each size is reproducible
      rng = int(mol%nat, int64)*97_int64 + 12345_int64
      accepted = 0
      do attempt = 1, max_attempts
         do ax = 1, 3
            cand(ax) = box_min(ax) + lcg_uniform(rng)*span(ax)
         end do
         if (.not. far_from_all_atoms(cand, mol%xyz, min_dist)) cycle
         accepted = accepted + 1
         points(:, accepted) = cand
         if (accepted == total) return
      end do
      error stop "get_test_points: not enough valid points"
   end subroutine get_test_points

   !> Build an iSwiG surface for `mol`
   !>
   !> @param[in]  mol           Structure to wrap
   !> @param[out] cavity        Constructed iSwiG cavity
   !> @param[out] error         Error handling
   !> @param[in]  nleb          Lebedev order; default `default_nleb`
   !> @param[in]  radius_model  Radius model; default CPCM-table per-atom radii
   !> @param[in]  cut_f         Switching-factor cutoff; cavity default if absent
   subroutine get_test_cavity_iswig(mol, cavity, error, nleb, radius_model, cut_f)
      !> Structure to wrap
      type(structure_type), intent(in) :: mol
      !> Constructed iSwiG cavity
      type(cavity_type_iswig), intent(out) :: cavity
      !> Error handling
      type(moist_error_type), allocatable, intent(out) :: error
      !> Optional Lebedev order
      integer, intent(in), optional :: nleb
      !> Optional radius model; CPCM-table per-atom radii if absent
      class(radius_type), intent(in), optional :: radius_model
      !> Optional switching-factor cutoff;  grid points with `f <= cut_f` are
      !> dropped at construction. Raise it to keep only the exposed surface
      real(wp), intent(in), optional :: cut_f

      !> Resolved constructor parameters
      type(moist_cavity_iswig_parameters_type) :: param

      !> Per-atom radii of the default model
      real(wp), allocatable :: radii(:)
      !> Default radius model, built only when none was supplied
      class(radius_type), allocatable :: default_model
      !> Resolved Lebedev order
      integer :: num_leb
      !> Run context
      type(moist_context_type), pointer :: ctx

      num_leb = default_nleb
      if (present(nleb)) num_leb = nleb

      param%num_leb = num_leb
      if (present(cut_f)) param%cut_f = cut_f

      allocate (ctx)
      call new_context(ctx)
      if (present(radius_model)) then
         call new_cavity_iswig( &
            self=cavity, &
            ctx=ctx, &
            radius_model=radius_model, &
            error=error, &
            param=param)
      else
         call get_test_radii(mol, radii)
         call new_radii_custom_atoms(radii, default_model, error)
         if (allocated(error)) return
         call new_cavity_iswig( &
            self=cavity, &
            ctx=ctx, &
            radius_model=default_model, &
            error=error, &
            param=param)
      end if
      if (allocated(error)) return

      call cavity%update(mol, error=error)
   end subroutine get_test_cavity_iswig

   !> Build a COSMO-radii iSwiG test cavity for a given molecule and Lebedev
   !> grid size
   !>
   !> Unlike `get_test_cavity_iswig`, the run context is created and owned by
   !> the caller so that it outlives the cavity borrowing it, so several
   !> cavities can share one context, and so the same context can be handed to
   !> a solvation-model component
   !>
   !> @param[in]  mol          Molecular structure
   !> @param[in]  nleb         Lebedev grid size
   !> @param[in]  ctx          Run context owned by the caller, borrowed by the cavity
   !> @param[out] radius_model Radius model storage
   !> @param[out] cavity       Constructed cavity
   !> @param[out] error        Error handling
   subroutine build_test_cavity(mol, nleb, ctx, radius_model, cavity, error)

      !> Molecular structure
      type(structure_type), intent(in) :: mol

      !> Lebedev grid size
      integer, intent(in) :: nleb

      !> Run context owned by the caller, borrowed by the cavity
      type(moist_context_type), intent(in), target :: ctx

      !> Radius model storage
      type(radius_type_static), intent(out) :: radius_model

      !> Constructed cavity
      type(cavity_type_iswig), intent(out) :: cavity

      !> Error handling
      type(moist_error_type), allocatable, intent(out) :: error

      call new_cosmo_radii(radius_model)
      call new_cavity_iswig(cavity, ctx, radius_model=radius_model, error=error, &
         param=moist_cavity_iswig_parameters_type(num_leb=nleb))
      if (allocated(error)) return

      call cavity%update(mol, error=error)

   end subroutine build_test_cavity

   !> Point-charge potential trace on the cavity grid, answered on the
   !> coupling's potential requests
   !>
   !> Fills the Coulomb potential of the point charges `qat` at the atom
   !> positions of `mol`, evaluated on the cavity grid the coupling was
   !> prepared with, and answers every potential request the host walk
   !> visits with it. A Gaussian request receives the potential a Gaussian
   !> probe of width `xi0` sees,
   !>
   !>    phi(i) = sum_j q_j erf(xi_i |r_i - R_j|) / |r_i - R_j|
   !>
   !> a point request the bare `sum_j q_j / |r_i - R_j|`. Only the potential
   !> is answered; every other output of `coupling` is left as is. One full
   !> pass: call it again after every `prepare_energy`
   !>
   !> @param[in]    cavity   Cavity the coupling was prepared with
   !> @param[inout] coupling Coupling whose potential requests are answered
   !> @param[in]    qat      Atomic point charges (nat)
   !> @param[in]    mol      Structure supplying the atom positions
   subroutine fill_point_charge_potential(cavity, coupling, qat, mol)

      !> Cavity the coupling was prepared with
      class(cavity_type), intent(in) :: cavity

      !> Coupling whose potential requests are answered
      type(coupling_type), intent(inout), target :: coupling

      !> Atomic point charges
      real(wp), intent(in) :: qat(:)

      !> Structure supplying the atom positions
      type(structure_type), intent(in) :: mol

      do while (coupling%next())
         select type (request => coupling%request())
         type is (gaussian_potential_request_type)
            call submit(coupling, "phi", point_charge_potential(cavity, qat, mol, .true.))
         type is (point_potential_request_type)
            call submit(coupling, "phi", point_charge_potential(cavity, qat, mol, .false.))
         end select
      end do

   end subroutine fill_point_charge_potential

   !> Coulomb potential of point charges on the cavity grid, bare or Gaussian-probed
   !>
   !> @param[in] cavity   Cavity supplying the grid and, when probed, the widths
   !> @param[in] qat      Atomic point charges (nat)
   !> @param[in] mol      Structure supplying the atom positions
   !> @param[in] gaussian Whether to screen with the cavity's `xi0`
   function point_charge_potential(cavity, qat, mol, gaussian) result(phi)
      class(cavity_type), intent(in) :: cavity
      real(wp), intent(in) :: qat(:)
      type(structure_type), intent(in) :: mol
      logical, intent(in) :: gaussian
      !> Potential trace on the grid
      real(wp), allocatable :: phi(:)
      !> Grid point and atom indices
      integer :: i, j
      !> Separation vector and its length
      real(wp) :: r_vec(3), r_dist
      !> Separation below which the singular self-term is skipped
      real(wp), parameter :: min_dist = 1.0e-10_wp

      allocate (phi(cavity%ngrid), source=0.0_wp)
      do i = 1, cavity%ngrid
         do j = 1, mol%nat
            r_vec(:) = cavity%xyz(:, i) - mol%xyz(:, j)
            r_dist = sqrt(sum(r_vec**2))
            if (r_dist < min_dist) cycle
            if (gaussian) then
               phi(i) = phi(i) + qat(j)*erf(cavity%xi0(i)*r_dist)/r_dist
            else
               phi(i) = phi(i) + qat(j)/r_dist
            end if
         end do
      end do
   end function point_charge_potential

   !> Raw point-charge potential derivatives on the cavity grid, answered on
   !> the coupling's potential requests
   !>
   !> The raw derivatives of the potential `fill_point_charge_potential`
   !> answers, evaluated on the cavity grid: the gradient with respect to the
   !> grid point,
   !>
   !>    dphi_dr(:, i) = - sum_j qat_j (r_i - R_j) / |r_i - R_j|^3
   !>
   !> screened by `erf(x) - 2 x exp(-x^2)/sqrt(pi)`, `x = xi_i |r_i - R_j|`,
   !> for a Gaussian request, and then also the width derivative
   !> `dphi_dxi(i) = sum_j qat_j 2 exp(-x^2)/sqrt(pi)`. The surface charge
   !> is not folded in: the component contracts the raw derivatives itself,
   !> and the FD gradient tests pin the result. Sign: `grad_r (1/|r - R|) =
   !> -(r - R)/|r - R|^3`, so the gradient points from the grid point towards
   !> a positive charge. Only the derivatives are answered, on every
   !> potential request the host walk visits; every other output is left as is
   !> Call it after every `prepare_gradient` (or `prepare_response` of a
   !> field-dependent cavity), and after `fill_missing_with_zeros` if both
   !> are used
   !>
   !> @param[in]    cavity   Cavity the coupling was prepared with
   !> @param[inout] coupling Coupling whose position weight request is answered
   !> @param[in]    qat      Atomic point charges (nat)
   !> @param[in]    mol      Structure supplying the atom positions
   subroutine fill_point_charge_field(cavity, coupling, qat, mol)

      !> Cavity the coupling was prepared with
      class(cavity_type), intent(in) :: cavity

      !> Coupling whose derivative outputs are answered
      type(coupling_type), intent(inout), target :: coupling

      !> Atomic point charges
      real(wp), intent(in) :: qat(:)

      !> Structure supplying the atom positions
      type(structure_type), intent(in) :: mol

      !> Raw derivatives on the grid
      real(wp), allocatable :: w_xyz(:, :), w_xi(:)

      do while (coupling%next())
         select type (request => coupling%request())
         type is (gaussian_potential_request_type)
            call point_charge_field(cavity, qat, mol, .true., w_xyz, w_xi)
            call submit(coupling, "dphi_dr", w_xyz)
            call submit(coupling, "dphi_dxi", w_xi)
         type is (point_potential_request_type)
            call point_charge_field(cavity, qat, mol, .false., w_xyz, w_xi)
            call submit(coupling, "dphi_dr", w_xyz)
         end select
      end do

   end subroutine fill_point_charge_field

   !> Raw derivatives of the point-charge potential on the cavity grid
   !>
   !> @param[in]  cavity   Cavity supplying the grid and, when probed, the widths
   !> @param[in]  qat      Atomic point charges (nat)
   !> @param[in]  mol      Structure supplying the atom positions
   !> @param[in]  gaussian Whether to screen with the cavity's `xi0`
   !> @param[out] w_xyz    dphi/dr (3, ngrid)
   !> @param[out] w_xi     dphi/dxi (ngrid), zero for the bare potential
   subroutine point_charge_field(cavity, qat, mol, gaussian, w_xyz, w_xi)
      class(cavity_type), intent(in) :: cavity
      real(wp), intent(in) :: qat(:)
      type(structure_type), intent(in) :: mol
      logical, intent(in) :: gaussian
      real(wp), allocatable, intent(out) :: w_xyz(:, :), w_xi(:)
      !> Grid point and atom indices
      integer :: i, j
      !> Separation vector, its length, and the potential gradient
      real(wp) :: r_vec(3), r_dist, grad_phi(3), x, screening
      !> Separation below which the singular self-term is skipped
      real(wp), parameter :: min_dist = 1.0e-10_wp

      allocate (w_xyz(3, cavity%ngrid), w_xi(cavity%ngrid), source=0.0_wp)
      do i = 1, cavity%ngrid
         grad_phi(:) = 0.0_wp
         do j = 1, mol%nat
            r_vec(:) = cavity%xyz(:, i) - mol%xyz(:, j)
            r_dist = sqrt(sum(r_vec**2))
            if (r_dist < min_dist) cycle
            screening = 1.0_wp
            if (gaussian) then
               x = cavity%xi0(i)*r_dist
               screening = erf(x) - 2.0_wp*x*exp(-x*x)/sqrt(acos(-1.0_wp))
               w_xi(i) = w_xi(i) + qat(j)*2.0_wp*exp(-x*x)/sqrt(acos(-1.0_wp))
            end if
            grad_phi(:) = grad_phi(:) - qat(j)*screening*r_vec(:)/(r_dist*r_dist*r_dist)
         end do
         w_xyz(:, i) = grad_phi(:)
      end do
   end subroutine point_charge_field

   !> Answer the potential request of `coupling` with a host-computed array
   !>
   !> For tests whose potential is not the point-charge trace; answers every
   !> potential request the host walk visits and stops on a refused answer
   !>
   !> @param[inout] coupling Coupling whose potential request is answered
   !> @param[in]    phi      Potential at the cavity grid points (ngrid)
   subroutine set_host_potential(coupling, phi)

      !> Coupling whose potential request is answered
      type(coupling_type), intent(inout), target :: coupling

      !> Potential at the cavity grid points
      real(wp), intent(in) :: phi(:)

      do while (coupling%next())
         select type (request => coupling%request())
         type is (gaussian_potential_request_type)
            call submit(coupling, "phi", phi)
         type is (point_potential_request_type)
            call submit(coupling, "phi", phi)
         end select
      end do

   end subroutine set_host_potential

   !> Answer the derivative outputs of every potential request the host walk
   !> visits with zeros
   !>
   !> For tests whose subject is not the host quantity in question: a host
   !> whose potential is a fixed array has no surface response and supplies
   !> zero weights. Answering (rather than skipping) also satisfies the
   !> mandatory requests of the phase. The potential and the Gaussian moments
   !> are left missing: zeros there would hide a missing fixture. So is the
   !> surface position weight in the gradient phase: there it is the total
   !> `q_i grad phi(r_i)`, never zero for a point-charge potential, and a
   !> forgotten `fill_point_charge_field` must fail by name rather than pass
   !> a silently wrong gradient (the forward-versus-reverse tests would not
   !> notice)
   !>
   !> @param[in]    cavity   Cavity the coupling was prepared with
   !> @param[inout] coupling Staged coupling
   subroutine fill_missing_with_zeros(cavity, coupling)
      !> Cavity the coupling was prepared with
      class(cavity_type), intent(in) :: cavity
      !> Staged coupling
      type(coupling_type), intent(inout), target :: coupling
      !> Scoped view, only to learn the staged phase
      type(coupling_view_type) :: view
      !> Staged phase
      integer :: phase
      call coupling_make_view(coupling, 0, view)
      phase = view%phase
      call coupling_close_view(coupling)
      if (phase == moist_phase_gradient) return
      do while (coupling%next())
         select type (request => coupling%request())
         type is (gaussian_potential_request_type)
            call submit(coupling, "dphi_dr", spread(spread(0.0_wp, 1, 3), 2, cavity%ngrid))
            call submit(coupling, "dphi_dxi", spread(0.0_wp, 1, cavity%ngrid))
         type is (point_potential_request_type)
            call submit(coupling, "dphi_dr", spread(spread(0.0_wp, 1, 3), 2, cavity%ngrid))
         end select
      end do
   end subroutine fill_missing_with_zeros

   !> Build a bare component's coupling and answer its energy phase with the
   !> point-charge potential of `qat`
   !>
   !> `new_coupling`, `prepare_energy` and `fill_point_charge_potential` in
   !> one call. A later SCF iteration, or a cavity rebuilt in place, is staged
   !> by calling `prepare_energy` and `fill_point_charge_potential` again on
   !> the same coupling (the coupling is created once per component and
   !> survives updates). Library errors become testdrive failures
   !>
   !> @param[out]   error     testdrive failure
   !> @param[inout] component Component whose requests the coupling declares
   !> @param[in]    cavity    Updated cavity the component was updated on
   !> @param[in]    qat       Atomic point charges (nat)
   !> @param[in]    mol       Structure supplying the atom positions
   !> @param[out]   coupling  Coupling built and staged for the energy phase
   subroutine stage_point_charge_energy(error, component, cavity, qat, mol, coupling)

      !> testdrive failure
      type(error_type), allocatable, intent(out) :: error

      !> Component whose requests the coupling declares
      class(solvation_model_component_type), intent(inout) :: component

      !> Updated cavity
      class(cavity_type), intent(in) :: cavity

      !> Atomic point charges
      real(wp), intent(in) :: qat(:)

      !> Structure supplying the atom positions
      type(structure_type), intent(in) :: mol

      !> Coupling built and staged for the energy phase
      type(coupling_type), intent(out), target :: coupling

      !> moist error
      type(moist_error_type), allocatable :: err

      call component%new_coupling(cavity, coupling, err)
      if (allocated(err)) then
         call test_failed(error, "coupling setup failed: "//err%message)
         return
      end if
      call component%prepare_energy(cavity, coupling, err)
      if (allocated(err)) then
         call test_failed(error, "energy-phase staging failed: "//err%message)
         return
      end if
      call fill_point_charge_potential(cavity, coupling, qat, mol)

   end subroutine stage_point_charge_energy

   !> Build a general model's coupling and answer its energy phase with the
   !> point-charge potential of `qat`
   !>
   !> The model-level counterpart of `stage_point_charge_energy`; the model
   !> must have been updated
   !>
   !> @param[out]   error    testdrive failure
   !> @param[inout] model    Updated general model
   !> @param[in]    qat      Atomic point charges (nat)
   !> @param[in]    mol      Structure supplying the atom positions
   !> @param[out]   coupling Coupling built and staged for the energy phase
   subroutine stage_model_point_charge_energy(error, model, qat, mol, coupling)

      !> testdrive failure
      type(error_type), allocatable, intent(out) :: error

      !> Updated general model
      class(solvation_model_general), intent(inout) :: model

      !> Atomic point charges
      real(wp), intent(in) :: qat(:)

      !> Structure supplying the atom positions
      type(structure_type), intent(in) :: mol

      !> Coupling built and staged for the energy phase
      type(coupling_type), pointer, intent(out) :: coupling

      !> moist error
      type(moist_error_type), allocatable :: err

      call model%new_coupling(coupling, err)
      if (allocated(err)) then
         call test_failed(error, "coupling setup failed: "//err%message)
         return
      end if
      call model%prepare_energy(coupling, err)
      if (allocated(err)) then
         call test_failed(error, "energy-phase staging failed: "//err%message)
         return
      end if
      call fill_point_charge_potential(model%cavity, coupling, qat, mol)

   end subroutine stage_model_point_charge_energy

   !> Five-carbon cross, converted to bohr
   !>
   !> Deliberately concave seams between the four outer atoms: at the
   !> unconditional-multistart projection level this geometry produces branched
   !> anchors, which is what the warm-start and gradient FD suites need
   !>
   !> @param[out] mol  Five-carbon cross structure, coordinates in bohr
   subroutine get_test_cross(mol)
      !> Resulting structure; coordinates in bohr
      type(structure_type), intent(out) :: mol

      call new(mol, [6, 6, 6, 6, 6], reshape([ &
                                             0.00_wp, 4.21_wp, 0.00_wp, &
                                             0.00_wp, 0.00_wp, 4.22_wp, &
                                             0.00_wp, -4.18_wp, 0.00_wp, &
                                             0.00_wp, 0.00_wp, -4.15_wp, &
                                             0.02_wp, 0.10_wp, -0.20_wp], &
                                             [3, 5])*aatoau)
   end subroutine get_test_cross

   !> Turn a moist error into a testdrive test failure
   !>
   !> Safe to call unconditionally: an unallocated `err` is a no-op, so the
   !> idiom at every call site collapses to
   !> `call check_moist_error(error, err); if (allocated(error)) return`
   !> This is the one assertion in an otherwise fixture-only module; it lives
   !> here because the alternative is the same four lines copied into every
   !> suite that touches a moist routine
   !>
   !> @param[out] error    Test failure, allocated only when `err` was
   !> @param[in]  err      Moist error to translate
   !> @param[in]  context  Optional prefix, e.g. the operation that failed
   subroutine check_moist_error(error, err, context)
      !> Test failure; allocated only when `err` is
      type(error_type), allocatable, intent(out) :: error
      !> Moist error to translate; unallocated means success
      type(moist_error_type), allocatable, intent(in) :: err
      !> Optional prefix describing what was being attempted
      character(len=*), intent(in), optional :: context

      if (.not. allocated(err)) return
      if (present(context)) then
         call test_failed(error, context//": "//trim(err%message))
      else
         call test_failed(error, trim(err%message))
      end if
   end subroutine check_moist_error

   !> 4-point central finite-difference formula:
   !>   f'(x) ~ (-f(x+2h) + 8 f(x+h) - 8 f(x-h) + f(x-2h)) / (12 h)
   !> Truncation O(h^4 f^(5)); useful for FD-checking analytic derivatives
   pure real(wp) function fd4_scalar(fpp, fp, fm, fmm, h) result(df)
      !> Value at x + 2h
      real(wp), intent(in) :: fpp
      !> Value at x + h
      real(wp), intent(in) :: fp
      !> Value at x - h
      real(wp), intent(in) :: fm
      !> Value at x - 2h
      real(wp), intent(in) :: fmm
      !> Step size h
      real(wp), intent(in) :: h

      df = (-fpp + 8.0_wp*fp - 8.0_wp*fm + fmm)/(12.0_wp*h)
   end function fd4_scalar

   !> Deviation of `a` from reference `b`, relative but safe near zero
   !>   |a - b| / (1 + |b|)
   !>
   !> @param[in] a  Value under test
   !> @param[in] b  Reference value
   elemental pure real(wp) function rel_deviation(a, b) result(dev)
      !> Value under test
      real(wp), intent(in) :: a
      !> Reference value
      real(wp), intent(in) :: b

      dev = abs(a - b)/(1.0_wp + abs(b))
   end function rel_deviation

   !* ===================================================================
   !*                          Private helpers
   !* ===================================================================

   !> Dispatch to the per-dataset records getter. Caller frees `records`
   subroutine load_dataset(name, records)
      character(len=*), intent(in) :: name
      type(record_type), allocatable, intent(out) :: records(:)

      select case (name)
      case ("MB16-43")
         call get_mb16_43_records(records)
      case ("Heavy28")
         call get_heavy28_records(records)
      case ("Amino20x4")
         call get_amino20x4_records(records)
      case ("But14diol")
         call get_but14diol_records(records)
      case ("UPU23")
         call get_upu23_records(records)
      case default
         error stop "load_dataset: unknown collection '"//trim(name)//"'"
      end select
   end subroutine load_dataset

   !> Fisher-Yates: produce a random permutation of [1..n] using `rng`
   !> Allocates `order(n)` on output. Stateful in `rng`, no globals touched
   subroutine shuffle_indices(n, rng, order)
      integer, intent(in) :: n
      integer(int64), intent(inout) :: rng
      integer, allocatable, intent(out) :: order(:)
      integer :: i, j, tmp

      allocate (order(n))
      do i = 1, n
         order(i) = i
      end do
      !* Standard Fisher-Yates: for i = n downto 2, swap order(i) with
      !* order(j) where j uniform in [1, i]
      do i = n, 2, -1
         j = 1 + int(lcg_uniform(rng)*real(i, wp))
         if (j > i) j = i   ! guard against rounding to exactly i
         tmp = order(i)
         order(i) = order(j)
         order(j) = tmp
      end do
   end subroutine shuffle_indices

   !> True iff `point` is at least `min_dist` from every column of `centers`
   logical function far_from_all_atoms(point, centers, min_dist) result(ok)
      real(wp), intent(in) :: point(3)
      real(wp), intent(in) :: centers(:, :)
      real(wp), intent(in) :: min_dist
      real(wp) :: min_dist_sq
      integer :: iat

      min_dist_sq = min_dist*min_dist
      ok = .true.
      do iat = 1, size(centers, dim=2)
         if (sum((point - centers(:, iat))**2) < min_dist_sq) then
            ok = .false.
            return
         end if
      end do
   end function far_from_all_atoms

   !> 64-bit LCG (Knuth MMIX constants from TAOCP vol 2). Self-contained;
   !> does not touch the global Fortran RNG state, so safe under parallel
   !> test execution
   !>
   !> @param[inout] state  LCG state, advanced by one step
   real(wp) function lcg_uniform(state) result(u)
      integer(int64), intent(inout) :: state
      state = state*6364136223846793005_int64 + 1442695040888963407_int64
      !* Top 31 bits as a nonneg integer; divide by 2^31 to land in [0, 1)
      u = real(ishft(state, -33), wp)/real(2_int64**31, wp)
   end function lcg_uniform

   !> Fill per-atom radii from the legacy per-element table
   !>
   !> The CPCM-flavoured cavity tests want the same radii the legacy code used,
   !> which is not what `get_test_radii` returns
   !>
   !> @param[in]  mol   Structure whose per-atom radii are filled
   !> @param[out] radii Allocated on exit to `mol%nat`
   !> @param[out] error Test failure, allocated if a lookup fails
   subroutine fill_legacy_radii(mol, radii, error)
      !> Structure whose per-atom radii are filled
      type(structure_type), intent(in) :: mol
      !> Allocated on exit to mol%nat
      real(wp), allocatable, intent(out) :: radii(:)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      !> Lookup failure from the radius table
      type(moist_error_type), allocatable :: err
      !> Atom index
      integer :: iat

      allocate (radii(mol%nat))
      do iat = 1, mol%nat
         radii(iat) = get_radius_func(mol%num(mol%id(iat)), err)
         if (allocated(err)) then
            call test_failed(error, "radius lookup failed: "//trim(err%message))
            return
         end if
      end do
   end subroutine fill_legacy_radii

   !> Invert a cavity's persistent grid numbering into an index lookup
   !>
   !> DROP filters and reorders its grid on every rebuild, so a finite-difference
   !> comparison between two builds cannot be keyed on the array index. It has to
   !> go through `cavity%numbering`, which survives the rebuild. `map(n)` is the
   !> current array index of the point numbered `n`, or 0 if that point is absent
   !> from this build
   !>
   !> @param[in]  numbering Persistent point numbers, one per current grid point
   !> @param[out] map       Numbering -> index, sized to the largest number seen
   subroutine build_numbering_map(numbering, map)
      !> Persistent point numbers
      integer, intent(in) :: numbering(:)
      !> Numbering -> current index
      integer, allocatable, intent(out) :: map(:)

      !> Grid index and the number it carries
      integer :: igrid, inum
      !> Largest number present
      integer :: max_num

      max_num = 0
      if (size(numbering) > 0) max_num = maxval(numbering)

      allocate (map(max(1, max_num)), source=0)
      do igrid = 1, size(numbering)
         inum = numbering(igrid)
         if (inum > 0 .and. inum <= max_num) map(inum) = igrid
      end do
   end subroutine build_numbering_map

   !> Borrow a component-local view for a single test invocation
   !>
   !> A coupling nobody staged, e.g. the default one of a request-free
   !> component, is armed for the energy phase first
   function component_view(coupling) result(view)
      type(coupling_type), target, intent(inout) :: coupling
      type(coupling_view_type) :: view
      type(moist_error_type), allocatable :: err
      call coupling_make_view(coupling, 1, view)
      if (view%phase /= 0) return
      call coupling_arm(coupling, moist_phase_energy, err)
      call coupling_make_view(coupling, 1, view)
   end function component_view

   !> Read the four moment outputs of a fixture with exactly one moment calculation
   subroutine read_fixture_moments(coupling, gt, pt, mt, rt, err)
      type(coupling_type), intent(inout), target :: coupling
      real(wp), allocatable, intent(out) :: gt(:), pt(:, :), mt(:, :, :), rt(:, :)
      type(moist_error_type), allocatable, intent(out) :: err
      type(coupling_view_type) :: view
      call coupling_make_view(coupling, 1, view)
      call view%read("moments", "gt", gt, err)
      if (.not. allocated(err)) call view%read("moments", "pt", pt, err)
      if (.not. allocated(err)) call view%read("moments", "mt", mt, err)
      if (.not. allocated(err)) call view%read("moments", "rt", rt, err)
      call coupling_close_view(coupling)
   end subroutine read_fixture_moments

   !> Answer one vector output of the current request or stop
   subroutine submit_1(coupling, name, values)
      type(coupling_type), intent(inout) :: coupling
      character(len=*), intent(in) :: name
      real(wp), intent(in) :: values(:)
      type(moist_error_type), allocatable :: err
      call coupling%answer(name, values, err)
      if (allocated(err)) error stop err%message
   end subroutine submit_1

   !> Answer one (lead, ngrid) output of the current request or stop
   subroutine submit_2(coupling, name, values)
      type(coupling_type), intent(inout) :: coupling
      character(len=*), intent(in) :: name
      real(wp), intent(in) :: values(:, :)
      type(moist_error_type), allocatable :: err
      call coupling%answer(name, values, err)
      if (allocated(err)) error stop err%message
   end subroutine submit_2

   !> Answer one (lead(1), lead(2), ngrid) output of the current request or stop
   subroutine submit_3(coupling, name, values)
      type(coupling_type), intent(inout) :: coupling
      character(len=*), intent(in) :: name
      real(wp), intent(in) :: values(:, :, :)
      type(moist_error_type), allocatable :: err
      call coupling%answer(name, values, err)
      if (allocated(err)) error stop err%message
   end subroutine submit_3

   !> Copy of the potential adjoint item, unallocated when the response has none
   !>
   !> Walks the whole pass, so the response cursor is rewound on return
   subroutine copy_potential_adjoint(response, item)
      type(response_type), intent(inout) :: response
      type(potential_adjoint_response_type), allocatable, intent(out) :: item
      do while (response%next())
         select type (found => response%item())
         type is (potential_adjoint_response_type)
            item = found
         end select
      end do
   end subroutine copy_potential_adjoint

   !> Copy of the density item, unallocated when the response has none
   !>
   !> Walks the whole pass, so the response cursor is rewound on return
   subroutine copy_density(response, item)
      type(response_type), intent(inout) :: response
      type(density_response_type), allocatable, intent(out) :: item
      do while (response%next())
         select type (found => response%item())
         type is (density_response_type)
            item = found
         end select
      end do
   end subroutine copy_density

   !> Copy of the GOSTSHYP amplitude item, unallocated when the response has none
   !>
   !> Walks the whole pass, so the response cursor is rewound on return
   subroutine copy_gostshyp_amplitude(response, item)
      type(response_type), intent(inout) :: response
      type(gostshyp_amplitude_response_type), allocatable, intent(out) :: item
      do while (response%next())
         select type (found => response%item())
         type is (gostshyp_amplitude_response_type)
            item = found
         end select
      end do
   end subroutine copy_gostshyp_amplitude
end module test_helpers
