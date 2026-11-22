!> Molecular frame: axes that follow the atoms
!>
!> Provides
!>
!> - [[sym_frame_type]]: axes, level and the data of their derivative
!> - [[new_sym_frame]]: frame of a set of coordinates
!> - [[sym_frame_vjp]]: reverse-mode derivative of the frame
!> - [[sym_threshold_type]]: thresholds of the discrete choices
!> - `frame_level_*` codes with [[frame_level_name]] and [[frame_level_list]]
!>
!> Concepts
!>
!> - frame tensor: the second moment T = sum_a d_a d_a^T of d_a = r_a -
!>   centroid, uniform weights; mixed with a fourth moment where T has a
!>   degenerate pair (`sym_mix_type`)
!> - level: how the axes are fixed, tested in the order of the table
!>
!> | level    | when                                  | axis                         | in-plane angle     |
!> |----------|---------------------------------------|------------------------------|--------------------|
!> | `none`   | fewer than two distinct points        | identity frame               |                    |
!> | `axial`  | atoms on one line                     | the line                     | lab frame          |
!> | `eigen`  | three separated eigenvalues           | all three eigenvectors       |                    |
!> | `plane`  | one separated eigenvalue              | its eigenvector              | harmonic           |
!> | `axial`  | as `plane`, no harmonic large enough  | its eigenvector              | lab frame          |
!> | `sphere` | no separated eigenvalue               | extremum of the cubic or     | harmonic, else lab |
!> |          |                                       | quartic moment, else an atom |                    |
!>
!> - separated: relative gap to both neighbours above `gap` (`pair_gaps`)
!> - harmonic: c_m = sum_a z_a^m of the in-plane coordinates as complex
!>   numbers, m = 3, ..., 12; fixes the angle modulo 2 pi/m, the m candidates
!>   are ranked by the other harmonics (`plane_frame`)
!> - axis of the level `sphere`: `sphere_axis`
!> - axis signs and order are free: a Lebedev grid is invariant under every
!>   signed permutation of its axes, so only the axes as lines matter
!> - the axis is placed on a cube face axis or body diagonal of the grid,
!>   whichever leaves fewer candidates (`place_*`)
!> - the frame follows rotation, translation and reordering of the atoms,
!>   except for an angle taken from the lab frame and for exact ties, which
!>   arise under a symmetry of the atom positions (other frames in `tied`)
!> - every choice switches at a threshold, where the frame jumps under a
!>   small displacement; the `near_*` flags mark a frame close to one
module moist_math_symmetry_frame
   use mctc_env, only: wp
   use mctc_io_constants, only: pi
   use mctc_io_utils, only: to_string
   use moist_math_lapack_kinds, only: lapack_ik
   use moist_math_lapack_syev, only: lapack_syev
   use moist_math_linalg_geometry, only: cross_product, setup_tangent_frame
   implicit none(type, external)
   private

   public :: sym_frame_type, sym_frame_vjp
   public :: sym_threshold_type, new_sym_frame, near_below, near_above
   public :: frame_level_none, frame_level_eigen, frame_level_plane, frame_level_sphere, &
      & frame_level_axial, frame_level_name, frame_level_list

   !> Fewer than two distinct points: identity frame
   integer, parameter :: frame_level_none = 0
   !> Three separated eigenvalues of the frame tensor
   integer, parameter :: frame_level_eigen = 1
   !> One separated eigenvalue, in-plane angle from a harmonic
   integer, parameter :: frame_level_plane = 2
   !> No separated eigenvalue, axis from a higher moment or the farthest atom
   integer, parameter :: frame_level_sphere = 3
   !> Atoms on a line, or one separated eigenvalue without a harmonic; in-plane angle from the lab frame
   integer, parameter :: frame_level_axial = 4

   !> Relative size above which a harmonic or higher moment is significant
   real(wp), parameter :: harm_tol = 1.0e-6_wp
   !> Relative size below which a harmonic is taken as exactly zero
   real(wp), parameter :: zero_tol = 1.0e-12_wp
   !> Relative in-plane second moment below which the atoms lie on the axis
   !>
   !> - relative to the full second moment tr(T)
   !> - transverse coordinates of a line are roundoff, about 1e-16 of the
   !>   extent (1e-32 here; 1e-24 for a molecule 1e4 extents from the origin)
   !> - a bend of 1e-6 rad gives about 1e-12, so it stays `plane`
   real(wp), parameter :: line_tol = 1.0e-20_wp
   !> The same for the two smallest eigenvalues of T, relative to tr(T)
   !>
   !> - the eigenvalues of a line carry the rounding of the eigensolver,
   !>   about 1e-16 tr(T), wherever the line points
   real(wp), parameter :: line_eig_tol = 64.0_wp*epsilon(1.0_wp)
   !> Highest in-plane harmonic
   integer, parameter :: max_harmonic = 12
   !> Factor per order on the size of a harmonic when the strongest is chosen
   !>
   !> - a regular n-gon has full harmonics n, 2n, ...: the lowest one wins,
   !>   which has the fewest candidates
   real(wp), parameter :: order_bias = 0.9_wp
   !> Lower edge of the band around a threshold, as a factor on it
   !>
   !> - a measure inside the band is close enough for a small displacement
   !>   to cross the threshold; an exact degeneracy lies far below it
   real(wp), parameter :: near_below = 0.1_wp
   !> Upper edge of the band around a threshold, as a factor on it
   real(wp), parameter :: near_above = 100.0_wp
   !> Difference of two scores or harmonic sizes below which they tie exactly
   real(wp), parameter :: tie_zero = 1.0e-10_wp
   !> Atom directions tried as starting points of the moment maximization
   integer, parameter :: max_starts = 12
   !> Atoms whose pairs give more starting points for a moment minimum
   integer, parameter :: max_pair_starts = 5
   !> Most exactly tied frames kept in `tied`
   integer, parameter :: max_tied = 64
   !> Lab direction fixing the node numbering among equivalent frames
   !>
   !> - settles the sign of an eigenvector axis and candidates of one class
   !> - generic, so that an axis along a lab axis or a lab diagonal is
   !>   neither parallel nor normal to it; the sphere does not depend on it
   real(wp), parameter :: number_dir(3) = [0.36_wp, 0.48_wp, 0.80_wp]
   !> Axis along the third grid axis: the frame axes are (p, q, n)
   !>
   !> - candidates a quarter turn apart about the axis give the same sphere
   integer, parameter :: place_face = 0
   !> Axis along the body diagonal of the grid: the frame axes are (p, q, n) U^T
   !>
   !> - candidates a third of a turn apart give the same sphere
   !> - used where that leaves fewer distinct candidates than a face axis,
   !>   harmonics m = 3, 6, 9
   integer, parameter :: place_diagonal = 1
   !> As `place_diagonal` for an axis with a sense, turned by 30 degrees
   integer, parameter :: place_tetra = 2
   !> As `place_face`, turned by 45 degrees
   integer, parameter :: place_edge = 3
   !> Grid directions U of the diagonal placement as columns
   !>
   !> - third: the body diagonal (1, 1, 1)
   !> - first: the twofold axis (1, -1, 0), so that the half turn about p,
   !>   which a sign change of the axis amounts to, maps the grid onto itself
   real(wp), parameter :: diag_basis(3, 3) = reshape([ &
      & 1.0_wp/sqrt(2.0_wp), -1.0_wp/sqrt(2.0_wp), 0.0_wp, &
      & 1.0_wp/sqrt(6.0_wp), 1.0_wp/sqrt(6.0_wp), -2.0_wp/sqrt(6.0_wp), &
      & 1.0_wp/sqrt(3.0_wp), 1.0_wp/sqrt(3.0_wp), 1.0_wp/sqrt(3.0_wp)], [3, 3])
   !> Grid directions U of the tetrahedral placement as columns
   !>
   !> - first: the projection of the first grid axis, (2, -1, -1): atoms at
   !>   the azimuth of p below the plane lie on body diagonals, so a
   !>   tetrahedron meets the grid in its own cube whichever of its
   !>   threefold axes is the axis
   real(wp), parameter :: tetra_basis(3, 3) = reshape([ &
      & 2.0_wp/sqrt(6.0_wp), -1.0_wp/sqrt(6.0_wp), -1.0_wp/sqrt(6.0_wp), &
      & 0.0_wp, 1.0_wp/sqrt(2.0_wp), -1.0_wp/sqrt(2.0_wp), &
      & 1.0_wp/sqrt(3.0_wp), 1.0_wp/sqrt(3.0_wp), 1.0_wp/sqrt(3.0_wp)], [3, 3])

   !> Second moment mixed with the contracted fourth moment
   !>
   !> S = T/tr(T) + M/tr(M), M_ij = sum_klm Q_iklm Q_jklm,
   !> Q = sum_a d_a d_a d_a d_a
   !>
   !> - frame tensor where a pair of eigenvalues of T is not separated and
   !>   the atoms are no line
   !> - M has the axis of a symmetric top as an eigenvector also when T is
   !>   isotropic; only a cubic point group (or an accident) leaves both
   !>   isotropic
   type :: sym_mix_type
      !> T/tr(T)
      real(wp) :: that(3, 3) = 0.0_wp
      !> tr(T)
      real(wp) :: trace_t = 1.0_wp
      !> M/tr(M)
      real(wp) :: mhat(3, 3) = 0.0_wp
      !> tr(M)
      real(wp) :: trace_m = 1.0_wp
      !> Fourth moment Q
      real(wp) :: q(3, 3, 3, 3) = 0.0_wp
   end type sym_mix_type

   !> Grid directions U of the edge placement as columns
   !>
   !> - first: the face diagonal (1, 1, 0): atoms at the azimuth of p lie
   !>   above face diagonals, so a cube of atoms meets the grid in its own
   !>   axes whichever of its fourfold axes is the axis
   real(wp), parameter :: edge_basis(3, 3) = reshape([ &
      & 1.0_wp/sqrt(2.0_wp), 1.0_wp/sqrt(2.0_wp), 0.0_wp, &
      & -1.0_wp/sqrt(2.0_wp), 1.0_wp/sqrt(2.0_wp), 0.0_wp, &
      & 0.0_wp, 0.0_wp, 1.0_wp], [3, 3])

   !> Frame axes with the data its reverse-mode derivative needs
   type :: sym_frame_type
      !> Axes as columns, a proper rotation (3, 3)
      real(wp) :: axes(3, 3) = reshape([1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp, 0.0_wp, &
         & 0.0_wp, 0.0_wp, 1.0_wp], [3, 3])
      !> Resolution level, one of the `frame_level_*` values
      integer :: level = frame_level_none
      !> Smallest relative gap of neighbouring frame-tensor eigenvalues
      !>
      !> - relative to their sum; the frame quality, zero for a line
      real(wp) :: gap = 0.0_wp
      !> Eigenvalues of the frame tensor, ascending
      real(wp) :: lambda(3) = 0.0_wp
      !> Eigenvectors of the frame tensor as columns
      real(wp) :: evec(3, 3) = 0.0_wp
      !> Mixed frame tensor; unallocated when it is the second moment T
      type(sym_mix_type), allocatable :: mix
      !> Index of the separated eigenvalue (`plane`, `axial`)
      integer :: iu = 0
      !> In-plane harmonic fixing the angle, 0 when none
      integer :: harmonic = 0
      !> Size of that harmonic
      real(wp) :: c_abs = 0.0_wp
      !> Moment order fixing the `sphere` axis, 3 or 4; -1 for the farthest atom
      integer :: axis_order = 0
      !> Atom giving the `sphere` axis when `axis_order` is -1
      integer :: axis_atom = 0
      !> Whether the `sphere` axis is the minimum of the quartic moment
      logical :: axis_low = .false.
      !> Grid direction along the axis, one of the `place_*` values
      integer :: place = place_face
      !> Mean square distance of the atoms from the long axis, bohr^2
      real(wp) :: tau = 0.0_wp
      !> Whether a pair gap lies near the gap threshold
      logical :: near_gap = .false.
      !> Pair gap closest to the threshold when `near_gap`
      real(wp) :: near_gap_value = 0.0_wp
      !> Whether `tau` lies near the line threshold
      logical :: near_linear = .false.
      !> Whether a harmonic size lies near its threshold
      !>
      !> - also set when two harmonics or two candidates nearly tie
      logical :: near_harmonic = .false.
      !> Harmonic size or difference closest to its threshold when `near_harmonic`
      real(wp) :: near_harmonic_value = 0.0_wp
      !> Axes of exactly tied frames giving another sphere (3, 3, n)
      !>
      !> - unallocated without a tie
      !> - a tie arises under a symmetry of the atom positions, which do not
      !>   tell these apart; whether the molecule does is for the caller, who
      !>   knows the species
      real(wp), allocatable :: tied(:, :, :)
   end type sym_frame_type

   !> Thresholds at which the discrete choices of the frame switch
   !>
   !> | choice                        | measure                                   | against    |
   !> |-------------------------------|-------------------------------------------|------------|
   !> | line or not                   | (lambda_1 + lambda_2)/nat of T            | `linear`   |
   !> | `eigen`, `plane` or `sphere`  | pair gaps of the frame tensor             | `gap`      |
   !> | harmonic m or the lab angle   | relative size abs(c_m)/s_m                | `harmonic` |
   !>
   !> - choices without a threshold of their own (which harmonic, which of
   !>   its candidates, which `sphere` axis) take the largest value; see
   !>   `plane_frame` and `sphere_axis`
   type :: sym_threshold_type
      !> Relative gap of two frame-tensor eigenvalues separating the pair
      !>
      !> - relative to their sum; separated above it
      real(wp) :: gap = 1.0e-4_wp
      !> Relative size of an in-plane harmonic above which it fixes the angle
      real(wp) :: harmonic = 1.0e-4_wp
      !> Mean square distance from the long axis making a line, bohr^2
      !>
      !> - atoms below it form a line
      !> - absolute: only atoms within about 0.01 bohr of one line count,
      !>   whatever the length of the line
      real(wp) :: linear = 1.0e-4_wp
   end type sym_threshold_type

contains

   !> Molecular frame of a set of atoms
   !>
   !> @param[out] self Frame, its derivative data and the `near_*` flags
   !> @param[in]  xyz  Coordinates (3, nat)
   !> @param[in]  thr  Thresholds of the choices
   subroutine new_sym_frame(self, xyz, thr)
      !> Frame
      type(sym_frame_type), intent(out) :: self
      !> Coordinates
      real(wp), intent(in) :: xyz(:, :)
      !> Thresholds
      type(sym_threshold_type), intent(in) :: thr

      !> Coordinates about the centroid
      real(wp), allocatable :: d(:, :)
      real(wp) :: tmat(3, 3), trace, lambda(3), evec(3, 3), n(3), work(64), g(2)
      !> Axes that tie with the chosen one at the sphere level, and one of their frames
      real(wp), allocatable :: others(:, :)
      type(sym_frame_type) :: other
      integer(lapack_ik) :: info
      integer :: nat, k, j

      nat = size(xyz, 2)
      trace = 0.0_wp
      if (nat >= 2) then
         d = xyz - spread(sum(xyz, dim=2)/real(nat, wp), 2, nat)
         tmat = matmul(d, transpose(d))
         trace = tmat(1, 1) + tmat(2, 2) + tmat(3, 3)
      end if
      info = 1
      if (trace > 0.0_wp) then
         evec = tmat
         call lapack_syev("V", "U", 3_lapack_ik, evec, 3_lapack_ik, lambda, work, &
            & int(size(work), lapack_ik), info)
      end if
      ! No frame: the identity
      if (info /= 0) return

      self%lambda = lambda
      self%evec = evec
      self%tau = (lambda(1) + lambda(2))/real(nat, wp)
      self%near_linear = near_threshold(self%tau, thr%linear)

      ! Line: the transverse moments are rounding, or below the threshold
      if (lambda(1) + lambda(2) <= line_eig_tol*trace .or. self%tau < thr%linear) then
         self%level = frame_level_plane
         self%iu = 3
         call plane_frame(self, d, lab_sense(evec(:, 3)), evec(:, 1), thr, .false.)
         return
      end if

      ! Where a pair of the second moment is not separated the contracted
      ! fourth moment joins the frame tensor
      g = pair_gaps(lambda)
      call note_gap(self, g, thr)
      if (.not. all(g > thr%gap)) then
         allocate (self%mix)
         self%mix%that = tmat/trace
         self%mix%trace_t = trace
         call fourth_moment(d, self%mix%q, self%mix%mhat, self%mix%trace_m)
         evec = self%mix%that + self%mix%mhat
         call lapack_syev("V", "U", 3_lapack_ik, evec, 3_lapack_ik, lambda, work, &
            & int(size(work), lapack_ik), info)
         if (info /= 0) then
            ! Not reached for a symmetric 3x3 matrix: keep the second moment
            deallocate (self%mix)
            lambda = self%lambda
            evec = self%evec
         else
            self%lambda = lambda
            self%evec = evec
            g = pair_gaps(lambda)
            call note_gap(self, g, thr)
         end if
      end if
      self%gap = minval(g)

      if (all(g > thr%gap)) then
         ! Three separated eigenvalues
         self%level = frame_level_eigen
         self%axes = evec
         self%axes(:, 3) = cross_product(evec(:, 1), evec(:, 2))
      else if (g(2) > thr%gap) then
         ! Largest eigenvalue separated: its eigenvector is the axis
         self%level = frame_level_plane
         self%iu = 3
         call plane_frame(self, d, lab_sense(evec(:, 3)), evec(:, 1), thr, .true.)
      else if (g(1) > thr%gap) then
         ! Smallest eigenvalue separated
         self%level = frame_level_plane
         self%iu = 1
         call plane_frame(self, d, lab_sense(evec(:, 1)), evec(:, 2), thr, .true.)
      else
         ! Isotropic frame tensor
         self%level = frame_level_sphere
         call sphere_axis(self, d, n, others)
         call plane_frame(self, d, n, lab_reference(n), thr, .true.)
         ! The frames about the axes that tie with the chosen one
         do k = 1, size(others, 2)
            other = self
            if (allocated(other%tied)) deallocate (other%tied)
            call plane_frame(other, d, others(:, k), lab_reference(others(:, k)), thr, .true.)
            call add_tied(self, other%axes)
            if (.not. allocated(other%tied)) cycle
            do j = 1, size(other%tied, 3)
               call add_tied(self, other%tied(:, :, j))
            end do
         end do
      end if

   end subroutine new_sym_frame

   !> Frame about one axis: the in-plane angle from a harmonic or the lab
   !>
   !> - harmonic c_m = sum_a z_a^m, s_m = sum_a abs(z_a)^m, z_a the in-plane
   !>   coordinate of atom a as a complex number; c_m fixes the angle modulo
   !>   2 pi/m
   !> - m from 3 to `max_harmonic`: c_1 vanishes about the centroid, c_2 is
   !>   the second moment
   !> - harmonic used: largest `order_bias`^(m - 3) abs(c_m)/s_m among those
   !>   with abs(c_m)/s_m above `thr%harmonic`
   !> - candidate used, of the m angles: highest score of the other
   !>   harmonics (c_2 included); candidates giving the same sphere form one
   !>   class, see `place_face`, `place_diagonal`
   !> - classes that tie within `tie_zero`: on the diagonal, where classes a
   !>   sixth of a turn apart tie whenever the atoms come in pairs through a
   !>   center, the threefold harmonic odd in the height decides; then the
   !>   candidate closest to `number_dir`, so rounding does not pick; the
   !>   other tied classes go to `tied`
   !> - no harmonic: angle from the lab frame (`lab_reference` at the level
   !>   `plane`, which becomes `axial`; p0 at the level `sphere`)
   !>
   !> @param[in,out] self     Frame with level and axis data, receives the axes
   !> @param[in]     d        Coordinates about the centroid (3, nat)
   !> @param[in]     n        Axis, unit
   !> @param[in]     p0       Reference in-plane direction, unit and normal to n
   !> @param[in]     thr      Thresholds of the choices
   !> @param[in]     harmonic Whether the harmonics may fix the angle (no line)
   subroutine plane_frame(self, d, n, p0, thr, harmonic)
      !> Frame
      type(sym_frame_type), intent(inout) :: self
      !> Coordinates about the centroid
      real(wp), intent(in) :: d(:, :)
      !> Axis
      real(wp), intent(in) :: n(3)
      !> Reference in-plane direction
      real(wp), intent(in) :: p0(3)
      !> Thresholds
      type(sym_threshold_type), intent(in) :: thr
      !> Whether the harmonics may fix the angle
      logical, intent(in) :: harmonic

      real(wp) :: q0(3), s(max_harmonic), rel(max_harmonic), pref, best, second, p_lab(3)
      real(wp) :: phi0, phik, score, class_best(max_harmonic), class_side(max_harmonic)
      real(wp) :: height, s_side, top
      complex(wp) :: c(max_harmonic), z, c_side
      integer :: ia, m, mm, k, kbest, kclass, nclass, mbest, kalt
      logical :: class_tied(max_harmonic)

      q0 = cross_product(n, p0)
      c = (0.0_wp, 0.0_wp)
      s = 0.0_wp
      c_side = (0.0_wp, 0.0_wp)
      s_side = 0.0_wp
      do ia = 1, size(d, 2)
         z = cmplx(dot_product(d(:, ia), p0), dot_product(d(:, ia), q0), kind=wp)
         do m = 1, max_harmonic
            c(m) = c(m) + z**m
            s(m) = s(m) + powi(abs(z), m)
         end do
         ! Threefold harmonic odd in the height above the plane
         height = dot_product(d(:, ia), n)
         c_side = c_side + height*z**3
         s_side = s_side + abs(height)*powi(abs(z), 3)
      end do

      ! Strongest harmonic above the threshold, lower orders preferred;
      ! atoms on the axis carry rounding in place of harmonics
      mbest = 0
      if (harmonic .and. s(2) > line_tol*sum(d**2)) then
         best = 0.0_wp
         second = 0.0_wp
         rel = 0.0_wp
         do m = 3, max_harmonic
            if (.not. s(m) > 0.0_wp) cycle
            if (abs(c(m)) > zero_tol*s(m)) rel(m) = abs(c(m))/s(m)
            if (.not. rel(m) > thr%harmonic) cycle
            pref = order_bias**(m - 3)*rel(m)
            if (pref > best) then
               second = best
               best = pref
               mbest = m
            else if (pref > second) then
               second = pref
            end if
         end do
         ! Only the size that decides counts: the chosen harmonic, or the
         ! largest one when none passes
         if (mbest > 0) then
            call note_harmonic(self, rel(mbest), thr%harmonic)
            ! A lone harmonic has no rival to tie with
            if (second > 0.0_wp) call note_tie(self, best - second, thr%harmonic)
         else
            call note_harmonic(self, maxval(rel(3:)), thr%harmonic)
         end if
      end if

      if (mbest > 0) then
         m = mbest
         ! Candidates a quarter turn apart about a face axis, or a third of
         ! a turn about the diagonal, give the same sphere
         self%place = place_face
         nclass = m/gcd(m, 4)
         if (m/gcd(m, 3) < nclass) then
            self%place = place_diagonal
            nclass = m/gcd(m, 3)
            ! The maximum of the cubic moment gives the axis a sense
            if (self%axis_order == 3) self%place = place_tetra
         else if (self%axis_low) then
            ! Atoms lie between the fourfold axes of a quartic minimum
            self%place = place_edge
         end if

         ! Candidates phi0 + 2 pi k/m by the remaining harmonics
         phi0 = atan2(aimag(c(m)), real(c(m), wp))/real(m, wp)
         class_best(1:nclass) = -huge(1.0_wp)
         do k = 0, m - 1
            phik = phi0 + 2.0_wp*pi*real(k, wp)/real(m, wp)
            score = 0.0_wp
            do mm = 2, max_harmonic
               if (mm == m .or. .not. s(mm) > 0.0_wp) cycle
               score = score + (real(c(mm), wp)*cos(real(mm, wp)*phik) &
                  & + aimag(c(mm))*sin(real(mm, wp)*phik))/s(mm)
            end do
            class_best(mod(k, nclass) + 1) = max(class_best(mod(k, nclass) + 1), score)
         end do
         ! Another class gives another sphere
         top = maxval(class_best(1:nclass))
         do k = 1, nclass
            call note_tie(self, top - class_best(k), thr%harmonic)
         end do
         class_tied(1:nclass) = class_best(1:nclass) >= top - tie_zero
         ! Classes on the diagonal a sixth of a turn apart are mirror images
         ! in the plane, which the in-plane harmonics cannot tell apart (atoms
         ! in pairs through a center): the threefold harmonic odd in the
         ! height can; it does not depend on the sense of the axis
         if (count(class_tied(1:nclass)) > 1 .and. self%place /= place_face .and. self%place /= place_edge &
            & .and. s_side > 0.0_wp) then
            do k = 1, nclass
               if (.not. class_tied(k)) cycle
               z = c_side*exp(cmplx(0.0_wp, -3.0_wp*(phi0 + 2.0_wp*pi*real(k - 1, wp)/real(m, wp)), kind=wp))
               class_side(k) = merge(real(z, wp), aimag(z), self%place == place_tetra)/s_side
            end do
            top = maxval(class_side(1:nclass), mask=class_tied(1:nclass))
            do k = 1, nclass
               if (.not. class_tied(k)) cycle
               call note_tie(self, top - class_side(k), thr%harmonic)
               class_tied(k) = class_side(k) >= top - tie_zero
            end do
         end if
         ! Among classes that still tie, and among the candidates of one
         ! class (the same sphere with another numbering of its nodes): the
         ! candidate closest to a fixed lab direction. Rounding then does not
         ! pick, and the choice does not follow the arbitrary p0
         p_lab = number_dir - n*dot_product(n, number_dir)
         kbest = -1
         best = -huge(1.0_wp)
         do k = 0, m - 1
            if (.not. class_tied(mod(k, nclass) + 1)) cycle
            phik = phi0 + 2.0_wp*pi*real(k, wp)/real(m, wp)
            score = dot_product(cos(phik)*p0 + sin(phik)*q0, p_lab)
            if (score > best) then
               best = score
               kbest = k
            end if
         end do
         kclass = mod(kbest, nclass) + 1

         self%harmonic = m
         self%c_abs = abs(c(m))
         self%axes = placed_axes(kbest)
         ! The other classes of the tie
         do kalt = 1, nclass
            if (kalt == kclass .or. .not. class_tied(kalt)) cycle
            call add_tied(self, placed_axes(kalt - 1))
         end do
         return
      end if

      ! No harmonic: the in-plane angle of the lab frame
      self%harmonic = 0
      self%place = place_face
      self%axes(:, 3) = n
      if (self%level == frame_level_plane) then
         self%level = frame_level_axial
         p_lab = lab_reference(n)
         self%axes(:, 1) = p_lab
      else
         self%axes(:, 1) = p0
      end if
      self%axes(:, 2) = cross_product(self%axes(:, 3), self%axes(:, 1))

   contains

      !> Frame axes of one candidate of the chosen harmonic
      !>
      !> @param[in] kc Candidate, from 0
      pure function placed_axes(kc) result(axes)
         !> Candidate
         integer, intent(in) :: kc
         !> Frame axes
         real(wp) :: axes(3, 3)

         real(wp) :: pqn(3, 3), phi

         phi = phi0 + 2.0_wp*pi*real(kc, wp)/real(m, wp)
         pqn(:, 3) = n
         pqn(:, 1) = cos(phi)*p0 + sin(phi)*q0
         pqn(:, 2) = cross_product(n, pqn(:, 1))
         select case (self%place)
         case (place_diagonal)
            axes = matmul(pqn, transpose(diag_basis))
         case (place_tetra)
            axes = matmul(pqn, transpose(tetra_basis))
         case (place_edge)
            axes = matmul(pqn, transpose(edge_basis))
         case default
            axes = pqn
         end select
      end function placed_axes

   end subroutine plane_frame

   !> Append frames to the exactly tied ones of a frame, up to `max_tied`
   !>
   !> @param[in,out] self Frame
   !> @param[in]     axes Frame axes (3, 3)
   pure subroutine add_tied(self, axes)
      !> Frame
      type(sym_frame_type), intent(inout) :: self
      !> Frame axes
      real(wp), intent(in) :: axes(3, 3)

      real(wp), allocatable :: grown(:, :, :)
      integer :: n

      n = 0
      if (allocated(self%tied)) n = size(self%tied, 3)
      if (n >= max_tied) return
      allocate (grown(3, 3, n + 1))
      if (n > 0) grown(:, :, 1:n) = self%tied
      grown(:, :, n + 1) = axes
      call move_alloc(grown, self%tied)
   end subroutine add_tied

   !> Eigenvector with a sense fixed in the lab frame
   !>
   !> - the sign of an eigenvector is arbitrary and both give the same sphere
   !>   with another numbering of its nodes: the side of a fixed lab
   !>   direction keeps the numbering under a small displacement
   !>
   !> @param[in] v Unit vector
   pure function lab_sense(v) result(n)
      !> Unit vector
      real(wp), intent(in) :: v(3)
      !> v or -v
      real(wp) :: n(3)

      n = v
      if (dot_product(v, number_dir) < 0.0_wp) n = -v
   end function lab_sense

   !> Pair gaps of three ascending eigenvalues
   !>
   !> - g_ij = (lambda_j - lambda_i)/(lambda_i + lambda_j)
   !> - equals abs(c_2)/s_2 about the third eigenvector: the in-plane
   !>   quadrupole on the scale of a harmonic size; the length of the
   !>   molecule along that axis does not enter
   !>
   !> @param[in] lambda Eigenvalues, ascending, the largest positive
   pure function pair_gaps(lambda) result(g)
      !> Eigenvalues
      real(wp), intent(in) :: lambda(3)
      !> Gaps of the pairs (1, 2) and (2, 3)
      real(wp) :: g(2)

      g = 0.0_wp
      if (lambda(1) + lambda(2) > 0.0_wp) g(1) = (lambda(2) - lambda(1))/(lambda(1) + lambda(2))
      g(2) = (lambda(3) - lambda(2))/(lambda(2) + lambda(3))
   end function pair_gaps

   !> Whether a measure lies in the band around its threshold
   !>
   !> @param[in] x   Measure
   !> @param[in] thr Threshold
   pure function near_threshold(x, thr) result(near)
      !> Measure
      real(wp), intent(in) :: x
      !> Threshold
      real(wp), intent(in) :: thr
      !> Whether the measure is near
      logical :: near

      near = x > near_below*thr .and. x < near_above*thr
   end function near_threshold

   !> Flag pair gaps near the gap threshold
   !>
   !> @param[in,out] self Frame
   !> @param[in]     g    Pair gaps
   !> @param[in]     thr  Thresholds
   pure subroutine note_gap(self, g, thr)
      !> Frame
      type(sym_frame_type), intent(inout) :: self
      !> Pair gaps
      real(wp), intent(in) :: g(2)
      !> Thresholds
      type(sym_threshold_type), intent(in) :: thr

      integer :: k

      do k = 1, 2
         if (.not. near_threshold(g(k), thr%gap)) cycle
         if (self%near_gap) then
            if (abs(log(g(k)/thr%gap)) >= abs(log(self%near_gap_value/thr%gap))) cycle
         end if
         self%near_gap = .true.
         self%near_gap_value = g(k)
      end do
   end subroutine note_gap

   !> Flag a harmonic size near its threshold
   !>
   !> @param[in,out] self Frame
   !> @param[in]     rel  Relative size of the harmonic
   !> @param[in]     thr  Harmonic threshold
   pure subroutine note_harmonic(self, rel, thr)
      !> Frame
      type(sym_frame_type), intent(inout) :: self
      !> Relative size
      real(wp), intent(in) :: rel
      !> Threshold
      real(wp), intent(in) :: thr

      if (.not. near_threshold(rel, thr)) return
      if (self%near_harmonic) then
         if (abs(log(rel/thr)) >= abs(log(self%near_harmonic_value/thr))) return
      end if
      self%near_harmonic = .true.
      self%near_harmonic_value = rel
   end subroutine note_harmonic

   !> Flag two harmonics or two candidates that nearly tie
   !>
   !> - an exact tie, below `tie_zero`, is a symmetry of the atoms
   !>
   !> @param[in,out] self Frame
   !> @param[in]     diff Lead of the chosen one, not negative
   !> @param[in]     thr  Harmonic threshold
   pure subroutine note_tie(self, diff, thr)
      !> Frame
      type(sym_frame_type), intent(inout) :: self
      !> Lead of the chosen one
      real(wp), intent(in) :: diff
      !> Threshold
      real(wp), intent(in) :: thr

      if (.not. (diff > tie_zero .and. diff < near_above*thr)) return
      if (self%near_harmonic .and. diff >= self%near_harmonic_value) return
      self%near_harmonic = .true.
      self%near_harmonic_value = diff
   end subroutine note_tie

   !> Fourth moment Q = sum_a d_a d_a d_a d_a and its contraction
   !>
   !> M_ij = sum_klm Q_iklm Q_jklm follows the atoms like the second moment
   !> and is isotropic only for a cubic point group or by accident
   !>
   !> @param[in]  d     Coordinates about the centroid (3, nat)
   !> @param[out] q     Fourth moment
   !> @param[out] mhat  M/tr(M)
   !> @param[out] trace tr(M)
   pure subroutine fourth_moment(d, q, mhat, trace)
      !> Coordinates about the centroid
      real(wp), intent(in) :: d(:, :)
      !> Fourth moment
      real(wp), intent(out) :: q(3, 3, 3, 3)
      !> Normalized contraction
      real(wp), intent(out) :: mhat(3, 3)
      !> Trace of the contraction
      real(wp), intent(out) :: trace

      real(wp) :: o(3, 3)
      integer :: ia, i, j, k, l

      q = 0.0_wp
      do ia = 1, size(d, 2)
         do j = 1, 3
            o(:, j) = d(:, ia)*d(j, ia)
         end do
         do l = 1, 3
            do k = 1, 3
               q(:, :, k, l) = q(:, :, k, l) + o*o(k, l)
            end do
         end do
      end do
      do j = 1, 3
         do i = 1, 3
            mhat(i, j) = sum(q(:, :, :, i)*q(:, :, :, j))
         end do
      end do
      trace = mhat(1, 1) + mhat(2, 2) + mhat(3, 3)
      mhat = mhat/trace
   end subroutine fourth_moment

   !> Reverse-mode derivative of the frame tensor
   !>
   !> Given the symmetric adjoint of the frame tensor S, accumulates dE/dd
   !>
   !> - second moment, S = T: 2 sbar d_a
   !> - mixed, S = T/tr(T) + M/tr(M): through T and through Q
   !>
   !> @param[in]     d    Coordinates about the centroid (3, nat)
   !> @param[in]     sbar Adjoint of S, symmetric
   !> @param[in,out] grad Gradient, accumulated
   !> @param[in]     mix  Mixed frame tensor, absent for the second moment
   pure subroutine tensor_vjp(d, sbar, grad, mix)
      !> Coordinates about the centroid
      real(wp), intent(in) :: d(:, :)
      !> Adjoint of the frame tensor
      real(wp), intent(in) :: sbar(3, 3)
      !> Gradient
      real(wp), intent(inout) :: grad(:, :)
      !> Mixed frame tensor
      type(sym_mix_type), intent(in), optional :: mix

      real(wp) :: tbar(3, 3), mbar(3, 3), qbar(3, 3, 3, 3), sb, st, x(3), g(3), o(3, 3)
      integer :: ia, i, j, k

      if (all(sbar == 0.0_wp)) return
      if (.not. present(mix)) then
         do ia = 1, size(d, 2)
            grad(:, ia) = grad(:, ia) + 2.0_wp*matmul(sbar, d(:, ia))
         end do
         return
      end if

      st = sum(sbar*mix%that)
      sb = sum(sbar*mix%mhat)
      tbar = sbar
      mbar = sbar
      do i = 1, 3
         tbar(i, i) = tbar(i, i) - st
         mbar(i, i) = mbar(i, i) - sb
      end do
      tbar = tbar/mix%trace_t
      mbar = mbar/mix%trace_m
      ! dE = sum mbar_ij dM_ij = sum qbar_iklm dQ_iklm
      do i = 1, 3
         qbar(i, :, :, :) = 0.0_wp
         do j = 1, 3
            qbar(i, :, :, :) = qbar(i, :, :, :) + 2.0_wp*mbar(i, j)*mix%q(j, :, :, :)
         end do
      end do
      do ia = 1, size(d, 2)
         x = d(:, ia)
         do j = 1, 3
            o(:, j) = x*x(j)
         end do
         do k = 1, 3
            ! First slot, and the three equivalent other slots
            g(k) = sum(qbar(k, :, :, :)*spread(o, 1, 3)*spread(spread(x, 2, 3), 3, 3))
            do i = 1, 3
               g(k) = g(k) + 3.0_wp*x(i)*sum(qbar(i, k, :, :)*o)
            end do
         end do
         grad(:, ia) = grad(:, ia) + 2.0_wp*matmul(tbar, x) + g
      end do
   end subroutine tensor_vjp

   !> Integer power by binary powering, in the operation order of a runtime `x**n`
   !>
   !> - the runtime power divides by its result whatever the sign of n: a
   !>   zero base traps there under `-ffpe-trap=zero`
   !>
   !> @param[in] x Base
   !> @param[in] n Exponent, n >= 0
   elemental function powi(x, n) result(p)
      !> Base
      real(wp), intent(in) :: x
      !> Exponent
      integer, intent(in) :: n
      !> Power
      real(wp) :: p

      real(wp) :: b
      integer :: k

      p = 1.0_wp
      b = x
      k = n
      do
         if (iand(k, 1) /= 0) p = p*b
         k = k/2
         if (k == 0) exit
         b = b*b
      end do
   end function powi

   !> Greatest common divisor of two positive integers
   !>
   !> @param[in] a First integer
   !> @param[in] b Second integer
   pure function gcd(a, b) result(g)
      !> First integer
      integer, intent(in) :: a
      !> Second integer
      integer, intent(in) :: b
      !> Greatest common divisor
      integer :: g

      integer :: r, t

      g = a
      r = b
      do while (r /= 0)
         t = mod(g, r)
         g = r
         r = t
      end do
   end function gcd

   !> Axis of an isotropic frame tensor: an extremum of a higher moment
   !>
   !> - cubic moment f3(u) = sum_a (d_a . u)^3 when its maximum is significant
   !>   (above `harm_tol` relative); tried first
   !> - else the quartic moment f4, when its anisotropy is significant: its
   !>   maximum or its minimum, whichever lies farther from the isotropic
   !>   value s4/5; both are fourfold axes for an octahedron and a cube
   !> - else the direction of the atom farthest from the centroid
   !>
   !> @param[in,out] self Frame, receives `axis_order` and `axis_atom`
   !> @param[in]     d    Coordinates about the centroid (3, nat)
   !> @param[out]    u      Axis
   !> @param[out]    others Axes that tie exactly with it (3, n)
   subroutine sphere_axis(self, d, u, others)
      !> Frame
      type(sym_frame_type), intent(inout) :: self
      !> Coordinates about the centroid
      real(wp), intent(in) :: d(:, :)
      !> Axis
      real(wp), intent(out) :: u(3)
      !> Axes that tie exactly with it (3, n)
      real(wp), allocatable, intent(out) :: others(:, :)

      !> Radial distances and their moments
      real(wp), allocatable :: r(:)
      real(wp), allocatable :: low_others(:, :)
      real(wp) :: s3, s4, fmax, fmin, ulow(3), best, score
      integer :: k, nother

      r = norm2(d, dim=1)
      s3 = sum(r*r*r)
      s4 = sum((r*r)**2)

      call maximize_moment(d, r, 3, .false., u, fmax, others)
      if (fmax > harm_tol*s3) then
         self%axis_order = 3
         return
      end if
      call maximize_moment(d, r, 4, .false., u, fmax, others)
      call maximize_moment(d, r, 4, .true., ulow, fmin, low_others)
      ! The isotropic part of f4 is s4/5 in every direction
      if (max(fmax - s4/5.0_wp, s4/5.0_wp - fmin) > harm_tol*s4) then
         self%axis_order = 4
         if (s4/5.0_wp - fmin > fmax - s4/5.0_wp) then
            self%axis_low = .true.
            u = ulow
            call move_alloc(low_others, others)
         end if
         return
      end if

      ! A farthest atom; among equally far ones that closest to a fixed lab
      ! direction, which neither rounding nor the numbering decides
      self%axis_order = -1
      best = -huge(1.0_wp)
      do k = 1, size(r)
         if (r(k) < maxval(r)*(1.0_wp - tie_zero)) cycle
         score = dot_product(d(:, k), number_dir)
         if (score > best) then
            best = score
            self%axis_atom = k
         end if
      end do
      u = d(:, self%axis_atom)/r(self%axis_atom)
      if (allocated(others)) deallocate (others)
      allocate (others(3, min(count(r >= maxval(r)*(1.0_wp - tie_zero)) - 1, max_tied)))
      nother = 0
      do k = 1, size(r)
         if (k == self%axis_atom .or. r(k) < maxval(r)*(1.0_wp - tie_zero)) cycle
         if (nother == size(others, 2)) exit
         nother = nother + 1
         others(:, nother) = d(:, k)/r(k)
      end do
   end subroutine sphere_axis

   !> Global maximum or minimum of f(u) = sum_a (d_a . u)^m on the unit sphere
   !>
   !> - local ascents (descents) from the directions of the `max_starts`
   !>   farthest atoms and, for an odd order or a minimum, more directions
   !>   between them: shifted power steps, then Newton on the sphere
   !> - the best local extremum wins; among exactly tied ones that closest
   !>   to a fixed lab direction, which neither rounding nor the numbering
   !>   of the atoms decides
   !>
   !> @param[in]  d      Coordinates about the centroid (3, nat)
   !> @param[in]  r      Their lengths (nat)
   !> @param[in]  m      Moment order, 3 or 4
   !> @param[in]  lowest Whether the minimum is searched
   !> @param[out] u      Extremal direction
   !> @param[out] fmax   Extremum
   !> @param[out] others Directions that tie exactly with it (3, n)
   subroutine maximize_moment(d, r, m, lowest, u, fmax, others)
      !> Coordinates about the centroid
      real(wp), intent(in) :: d(:, :)
      !> Their lengths
      real(wp), intent(in) :: r(:)
      !> Moment order
      integer, intent(in) :: m
      !> Whether the minimum is searched
      logical, intent(in) :: lowest
      !> Extremal direction
      real(wp), intent(out) :: u(3)
      !> Extremum
      real(wp), intent(out) :: fmax
      !> Tied directions
      real(wp), allocatable, intent(out) :: others(:, :)

      !> Local extrema found and their values
      real(wp), allocatable :: found(:, :), value(:)
      real(wp) :: best, score
      integer :: nfound, ibest
      logical :: known
      !> Atoms by descending distance
      integer, allocatable :: order(:)
      !> Starting directions
      real(wp), allocatable :: starts(:, :)
      !> Candidate, its gradient, Hessian and value
      real(wp) :: v(3), g(3), h(3, 3), f, alpha, sgn
      !> Tangent basis and the projected Newton system
      real(wp) :: b1(3), b2(3), a(2, 2), rhs(2), det, step(2)
      logical, allocatable :: used(:)
      integer :: nat, istart, iter, k, j, nstart, npair, ns

      nat = size(d, 2)
      sgn = merge(-1.0_wp, 1.0_wp, lowest)
      fmax = -huge(1.0_wp)
      u = [0.0_wp, 0.0_wp, 1.0_wp]
      alpha = real(m*(m - 1), wp)*sum(powi(r, m))

      allocate (order(nat), used(nat))
      used = .false.
      nstart = min(max_starts, nat)
      do k = 1, nstart
         order(k) = maxloc(r, dim=1, mask=.not. used)
         used(order(k)) = .true.
      end do

      ! Atom directions, their opposites for an odd order; a minimum lies
      ! between the atoms: unequal sums of two directions as well
      npair = merge(min(max_pair_starts, nstart), 0, lowest)
      allocate (starts(3, 2*nstart + npair*npair))
      allocate (found(3, size(starts, 2)), value(size(starts, 2)))
      nfound = 0
      ns = 0
      do k = 1, nstart
         if (.not. r(order(k)) > 0.0_wp) cycle
         ns = ns + 1
         starts(:, ns) = d(:, order(k))/r(order(k))
         if (mod(m, 2) == 1) then
            ns = ns + 1
            starts(:, ns) = -starts(:, ns - 1)
         end if
      end do
      do k = 1, npair
         do j = 1, npair
            if (j == k .or. .not. min(r(order(k)), r(order(j))) > 0.0_wp) cycle
            v = d(:, order(k))/r(order(k)) + 0.5_wp*d(:, order(j))/r(order(j))
            if (.not. norm2(v) > 1.0e-6_wp) cycle
            ns = ns + 1
            starts(:, ns) = v/norm2(v)
         end do
      end do

      do istart = 1, ns
         v = starts(:, istart)
         do iter = 1, 200
            call moment_jet(d, m, v, f, g, h)
            v = sgn*g + alpha*v
            v = v/norm2(v)
         end do
         do iter = 1, 20
            call moment_jet(d, m, v, f, g, h)
            call setup_tangent_frame(v, b1, b2)
            a(1, 1) = dot_product(b1, matmul(h, b1)) - real(m, wp)*f
            a(2, 2) = dot_product(b2, matmul(h, b2)) - real(m, wp)*f
            a(1, 2) = dot_product(b1, matmul(h, b2))
            a(2, 1) = a(1, 2)
            rhs = -[dot_product(b1, g), dot_product(b2, g)]
            det = a(1, 1)*a(2, 2) - a(1, 2)*a(2, 1)
            if (abs(det) <= tiny(1.0_wp)) exit
            step = [a(2, 2)*rhs(1) - a(1, 2)*rhs(2), a(1, 1)*rhs(2) - a(2, 1)*rhs(1)]/det
            if (norm2(step) < 1.0e-15_wp) exit
            v = v + step(1)*b1 + step(2)*b2
            v = v/norm2(v)
         end do
         call moment_jet(d, m, v, f, g, h)
         ! An even moment does not tell a direction from its opposite
         if (mod(m, 2) == 0) v = lab_sense(v)
         known = .false.
         do k = 1, nfound
            known = dot_product(found(:, k), v) > 1.0_wp - 1.0e-8_wp
            if (known) exit
         end do
         if (known) then
            ! Keep the better converged copy of a direction
            if (sgn*f > value(k)) then
               found(:, k) = v
               value(k) = sgn*f
            end if
            cycle
         end if
         nfound = nfound + 1
         found(:, nfound) = v
         value(nfound) = sgn*f
      end do

      allocate (others(3, 0))
      if (nfound == 0) then
         fmax = sgn*fmax
         return
      end if
      ! Largest kept value: some direction always passes the tie test below
      fmax = maxval(value(1:nfound))
      ibest = 0
      best = -huge(1.0_wp)
      do k = 1, nfound
         if (value(k) < fmax - tie_zero*sum(powi(r, m))) cycle
         score = dot_product(found(:, k), number_dir)
         if (score > best) then
            best = score
            ibest = k
         end if
      end do
      u = found(:, ibest)
      do k = 1, nfound
         if (k == ibest .or. value(k) < fmax - tie_zero*sum(powi(r, m))) cycle
         if (size(others, 2) >= max_tied) exit
         others = reshape([others, found(:, k)], [3, size(others, 2) + 1])
      end do
      fmax = sgn*fmax
   end subroutine maximize_moment

   !> Value, gradient and Hessian of f(u) = sum_a (d_a . u)^m
   !>
   !> @param[in]  d Coordinates about the centroid (3, nat)
   !> @param[in]  m Moment order
   !> @param[in]  u Direction
   !> @param[out] f Value
   !> @param[out] g Gradient
   !> @param[out] h Hessian
   pure subroutine moment_jet(d, m, u, f, g, h)
      !> Coordinates about the centroid
      real(wp), intent(in) :: d(:, :)
      !> Moment order
      integer, intent(in) :: m
      !> Direction
      real(wp), intent(in) :: u(3)
      !> Value
      real(wp), intent(out) :: f
      !> Gradient
      real(wp), intent(out) :: g(3)
      !> Hessian
      real(wp), intent(out) :: h(3, 3)

      real(wp) :: x
      integer :: ia, i

      f = 0.0_wp
      g = 0.0_wp
      h = 0.0_wp
      do ia = 1, size(d, 2)
         x = dot_product(d(:, ia), u)
         f = f + powi(x, m)
         g = g + real(m, wp)*powi(x, m - 1)*d(:, ia)
         do i = 1, 3
            h(:, i) = h(:, i) + real(m*(m - 1), wp)*powi(x, m - 2)*d(i, ia)*d(:, ia)
         end do
      end do
   end subroutine moment_jet

   !> Reverse-mode derivative of the frame
   !>
   !> The frame rotates by delta omega under a displacement of the atoms,
   !> delta axes = [delta omega]_x axes; given the torque tau = dE/d(omega),
   !> returns dE/dxyz through the frame alone
   !>
   !> - `none` contributes nothing
   !> - without a harmonic (`axial`, or `sphere` without one) the in-plane
   !>   angle follows the lab reference direction, which twists about the
   !>   axis when the axis tilts
   !> - zero at an eigenvalue degeneracy of a level that uses eigenvectors is
   !>   never reached: such a degeneracy selects the next level
   !>
   !> @param[in]  self   Frame from [[new_sym_frame]] at xyz
   !> @param[in]  xyz    Coordinates (3, nat)
   !> @param[in]  torque Torque on the frame, dE/d(omega)
   !> @param[out] grad   Frame contribution to dE/dxyz (3, nat)
   subroutine sym_frame_vjp(self, xyz, torque, grad)
      !> Frame
      type(sym_frame_type), intent(in) :: self
      !> Coordinates
      real(wp), intent(in) :: xyz(:, :)
      !> Torque on the frame
      real(wp), intent(in) :: torque(3)
      !> Frame contribution to the gradient
      real(wp), intent(out) :: grad(:, :)

      !> Coordinates about the centroid
      real(wp), allocatable :: d(:, :)
      !> In-plane axes and axis as columns, and torque components
      real(wp) :: pqn(3, 3), n(3), p(3), q(3), tau_n, v(3)
      !> Eigenvectors in cyclic order and coefficient
      real(wp) :: ej(3), ek(3), coef
      !> Adjoint of the frame tensor
      real(wp) :: sbar(3, 3)
      !> In-plane weights
      complex(wp) :: z, w
      real(wp) :: h(3)
      integer :: nat, ia, i, j, k, m

      nat = size(xyz, 2)
      grad = 0.0_wp
      if (self%level == frame_level_none) return
      d = xyz - spread(sum(xyz, dim=2)/real(nat, wp), 2, nat)

      sbar = 0.0_wp
      if (self%level == frame_level_eigen) then
         ! omega_i = e_j^T dS e_k/(lambda_j - lambda_k), (i, j, k) cyclic
         do i = 1, 3
            j = modulo(i, 3) + 1
            k = modulo(i + 1, 3) + 1
            ej = self%axes(:, j)
            ek = self%axes(:, k)
            coef = 0.5_wp*dot_product(torque, self%axes(:, i))/(self%lambda(j) - self%lambda(k))
            do m = 1, 3
               sbar(:, m) = sbar(:, m) + coef*(ej*ek(m) + ek*ej(m))
            end do
         end do
         call tensor_vjp(d, sbar, grad, self%mix)
         do i = 1, 3
            grad(i, :) = grad(i, :) - sum(grad(i, :))/real(nat, wp)
         end do
         return
      end if

      select case (self%place)
      case (place_diagonal)
         pqn = matmul(self%axes, diag_basis)
      case (place_tetra)
         pqn = matmul(self%axes, tetra_basis)
      case (place_edge)
         pqn = matmul(self%axes, edge_basis)
      case default
         pqn = self%axes
      end select
      p = pqn(:, 1)
      q = pqn(:, 2)
      n = pqn(:, 3)
      tau_n = dot_product(torque, n)
      v = torque - tau_n*n

      ! In-plane angle: delta phi = sum_a g_a . delta d_a + h . omega_perp
      if (self%harmonic > 0 .and. self%level /= frame_level_axial) then
         m = self%harmonic
         h = 0.0_wp
         do ia = 1, nat
            z = cmplx(dot_product(d(:, ia), p), dot_product(d(:, ia), q), kind=wp)
            w = z**(m - 1)/cmplx(self%c_abs, 0.0_wp, kind=wp)
            grad(:, ia) = grad(:, ia) + tau_n*(aimag(w)*p + real(w, wp)*q)
            h = h + aimag(w)*cross_product(p, d(:, ia)) + real(w, wp)*cross_product(q, d(:, ia))
         end do
         v = v + tau_n*h
      end if

      ! Lab reference of the in-plane angle: p = lab_reference(n) stays in the
      ! plane of n and the lab axis e_k, so a tilt of n twists the frame about
      ! n by delta phi = (n_k/p_k) (omega_perp . p)
      if (self%harmonic == 0) then
         k = minloc(abs(n), dim=1)
         v = v + tau_n*(n(k)/p(k))*p
      end if

      ! Tilt of the axis: v . omega_perp with omega_perp = n x delta n
      select case (self%level)
      case (frame_level_plane, frame_level_axial)
         call axis_eigen_vjp(self, n, v, sbar)
         call tensor_vjp(d, sbar, grad, self%mix)
      case (frame_level_sphere)
         call axis_sphere_vjp(self, d, n, v, grad)
      case default
         ! `none` and `eigen` returned above
      end select

      ! d_a depends on every atom through the centroid
      do i = 1, 3
         grad(i, :) = grad(i, :) - sum(grad(i, :))/real(nat, wp)
      end do

   end subroutine sym_frame_vjp

   !> Reverse-mode tilt of an eigenvector axis
   !>
   !> delta n = sum_{j /= u} e_j (e_j^T dS n)/(lambda_u - lambda_j)
   !>
   !> @param[in]     self Frame
   !> @param[in]     n    Axis
   !> @param[in]     v    Adjoint of omega_perp
   !> @param[in,out] sbar Adjoint of the frame tensor, accumulated
   pure subroutine axis_eigen_vjp(self, n, v, sbar)
      !> Frame
      type(sym_frame_type), intent(in) :: self
      !> Axis
      real(wp), intent(in) :: n(3)
      !> Adjoint of omega_perp
      real(wp), intent(in) :: v(3)
      !> Adjoint of the frame tensor
      real(wp), intent(inout) :: sbar(3, 3)

      real(wp) :: ej(3), coef
      integer :: j, m

      do j = 1, 3
         if (j == self%iu) cycle
         ej = self%evec(:, j)
         coef = 0.5_wp*dot_product(v, cross_product(n, ej))/(self%lambda(self%iu) - self%lambda(j))
         do m = 1, 3
            sbar(:, m) = sbar(:, m) + coef*(ej*n(m) + n*ej(m))
         end do
      end do
   end subroutine axis_eigen_vjp

   !> Reverse-mode tilt of a `sphere` axis
   !>
   !> - moment maximum: (P H P - m f P) delta u = -P delta_d g
   !> - farthest atom k: delta u = P delta d_k / |d_k|
   !>
   !> @param[in]     self Frame
   !> @param[in]     d    Coordinates about the centroid (3, nat)
   !> @param[in]     u    Axis
   !> @param[in]     v    Adjoint of omega_perp
   !> @param[in,out] grad Gradient, accumulated
   subroutine axis_sphere_vjp(self, d, u, v, grad)
      !> Frame
      type(sym_frame_type), intent(in) :: self
      !> Coordinates about the centroid
      real(wp), intent(in) :: d(:, :)
      !> Axis
      real(wp), intent(in) :: u(3)
      !> Adjoint of omega_perp
      real(wp), intent(in) :: v(3)
      !> Gradient
      real(wp), intent(inout) :: grad(:, :)

      !> s = v x u, the adjoint of delta u, and t = A^+ s
      real(wp) :: s(3), t(3), b1(3), b2(3), a(2, 2), rhs(2), det, x(2)
      real(wp) :: f, g(3), h(3, 3), dk(3), r, proj
      integer :: m, ia, k

      s = cross_product(v, u)
      if (self%axis_order < 0) then
         k = self%axis_atom
         dk = d(:, k)
         r = norm2(dk)
         grad(:, k) = grad(:, k) + (s - u*dot_product(u, s))/r
         return
      end if

      m = self%axis_order
      call moment_jet(d, m, u, f, g, h)
      call setup_tangent_frame(u, b1, b2)
      a(1, 1) = dot_product(b1, matmul(h, b1)) - real(m, wp)*f
      a(2, 2) = dot_product(b2, matmul(h, b2)) - real(m, wp)*f
      a(1, 2) = dot_product(b1, matmul(h, b2))
      a(2, 1) = a(1, 2)
      rhs = [dot_product(b1, s), dot_product(b2, s)]
      det = a(1, 1)*a(2, 2) - a(1, 2)*a(2, 1)
      if (abs(det) <= tiny(1.0_wp)) return
      x = [a(2, 2)*rhs(1) - a(1, 2)*rhs(2), a(1, 1)*rhs(2) - a(2, 1)*rhs(1)]/det
      t = x(1)*b1 + x(2)*b2

      do ia = 1, size(d, 2)
         proj = dot_product(d(:, ia), u)
         grad(:, ia) = grad(:, ia) - real(m, wp)*(real(m - 1, wp)*powi(proj, m - 2) &
                                                 & *dot_product(t, d(:, ia))*u + powi(proj, m - 1)*t)
      end do
   end subroutine axis_sphere_vjp

   !> Lab axis least aligned with n, projected into the plane normal to n
   !>
   !> @param[in] n Unit vector
   pure function lab_reference(n) result(p)
      !> Unit vector
      real(wp), intent(in) :: n(3)
      !> Unit vector normal to n
      real(wp) :: p(3)

      integer :: k

      k = minloc(abs(n), dim=1)
      p = 0.0_wp
      p(k) = 1.0_wp
      p = p - n*dot_product(n, p)
      p = p/norm2(p)
   end function lab_reference

   !> Name of a molecular frame level
   !>
   !> @param[in] level One of the `frame_level_*` values
   pure function frame_level_name(level) result(name)
      !> Frame level
      integer, intent(in) :: level
      !> Name
      character(len=:), allocatable :: name

      select case (level)
      case (frame_level_none)
         name = "none"
      case (frame_level_eigen)
         name = "eigen"
      case (frame_level_plane)
         name = "plane"
      case (frame_level_sphere)
         name = "sphere"
      case (frame_level_axial)
         name = "axial"
      case default
         name = "unknown"
      end select
   end function frame_level_name

   !> Codes and names of the molecular frame levels, "0 none, 1 eigen, ..."
   pure function frame_level_list() result(text)
      !> Rendering
      character(len=:), allocatable :: text

      integer, parameter :: levels(5) = [frame_level_none, frame_level_eigen, frame_level_plane, &
         & frame_level_sphere, frame_level_axial]
      integer :: k

      text = ""
      do k = 1, size(levels)
         if (k > 1) text = text//", "
         text = text//to_string(levels(k))//" "//frame_level_name(levels(k))
      end do
   end function frame_level_list

end module moist_math_symmetry_frame
