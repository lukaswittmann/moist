!> Tangent of the DROP field contraction, and the jet tensors behind it
!>
!>     res(s, i) = w0 * T0(s, i) + sum_a w1(a) * T1(a, s, i)
!>                               + sum_ab w2(a, b) * T2(a, b, s, i)
!>
!> - `T0 = lsf1_rA = dS/dR_A`, `T1 = lsf2_r_rA = d^2S/(dr dR_A)`,
!>   `T2 = lsf3_rr_rA = d^3S/(dr^2 dR_A)`; slot `i` is `active_atom(i)`
!> - Tensors filled once per point by [[drop_field_jet_point]]:
!>   `39 n_active` doubles, one `f3_rr_rA` call
!> - Readers: [[drop_field_jet_contract]] (the row, as `vjp_f1_rA`),
!>   [[drop_field_jet_tangent]] (`sum_i v_i . T_k(.., i)`),
!>   [[drop_field_jet_column_packed]] (one unit direction, packed)
!> - Row tangent along `v`: weight tangents `(dw0, dw1, dw2)`, explicit
!>   nuclear motion (`hvp_jet_rA`), motion `dr` of the projected point
!> - Weighted explicit terms of all directions are the columns of the level
!>   set's `vjp_f2_rArB`; a caller holding them passes `explicit = .false.`
!> - Folding: weight tangents and two point-motion terms are one contraction,
!>   `row(dw0, dw1 + w0*dr, dw2 + outer(w1, dr))`
!> - `outer(w1, dr)(a, b) = w1(a) * dr(b)`, never symmetrised: the row
!>   contracts all nine `w2` entries, no symmetry, no factor of two
!> - Third term `w2_ab T3(a, b, k) dr_k` does not fold: `T3` is `f4_rrr_rA`,
!>   weights folded in once per point by [[drop_field_f4_fold]]
!> - Point half: [[drop_field_tangent_point]], then [[drop_field_f4_fold]];
!>   direction half: [[drop_field_tangent_dir]]
!> - [[drop_field_jet_point]] alone when no fourth derivative is read
!> - [[drop_field_tangent]] runs both halves for a single direction
!> - Fills are unconditional: buffers are reused across grid points; readers
!>   abort on a buffer not filled at this active count
!> - Tangent needs `hvp_jet_rA` and `f4_rrr_rA` (SvdW, CFC); isodensity LSFs
!>   abort inside the accessor, as intended
!> - Jet fill and its readers need `f3_rr_rA` only, provided by every DROP LSF
module moist_cavity_drop_derivatives_field_tangent
   use mctc_env_accuracy, only: wp
   use moist_cavity_drop_lsf_base, only: moist_cavity_drop_lsf_type, lsf_jet_row_entry

   implicit none(type, external)
   private

   public :: drop_field_tangent
   public :: drop_field_tangent_point, drop_field_f4_fold, drop_field_tangent_dir
   public :: drop_field_jet_point, drop_field_jet_contract
   public :: drop_field_jet_tangent, drop_field_jet_column_packed
   public :: drop_field_tangent_work_type
   public :: drop_n_sym2, drop_n_sym3, drop_n_jet_coef, drop_sym2_idx, drop_sym3_idx

   !> Spatial dimension
   integer, parameter :: ndim = 3

   !> Independent entries of symmetric `(3, 3)` and `(3, 3, 3)` spatial tensors
   integer, parameter :: drop_n_sym2 = 6, drop_n_sym3 = 10
   !> Packed jet coefficient count, `1 + 3 + 6 + 10`
   integer, parameter :: drop_n_jet_coef = 1 + ndim + drop_n_sym2 + drop_n_sym3

   !> Symmetry class representatives `a <= b`, lexicographic
   !>
   !> - With [[drop_sym3_idx]] the single definition of the packed order
   !> - Written by [[drop_field_jet_column_packed]], read by callers building
   !>   the matching symmetrised basis of directional jets
   integer, parameter :: drop_sym2_idx(2, drop_n_sym2) = reshape([ &
      & 1, 1, 1, 2, 1, 3, 2, 2, 2, 3, 3, 3], [2, drop_n_sym2])
   !> Symmetry class representatives `a <= b <= c`, lexicographic
   integer, parameter :: drop_sym3_idx(3, drop_n_sym3) = reshape([ &
      & 1, 1, 1, 1, 1, 2, 1, 1, 3, 1, 2, 2, 1, 2, 3, &
      & 1, 3, 3, 2, 2, 2, 2, 2, 3, 2, 3, 3, 3, 3, 3], [3, drop_n_sym3])

   !> Per-point jet tensors and scratch of the field contraction
   !>
   !> - Active-indexed, sized to the largest active count seen so far
   !> - Heap buffers, grown once and reused across grid points: no allocation
   !>   in steady state, no large automatic arrays inside the OpenMP grid loop
   !> - Slot markers are not cache keys: two points can share an active count,
   !>   so fills always refill
   !> - Markers only catch a reader whose fill never ran
   type :: drop_field_tangent_work_type
      !> Active slots the buffers are sized for
      integer :: capacity = 0
      !> Active slots of the last jet fill, `-1` before the first
      integer :: jet_slots = -1
      !> Active slots of the last `f4` fill, `-1` before the first
      integer :: f4_slots = -1
      !> Active slots of the last `f4w` fold, `-1` before the first
      integer :: f4w_slots = -1
      !> `dS/dR_A` [3, capacity]
      real(wp), allocatable :: t0(:, :)
      !> `d^2S/(dr dR_A)` [3, 3, capacity]
      real(wp), allocatable :: t1(:, :, :)
      !> `d^3S/(dr^2 dR_A)` [3, 3, 3, capacity]
      real(wp), allocatable :: t2(:, :, :, :)
      !> `sum_B v_B . d^2S/(dR_A dR_B)` [3, capacity]
      real(wp), allocatable :: hvp1(:, :)
      !> `sum_B v_B . d^3S/(dr dR_A dR_B)` [3, 3, capacity]
      real(wp), allocatable :: hvp2(:, :, :)
      !> `sum_B v_B . d^4S/(dr^2 dR_A dR_B)` [3, 3, 3, capacity]
      real(wp), allocatable :: hvp3(:, :, :, :)
      !> `d^4S/(dr^3 dR_A)` [3, 3, 3, 3, capacity]
      real(wp), allocatable :: f4(:, :, :, :, :)
      !> `sum_ab w2(a, b) d^4S/(dr_a dr_b dr_k dR_A)`, order `(k, s, i)` [3, 3, capacity]
      real(wp), allocatable :: f4w(:, :, :)
   contains
      !> Grow the buffers to hold at least `n` active slots
      procedure :: ensure => field_tangent_work_ensure
   end type drop_field_tangent_work_type

contains

   !* ================================================================================= *!
   !*                                  Scratch buffers                                  *!
   !* ================================================================================= *!

   !> Grow the scratch buffers to hold at least `n` active slots
   !>
   !> - Never shrinks: active counts fluctuate over a grid sweep, so
   !>   reallocation stays O(1) per sweep
   !>
   !> @param[in,out] self scratch buffers
   !> @param[in]     n    active slots the buffers must hold
   pure subroutine field_tangent_work_ensure(self, n)
      !> Scratch buffers
      class(drop_field_tangent_work_type), intent(inout) :: self
      !> Active slots the buffers must hold
      integer, intent(in) :: n

      if (allocated(self%f4) .and. self%capacity >= n) return

      if (allocated(self%t0)) deallocate (self%t0)
      if (allocated(self%t1)) deallocate (self%t1)
      if (allocated(self%t2)) deallocate (self%t2)
      if (allocated(self%hvp1)) deallocate (self%hvp1)
      if (allocated(self%hvp2)) deallocate (self%hvp2)
      if (allocated(self%hvp3)) deallocate (self%hvp3)
      if (allocated(self%f4)) deallocate (self%f4)
      if (allocated(self%f4w)) deallocate (self%f4w)

      allocate (self%t0(ndim, n))
      allocate (self%t1(ndim, ndim, n))
      allocate (self%t2(ndim, ndim, ndim, n))
      allocate (self%hvp1(ndim, n))
      allocate (self%hvp2(ndim, ndim, n))
      allocate (self%hvp3(ndim, ndim, ndim, n))
      allocate (self%f4(ndim, ndim, ndim, ndim, n))
      allocate (self%f4w(ndim, ndim, n))
      self%capacity = n
   end subroutine field_tangent_work_ensure

   !* ================================================================================= *!
   !*                              Point half: the fills                                *!
   !* ================================================================================= *!

   !> Materialise the jet tensors `T0`, `T1`, `T2` of the prepared point
   !>
   !> - One `f3_rr_rA` call at the point `lsf` is prepared at; reads no
   !>   direction, adjoint or point motion
   !> - Unconditional: rerun after every re-preparation of the LSF, before the
   !>   first contraction of that point; see the module header
   !> - Prepared order as `f3_rr_rA` requires: 2 for SvdW, 3 for CFC
   !>
   !> @param[in]     lsf  LSF instance, prepared at the evaluation point
   !> @param[in,out] work scratch buffers, reused across points
   subroutine drop_field_jet_point(lsf, work)
      !> LSF instance, prepared at the evaluation point
      class(moist_cavity_drop_lsf_type), intent(in) :: lsf
      !> Scratch buffers, reused across points
      type(drop_field_tangent_work_type), intent(inout) :: work

      !> Active slots of the prepared point
      integer :: n_active

      n_active = lsf%active_count()
      work%jet_slots = n_active
      if (n_active == 0) return

      call work%ensure(n_active)
      call lsf%f3_rr_rA(work%t0, work%t1, work%t2)
   end subroutine drop_field_jet_point

   !> Point half of the tangent: the jet tensors and `d^4S/(dr^3 dR_A)`
   !>
   !> - [[drop_field_jet_point]] plus the fill of `work%f4`: all
   !>   direction-free input of the direction half
   !> - Same unconditional-refill rule
   !> - Prepared order as `f4_rrr_rA` requires: 3 for SvdW, 4 for CFC; the
   !>   accessor checks and aborts otherwise
   !>
   !> @param[in]     lsf  LSF instance, prepared at the evaluation point
   !> @param[in,out] work scratch buffers, reused across points
   subroutine drop_field_tangent_point(lsf, work)
      !> LSF instance, prepared at the evaluation point
      class(moist_cavity_drop_lsf_type), intent(in) :: lsf
      !> Scratch buffers, reused across points
      type(drop_field_tangent_work_type), intent(inout) :: work

      !> Active slots of the prepared point
      integer :: n_active

      call drop_field_jet_point(lsf, work)

      n_active = lsf%active_count()
      work%f4_slots = n_active
      if (n_active == 0) return

      call lsf%f4_rrr_rA(work%f4)
   end subroutine drop_field_tangent_point

   !> Fold the point's Hessian weights into the mixed fourth derivative
   !>
   !> - `f4w(k, s, i) = sum_ab w2(a, b) f4(a, b, k, s, i)`, the point-motion
   !>   term the row cannot absorb; see the module header
   !> - Each direction then reads three numbers per slot instead of 27
   !> - Point-half fill: the weights belong to the point, not to a direction
   !> - Reads `f4`: run after [[drop_field_tangent_point]]
   !> - Unconditional refill; a new weight set at the same point needs a new fold
   !>
   !> @param[in,out] work     scratch buffers, `f4` filled at this point
   !> @param[in]     n_active active slots of the prepared point
   !> @param[in]     w2       adjoint weights of the spatial Hessian [3, 3]
   pure subroutine drop_field_f4_fold(work, n_active, w2)
      !> Scratch buffers, `f4` filled at this point
      type(drop_field_tangent_work_type), intent(inout) :: work
      !> Active slots of the prepared point
      integer, intent(in) :: n_active
      !> Adjoint weights of the spatial Hessian
      real(wp), intent(in) :: w2(3, 3)

      !> Fold accumulator
      real(wp) :: acc
      !> Active slot, nuclear axis and spatial axes
      integer :: i, s, k, a, b

      work%f4w_slots = n_active
      if (n_active == 0) return
      call assert_f4_filled(work, n_active, "drop_field_f4_fold")

      do i = 1, n_active
         do s = 1, ndim
            do k = 1, ndim
               acc = 0.0_wp
               do b = 1, ndim
                  do a = 1, ndim
                     acc = acc + w2(a, b)*work%f4(a, b, k, s, i)
                  end do
               end do
               work%f4w(k, s, i) = acc
            end do
         end do
      end do
   end subroutine drop_field_f4_fold

   !* ================================================================================= *!
   !*                          Readings of the jet tensors                              *!
   !* ================================================================================= *!

   !> Weighted nuclear-gradient row, `vjp_f1_rA` read off the buffer
   !>
   !> - Same contract as `vjp_f1_rA`, interchangeable at a call site
   !> - Columns `1 .. n_active` of `res` overwritten, the rest untouched;
   !>   nothing written for an empty active list
   !> - `w2` is a general `3 x 3`: all nine entries contracted, no symmetry
   !>   assumed, no factor of two folded in
   !>
   !> @param[in]     work     scratch buffers, jet tensors filled at this point
   !> @param[in]     n_active active slots of the prepared point
   !> @param[in]     w0       adjoint weight of the level-set value
   !> @param[in]     w1       adjoint weights of the spatial gradient [3]
   !> @param[in]     w2       adjoint weights of the spatial Hessian [3, 3]
   !> @param[in,out] res      weighted nuclear-gradient row [3, >= n_active]
   pure subroutine drop_field_jet_contract(work, n_active, w0, w1, w2, res)
      !> Scratch buffers, jet tensors filled at this point
      type(drop_field_tangent_work_type), intent(in) :: work
      !> Active slots of the prepared point
      integer, intent(in) :: n_active
      !> Adjoint weight of the level-set value
      real(wp), intent(in) :: w0
      !> Adjoint weights of the spatial gradient
      real(wp), intent(in) :: w1(3)
      !> Adjoint weights of the spatial Hessian
      real(wp), intent(in) :: w2(3, 3)
      !> Weighted nuclear-gradient row
      real(wp), intent(inout) :: res(:, :)

      !> Active slot and nuclear axis
      integer :: i, s

      if (n_active == 0) return
      call assert_jet_filled(work, n_active, "drop_field_jet_contract")

      do i = 1, n_active
         do s = 1, ndim
            res(s, i) = lsf_jet_row_entry(w0, w1, w2, work%t0(s, i), work%t1(:, s, i), &
                                          work%t2(:, :, s, i))
         end do
      end do
   end subroutine drop_field_jet_contract

   !> Directional nuclear derivative of the jet at the frozen point
   !>
   !> - `dv_k = sum_i sum_s v_act(s, i) T_k(.., s, i)`, as returned by
   !>   `tangent_f0`, `tangent_f1_r`, `tangent_f2_rr` and, through `f4`,
   !>   `tangent_f3_rrr`
   !> - Direction is slot indexed: `v_act(:, i)` displaces `active_atom(i)`
   !> - Caller owns the gather, the only use of the active-slot index space
   !>   in the traversal
   !> - Direction components off the active set contribute nothing, as in
   !>   the accessors
   !> - Optional `dv3` reads `f4`: only after [[drop_field_tangent_point]]; the
   !>   other three need [[drop_field_jet_point]]
   !> - All outputs fully written, zero for an empty active list
   !>
   !> @param[in]  work     scratch buffers, filled at this point
   !> @param[in]  n_active active slots of the prepared point
   !> @param[in]  v_act    direction restricted to the active slots [3, >= n_active]
   !> @param[out] dv0      directional derivative of the value
   !> @param[out] dv1      directional derivative of the spatial gradient [3]
   !> @param[out] dv2      directional derivative of the spatial Hessian [3, 3]
   !> @param[out] dv3      directional derivative of the third derivative [3, 3, 3]
   pure subroutine drop_field_jet_tangent(work, n_active, v_act, dv0, dv1, dv2, dv3)
      !> Scratch buffers, filled at this point
      type(drop_field_tangent_work_type), intent(in) :: work
      !> Active slots of the prepared point
      integer, intent(in) :: n_active
      !> Direction restricted to the active slots
      real(wp), intent(in) :: v_act(:, :)
      !> Directional derivative of the value
      real(wp), intent(out) :: dv0
      !> Directional derivative of the spatial gradient
      real(wp), intent(out) :: dv1(3)
      !> Directional derivative of the spatial Hessian
      real(wp), intent(out) :: dv2(3, 3)
      !> Directional derivative of the third spatial derivative
      real(wp), intent(out), optional :: dv3(3, 3, 3)

      !> Hoisted direction component
      real(wp) :: vs
      !> Active slot, nuclear axis and spatial axes
      integer :: i, s, a, b, c

      dv0 = 0.0_wp
      dv1 = 0.0_wp
      dv2 = 0.0_wp
      if (present(dv3)) dv3 = 0.0_wp
      if (n_active == 0) return

      call assert_jet_filled(work, n_active, "drop_field_jet_tangent")
      if (present(dv3)) call assert_f4_filled(work, n_active, "drop_field_jet_tangent")

      do i = 1, n_active
         do s = 1, ndim
            vs = v_act(s, i)
            dv0 = dv0 + vs*work%t0(s, i)
            do a = 1, ndim
               dv1(a) = dv1(a) + vs*work%t1(a, s, i)
            end do
            do b = 1, ndim
               do a = 1, ndim
                  dv2(a, b) = dv2(a, b) + vs*work%t2(a, b, s, i)
               end do
            end do
            if (present(dv3)) then
               do c = 1, ndim
                  do b = 1, ndim
                     do a = 1, ndim
                        dv3(a, b, c) = dv3(a, b, c) + vs*work%f4(a, b, c, s, i)
                     end do
                  end do
               end do
            end if
         end do
      end do
   end subroutine drop_field_jet_tangent

   !> Jet tangent of one unit direction, packed by spatial symmetry class
   !>
   !> - Jet tangent `(dv0, dv1, dv2, dv3)` of the unit direction `e_(s,i)`,
   !>   one column of the buffer
   !> - `(3, 3)` and `(3, 3, 3)` tangents reduced to one entry per symmetry
   !>   class, in the order of [[drop_sym2_idx]] and [[drop_sym3_idx]]
   !>
   !>     coef(1)      = dv0
   !>     coef(2:4)    = dv1(a)
   !>     coef(5:10)   = dv2(a, b),     a <= b
   !>     coef(11:20)  = dv3(a, b, c),  a <= b <= c
   !>
   !> - Coordinates in the symmetrised basis: a class element holds a one at
   !>   every index permutation of its representative
   !> - A linear map of the jet tangent, evaluated on that basis and
   !>   contracted with `coef`, gives the map of the full tangent
   !> - Valid since `t2` and `f4` are symmetric in the spatial indices: exact
   !>   for SvdW, to round-off for CFC
   !> - Reads `f4`: requires [[drop_field_tangent_point]]
   !>
   !> @param[in]  work     scratch buffers, filled at this point
   !> @param[in]  n_active active slots of the prepared point
   !> @param[in]  s        Cartesian axis of the direction
   !> @param[in]  i        active slot of the direction
   !> @param[out] coef     packed jet tangent [drop_n_jet_coef]
   pure subroutine drop_field_jet_column_packed(work, n_active, s, i, coef)
      !> Scratch buffers, filled at this point
      type(drop_field_tangent_work_type), intent(in) :: work
      !> Active slots of the prepared point
      integer, intent(in) :: n_active
      !> Cartesian axis and active slot of the direction
      integer, intent(in) :: s, i
      !> Packed jet tangent
      real(wp), intent(out) :: coef(drop_n_jet_coef)

      !> Symmetry class
      integer :: k

      call assert_jet_filled(work, n_active, "drop_field_jet_column_packed")
      call assert_f4_filled(work, n_active, "drop_field_jet_column_packed")

      coef(1) = work%t0(s, i)
      coef(2:1 + ndim) = work%t1(:, s, i)
      do k = 1, drop_n_sym2
         coef(1 + ndim + k) = work%t2(drop_sym2_idx(1, k), drop_sym2_idx(2, k), s, i)
      end do
      do k = 1, drop_n_sym3
         coef(1 + ndim + drop_n_sym2 + k) = work%f4(drop_sym3_idx(1, k), drop_sym3_idx(2, k), &
                                                    drop_sym3_idx(3, k), s, i)
      end do
   end subroutine drop_field_jet_column_packed

   !* ================================================================================= *!
   !*                            Field-contraction tangent                              *!
   !* ================================================================================= *!

   !> Directional derivative of the jet-contracted nuclear-gradient row
   !>
   !> - `d_v` of the weighted row along the nuclear direction `v`, from the
   !>   weight tangents `(dw0, dw1, dw2)` and the point motion `dr`
   !> - Terms and folding identity in the module header
   !> - `dr` is the caller's: tangent of the projected point `r*(p)` for the
   !>   cavity Hessian, `dr = 0` for the partial at a frozen point
   !> - `res` as in `vjp_f1_rA`: columns `1 .. active_count()` overwritten, the
   !>   rest untouched; nothing written for an empty active list
   !> - `res` is one point's row, not accumulated into; the caller scatter-adds
   !>   via `active_atom(i)`
   !> - Single-direction form: [[drop_field_tangent_point]],
   !>   [[drop_field_f4_fold]], then [[drop_field_tangent_dir]]
   !> - A sweep over a direction basis calls those itself, point half once
   !>
   !> @param[in]     lsf  LSF instance, prepared at the evaluation point
   !> @param[in]     w0   adjoint weight of the level-set value
   !> @param[in]     w1   adjoint weights of the spatial gradient [3]
   !> @param[in]     w2   adjoint weights of the spatial Hessian [3, 3]
   !> @param[in]     dw0  tangent of `w0` along `v`
   !> @param[in]     dw1  tangent of `w1` along `v` [3]
   !> @param[in]     dw2  tangent of `w2` along `v` [3, 3]
   !> @param[in]     dr   induced motion of the evaluation point along `v` [3]
   !> @param[in]     v    nuclear displacement directions [3, ncenters]
   !> @param[in,out] work scratch buffers, reused across points
   !> @param[in,out] res  tangent of the nuclear-gradient row [3, >= n_active]
   subroutine drop_field_tangent(lsf, w0, w1, w2, dw0, dw1, dw2, dr, v, work, res)
      !> LSF instance, prepared at the evaluation point
      class(moist_cavity_drop_lsf_type), intent(in) :: lsf
      !> Adjoint weight of the level-set value
      real(wp), intent(in) :: w0
      !> Adjoint weights of the spatial gradient
      real(wp), intent(in) :: w1(3)
      !> Adjoint weights of the spatial Hessian
      real(wp), intent(in) :: w2(3, 3)
      !> Tangent of `w0` along `v`
      real(wp), intent(in) :: dw0
      !> Tangent of `w1` along `v`
      real(wp), intent(in) :: dw1(3)
      !> Tangent of `w2` along `v`
      real(wp), intent(in) :: dw2(3, 3)
      !> Induced motion of the evaluation point along `v`
      real(wp), intent(in) :: dr(3)
      !> Nuclear displacement directions
      real(wp), intent(in) :: v(:, :)
      !> Scratch buffers, reused across points
      type(drop_field_tangent_work_type), intent(inout) :: work
      !> Tangent of the nuclear-gradient row
      real(wp), intent(inout) :: res(:, :)

      call drop_field_tangent_point(lsf, work)
      call drop_field_f4_fold(work, lsf%active_count(), w2)
      call drop_field_tangent_dir(lsf, w0, w1, w2, dw0, dw1, dw2, dr, v, .true., work, res)
   end subroutine drop_field_tangent

   !> Direction half: everything in the tangent that reads `v`, `dr` or a weight
   !>
   !> - Body of [[drop_field_tangent]] less the fills, same contract for `res`
   !> - [[drop_field_tangent_point]] and [[drop_field_f4_fold]] must have run
   !>   at this evaluation point; checked through the slot markers
   !> - Only LSF call is `hvp_jet_rA`, the explicit nuclear motion
   !> - `explicit = .false.`: explicit term left out and `v` not read; the
   !>   caller adds its column of the level set's `vjp_f2_rArB`
   !>
   !> @param[in]     lsf      LSF instance, prepared at the evaluation point
   !> @param[in]     w0       adjoint weight of the level-set value
   !> @param[in]     w1       adjoint weights of the spatial gradient [3]
   !> @param[in]     w2       adjoint weights of the spatial Hessian [3, 3]
   !> @param[in]     dw0      tangent of `w0` along `v`
   !> @param[in]     dw1      tangent of `w1` along `v` [3]
   !> @param[in]     dw2      tangent of `w2` along `v` [3, 3]
   !> @param[in]     dr       induced motion of the evaluation point along `v` [3]
   !> @param[in]     v        nuclear displacement directions [3, ncenters]
   !> @param[in]     explicit include the explicit term, `hvp_jet_rA` along `v`
   !> @param[in,out] work     scratch buffers, filled and folded at this point
   !> @param[in,out] res      tangent of the nuclear-gradient row [3, >= n_active]
   subroutine drop_field_tangent_dir(lsf, w0, w1, w2, dw0, dw1, dw2, dr, v, explicit, work, res)
      !> LSF instance, prepared at the evaluation point
      class(moist_cavity_drop_lsf_type), intent(in) :: lsf
      !> Adjoint weight of the level-set value
      real(wp), intent(in) :: w0
      !> Adjoint weights of the spatial gradient
      real(wp), intent(in) :: w1(3)
      !> Adjoint weights of the spatial Hessian
      real(wp), intent(in) :: w2(3, 3)
      !> Tangent of `w0` along `v`
      real(wp), intent(in) :: dw0
      !> Tangent of `w1` along `v`
      real(wp), intent(in) :: dw1(3)
      !> Tangent of `w2` along `v`
      real(wp), intent(in) :: dw2(3, 3)
      !> Induced motion of the evaluation point along `v`
      real(wp), intent(in) :: dr(3)
      !> Nuclear displacement directions
      real(wp), intent(in) :: v(:, :)
      !> Include the explicit nuclear-motion term
      logical, intent(in) :: explicit
      !> Scratch buffers, filled and folded at this point
      type(drop_field_tangent_work_type), intent(inout) :: work
      !> Tangent of the nuclear-gradient row
      real(wp), intent(inout) :: res(:, :)

      !> Shifted weights of the folded contraction
      real(wp) :: w1_fold(ndim), w2_fold(ndim, ndim)
      !> Row accumulator
      real(wp) :: acc
      !> Active slot, nuclear axis and spatial axes
      integer :: i, s, a, b, k
      !> Active slots of the prepared point
      integer :: n_active

      n_active = lsf%active_count()
      if (n_active == 0) return
      call assert_jet_filled(work, n_active, "drop_field_tangent_dir")
      call assert_f4w_filled(work, n_active, "drop_field_tangent_dir")

      !* ---------------------- Weight tangents and folded motion --------------------- *!

      ! `w0 * T1_k dr_k` joins the gradient slot, `w1_a * T2_ak dr_k` the
      ! Hessian slot: one contraction with the weight tangents
      ! `outer(w1, dr)` is not symmetrised; see the module header
      w1_fold = dw1 + w0*dr
      do b = 1, ndim
         do a = 1, ndim
            w2_fold(a, b) = dw2(a, b) + w1(a)*dr(b)
         end do
      end do
      call drop_field_jet_contract(work, n_active, dw0, w1_fold, w2_fold, res)

      !* --------------------------- Explicit nuclear motion -------------------------- *!

      if (explicit) then
         call lsf%hvp_jet_rA(v, work%hvp1, work%hvp2, work%hvp3)
         do i = 1, n_active
            do s = 1, ndim
               res(s, i) = res(s, i) + lsf_jet_row_entry(w0, w1, w2, work%hvp1(s, i), &
                                                         work%hvp2(:, s, i), work%hvp3(:, :, s, i))
            end do
         end do
      end if

      !* ------------------------ Unfoldable point-motion term ------------------------ *!

      ! `w2_ab T3_abk dr_k` has three spatial indices, one more than the row
      ! absorbs; the point half folded the weights into `f4w`
      do i = 1, n_active
         do s = 1, ndim
            acc = 0.0_wp
            do k = 1, ndim
               acc = acc + dr(k)*work%f4w(k, s, i)
            end do
            res(s, i) = res(s, i) + acc
         end do
      end do
   end subroutine drop_field_tangent_dir

   !* ================================================================================= *!
   !*                                  Fill assertions                                  *!
   !* ================================================================================= *!

   !> Abort on a jet buffer whose fill did not run at this active count
   !>
   !> @param[in] work     scratch buffers
   !> @param[in] n_active active slots the reader is about to assume
   !> @param[in] caller   reader named in the diagnostic
   pure subroutine assert_jet_filled(work, n_active, caller)
      !> Scratch buffers
      type(drop_field_tangent_work_type), intent(in) :: work
      !> Active slots the reader is about to assume
      integer, intent(in) :: n_active
      !> Reader named in the diagnostic
      character(len=*), intent(in) :: caller

      if (work%jet_slots /= n_active) then
         error stop "moist DROP field tangent: "//caller//" ran without a "// &
            "drop_field_jet_point at this evaluation point"
      end if
   end subroutine assert_jet_filled

   !> Abort on an `f4` buffer whose fill did not run at this active count
   !>
   !> @param[in] work     scratch buffers
   !> @param[in] n_active active slots the reader is about to assume
   !> @param[in] caller   reader named in the diagnostic
   pure subroutine assert_f4_filled(work, n_active, caller)
      !> Scratch buffers
      type(drop_field_tangent_work_type), intent(in) :: work
      !> Active slots the reader is about to assume
      integer, intent(in) :: n_active
      !> Reader named in the diagnostic
      character(len=*), intent(in) :: caller

      if (work%f4_slots /= n_active) then
         error stop "moist DROP field tangent: "//caller//" ran without a "// &
            "drop_field_tangent_point at this evaluation point"
      end if
   end subroutine assert_f4_filled

   !> Abort on an `f4w` buffer whose fold did not run at this active count
   !>
   !> @param[in] work     scratch buffers
   !> @param[in] n_active active slots the reader is about to assume
   !> @param[in] caller   reader named in the diagnostic
   pure subroutine assert_f4w_filled(work, n_active, caller)
      !> Scratch buffers
      type(drop_field_tangent_work_type), intent(in) :: work
      !> Active slots the reader is about to assume
      integer, intent(in) :: n_active
      !> Reader named in the diagnostic
      character(len=*), intent(in) :: caller

      if (work%f4w_slots /= n_active) then
         error stop "moist DROP field tangent: "//caller//" ran without a "// &
            "drop_field_f4_fold at this evaluation point"
      end if
   end subroutine assert_f4w_filled

end module moist_cavity_drop_derivatives_field_tangent
