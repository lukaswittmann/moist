!> Tangent of the DROP field contraction, and the jet tensors behind it
!>
!> The reverse-mode nuclear gradient contracts each grid point's level-set
!> adjoints `(w0, w1, w2)` into a nuclear-gradient row (the field channel of
!> `derivatives/nuclear.f90`):
!>
!>     res(s, i) = w0 * T0(s, i) + sum_a w1(a) * T1(a, s, i)
!>                               + sum_ab w2(a, b) * T2(a, b, s, i)
!>
!>     T0(s, i)       = lsf1_rA(s, i)      = dS / dR_(s,A)
!>     T1(a, s, i)    = lsf2_r_rA(a, s, i) = d^2 S / (dr_a dR_(s,A))
!>     T2(a, b, s, i) = lsf3_rr_rA(a,b,s,i) = d^3 S / (dr_a dr_b dR_(s,A))
!>
!> with the LSF's compact active index: slot `i` belongs to `active_atom(i)`.
!> The gradient path has the level set contract the jet indices itself, once
!> per point, through `vjp_f1_rA`, and never materialises the three tensors.
!>
!> The second-order traversals are different: at one grid point they contract
!> the *same* three tensors against many weight sets -- one per basis direction
!> of the fixed channel, one per nuclear direction of the response channel --
!> and they need the tensors' directional derivatives as well. For them a
!> contracted accessor is the wrong primitive, because every call re-runs the
!> `O(n_active)` kernel that produces the tensors. This module therefore
!> materialises `T0`, `T1` and `T2` once per point (`39 n_active` doubles, one
!> `f3_rr_rA` call, [[drop_field_jet_point]]) and offers the three readings the
!> traversals need as plain contractions of that buffer:
!>
!>   * [[drop_field_jet_contract]] -- the weighted row above, `vjp_f1_rA` read
!>     off the buffer;
!>   * [[drop_field_jet_tangent]] -- the directional nuclear derivative of the
!>     jet at the frozen point, `sum_i v_i . T_k(.., i)`, which is what the
!>     `tangent_f0 / f1_r / f2_rr / f3_rrr` accessors compute;
!>   * [[drop_field_jet_column]] -- the same for a Cartesian unit direction,
!>     which is a single column of the buffer.
!>
!> Directional derivative of the row
!> ---------------------------------
!> Two things move when the nuclei move along `v`: the adjoint weights, whose
!> tangents `(dw0, dw1, dw2)` the caller hands in, and the level-set jet
!> itself. The jet moves for two reasons at once, because the row is evaluated
!> at the *projected* point `r*(p)` and that point rides along:
!>
!>     d_v[T0]    = d^2S/(dR_A dR_B) v_B          + sum_k T1(k, ..) dr_k
!>     d_v[T1_a]  = d^3S/(dr_a dR_A dR_B) v_B     + sum_k T2(a, k, ..) dr_k
!>     d_v[T2_ab] = d^4S/(dr_a dr_b dR_A dR_B) v_B + sum_k T3(a, b, k, ..) dr_k
!>
!> The explicit halves are the nuclear Hessian-vector family the LSF returns
!> in one pass as `hvp_jet_rA`; the fourth-order spatial tensor `T3` is
!> `f4_rrr_rA`. Weighted, the explicit halves of *every* direction at a point
!> are the columns of one matrix, the level set's `vjp_f2_rArB`; a caller
!> holding that block tells [[drop_field_tangent_dir]] to leave the explicit
!> term out and adds its column itself.
!>
!> The folding identity
!> --------------------
!> Raising the spatial order of a term by one is precisely what the `w1` and
!> `w2` slots of the row already do, so two of the three point-motion terms are
!> not separate contractions at all:
!>
!>     w0 * sum_k T1(k, s, i) dr_k       = row(0, w0*dr, 0)
!>     sum_a w1(a) sum_k T2(a, k, s, i) dr_k = row(0, 0, outer(w1, dr))
!>
!> Both share the weight slots of the explicit weight-tangent term, so the
!> whole block collapses into a *single* contraction with shifted weights,
!>
!>     row(dw0, dw1 + w0*dr, dw2 + outer(w1, dr))
!>
!> where `outer(w1, dr)(a, b) = w1(a) * dr(b)`. That matrix is not symmetric,
!> and must not be symmetrised: the row contracts all nine entries of `w2` with
!> no symmetry assumption and no folded factor of two, so any redistribution
!> across the diagonal would silently change the answer for a caller whose own
!> `dw2` is asymmetric. (`T2` happens to be symmetric in `(a, b)`, so the
!> ordering of the outer product is numerically inert today; it is written in
!> the order the contraction derives it in anyway.)
!>
!> The third point-motion term does *not* fold. The row stops at two spatial
!> indices, and `w2_ab T3(a,b,k) dr_k` needs three, so it goes through
!> `f4_rrr_rA` directly. That is the one reason this module holds a
!> `(3, 3, 3, 3, n_active)` tensor. Its weights `w2` belong to the point, not
!> to the direction, so [[drop_field_f4_fold]] contracts them in once per
!> point and every direction reads `dr_k f4w(k, s, i)`: three numbers per slot
!> instead of twenty-seven.
!>
!> Splitting the point from the direction
!> --------------------------------------
!> `f3_rr_rA` and `f4_rrr_rA` take no `v`: both are functions of the prepared
!> point alone, and `f4_rrr_rA` is by a wide margin the most expensive accessor
!> called here -- for SvdW it re-runs the fourth-order atom tensors and the
!> third-order nuclear evaluation of every active atom. A caller sweeping a
!> whole basis of nuclear directions at one point would rebuild bit-identical
!> tensors once per direction, so the fills are [[drop_field_tangent_point]]
!> (both tensors), [[drop_field_jet_point]] (the jet alone, for a caller that
!> never reads a fourth derivative) and [[drop_field_f4_fold]] (the weighted
!> fourth derivative, once the point's weights are known), and the
!> direction-dependent remainder is [[drop_field_tangent_dir]].
!> [[drop_field_tangent]] is the fills and the direction in a row, for a
!> caller with one direction.
!>
!> **The fills are unconditional.** The buffers are reused across grid points,
!> so a fill skipped at one point would serve the previous point's tensors at
!> the next; every reader checks the slot markers and aborts on a buffer whose
!> fill did not run at this active count, which is the cheap part of that
!> invariant, not the whole of it.
!>
!> Applicability
!> -------------
!> This is a primitive, not a driver: it neither reads cavity state nor
!> reconstructs the adjoints, and it does not know where `dr` came from. The
!> tangent needs the full optional derivative set (`hvp_jet_rA` and
!> `f4_rrr_rA`), which SvdW and CFC implement and the isodensity LSFs do not --
!> passing one of those aborts inside the accessor rather than returning a zero
!> row, which is the intended behaviour. The jet fill and its contractions need
!> `f3_rr_rA` only, which every DROP level set provides.
module moist_cavity_drop_derivatives_field_tangent
   use mctc_env_accuracy, only: wp
   use moist_cavity_drop_lsf_base, only: moist_cavity_drop_lsf_type

   implicit none(type, external)
   private

   public :: drop_field_tangent
   public :: drop_field_tangent_point, drop_field_f4_fold, drop_field_tangent_dir
   public :: drop_field_jet_point, drop_field_jet_contract
   public :: drop_field_jet_tangent, drop_field_jet_column
   public :: drop_field_tangent_work_type

   !> Spatial dimension
   integer, parameter :: ndim = 3

   !> Per-point jet tensors and scratch of the field contraction
   !>
   !> Everything here is active-indexed and sized to the largest active count
   !> seen so far; holding the tensors as automatic arrays would put well over a
   !> hundred doubles per active atom on the stack of every call, inside what is
   !> ultimately an OpenMP grid loop. Held here the buffers grow once and are
   !> reused for every point of the grid, so the primitives are allocation-free
   !> in steady state.
   !>
   !> The two slot markers are deliberately not cache keys. The fills refill
   !> unconditionally, because two different evaluation points can share an
   !> active count and a fill keyed on that count would silently serve the wrong
   !> tensor. They are read only to catch a reader whose fill never ran.
   type :: drop_field_tangent_work_type
      !> Active slots the buffers are currently sized for
      integer :: capacity = 0
      !> Active slots the jet tensors hold a fill for; `-1` before the first fill
      integer :: jet_slots = -1
      !> Active slots `f4` holds a fill for; `-1` before the first fill
      integer :: f4_slots = -1
      !> Active slots `f4w` holds a fold for; `-1` before the first fold
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
      !> `sum_ab w2(a, b) d^4S/(dr_a dr_b dr_k dR_A)`, the point's Hessian weights
      !> folded into `f4`; index order `(k, s, i)` [3, 3, capacity]
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
   !> Never shrinks: the caller sweeps a grid whose per-point active counts
   !> fluctuate, and a buffer that shrank would reallocate on the next larger
   !> point. Reallocation therefore happens O(1) times per grid sweep.
   !>
   !> @param[inout] self Scratch buffers
   !> @param[in]    n    Active slots the buffers must hold
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
   !> One `f3_rr_rA` call at whatever point `lsf` is currently prepared at.
   !> Nothing here reads a direction, an adjoint or a point motion. It must be
   !> called again after *every* re-preparation of the LSF and before the first
   !> contraction of that point, and it is not conditional on anything -- see
   !> the module header. The prepared order must be what the level set requires
   !> for `f3_rr_rA` (2 for SvdW, 3 for CFC).
   !>
   !> @param[in]    lsf  LSF instance, prepared at the evaluation point
   !> @param[inout] work Scratch buffers, reused across points
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
   !> [[drop_field_jet_point]] plus the fill of `work%f4`, i.e. everything the
   !> direction half reads that does not depend on the direction. Same
   !> unconditional-refill rule. The LSF must be prepared at an order high
   !> enough for `f4_rrr_rA` (3 for SvdW, 4 for CFC); the accessor checks that
   !> itself and aborts if not.
   !>
   !> @param[in]    lsf  LSF instance, prepared at the evaluation point
   !> @param[inout] work Scratch buffers, reused across points
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
   !> `f4w(k, s, i) = sum_ab w2(a, b) f4(a, b, k, s, i)`: the one point-motion
   !> term the row cannot absorb (see the module header), contracted with the
   !> weights once so that every direction at the point reads three numbers per
   !> slot instead of twenty-seven. The weights are a property of the point and
   !> of the adjoint held fixed, not of a direction, which is what makes this a
   !> point-half fill. It reads `f4`, so it follows [[drop_field_tangent_point]],
   !> and it obeys the same unconditional-refill rule with one addition: a new
   !> weight set at the same point needs a new fold.
   !>
   !> @param[inout] work     Scratch buffers, `f4` filled at this point
   !> @param[in]    n_active Active slots of the prepared point
   !> @param[in]    w2       Adjoint weights of the spatial Hessian [3, 3]
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

   !> The weighted nuclear-gradient row, `vjp_f1_rA` read off the buffer
   !>
   !> Follows `vjp_f1_rA`'s contract exactly, so the two are interchangeable at
   !> a call site: columns `1 .. n_active` of `res` are *overwritten*, columns
   !> beyond that are left untouched, and nothing is written when the active
   !> list is empty. `w2` is a general `3 x 3`: all nine entries are contracted,
   !> no symmetry is assumed and no factor of two is folded in.
   !>
   !> @param[in]    work     Scratch buffers, jet tensors filled at this point
   !> @param[in]    n_active Active slots of the prepared point
   !> @param[in]    w0       Adjoint weight of the level-set value
   !> @param[in]    w1       Adjoint weights of the spatial gradient [3]
   !> @param[in]    w2       Adjoint weights of the spatial Hessian [3, 3]
   !> @param[inout] res      Weighted nuclear-gradient row [3, >= n_active]
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

      !> Row accumulator
      real(wp) :: acc
      !> Active slot, nuclear axis and spatial axes
      integer :: i, s, a, b

      if (n_active == 0) return
      call assert_jet_filled(work, n_active, "drop_field_jet_contract")

      do i = 1, n_active
         do s = 1, ndim
            acc = w0*work%t0(s, i)
            do a = 1, ndim
               acc = acc + w1(a)*work%t1(a, s, i)
            end do
            do b = 1, ndim
               do a = 1, ndim
                  acc = acc + w2(a, b)*work%t2(a, b, s, i)
               end do
            end do
            res(s, i) = acc
         end do
      end do
   end subroutine drop_field_jet_contract

   !> Directional nuclear derivative of the jet at the frozen point
   !>
   !> `dv_k = sum_i sum_s v_act(s, i) T_k(.., s, i)`: what `tangent_f0`,
   !> `tangent_f1_r`, `tangent_f2_rr` and -- through `f4` -- `tangent_f3_rrr`
   !> return for the direction whose active-slot restriction is `v_act`. The
   !> direction is handed in *slot indexed*, `v_act(:, i)` being the
   !> displacement of `active_atom(i)`; the caller owns that gather, and the
   !> active-slot index space appears in no other place of the traversal.
   !> Components of the direction on atoms outside the active set contribute
   !> nothing here, exactly as they contribute nothing to the accessors.
   !>
   !> `dv3` is optional and reads `f4`, so it may be asked for only after
   !> [[drop_field_tangent_point]]; the other three need [[drop_field_jet_point]].
   !> All outputs are fully written, zero when the active list is empty.
   !>
   !> @param[in]  work     Scratch buffers, filled at this point
   !> @param[in]  n_active Active slots of the prepared point
   !> @param[in]  v_act    Direction restricted to the active slots [3, >= n_active]
   !> @param[out] dv0      Directional derivative of the value
   !> @param[out] dv1      Directional derivative of the spatial gradient [3]
   !> @param[out] dv2      Directional derivative of the spatial Hessian [3, 3]
   !> @param[out] dv3      Directional derivative of the third derivative [3, 3, 3]
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

   !> The jet tangent along the Cartesian unit direction of one active slot
   !>
   !> [[drop_field_jet_tangent]] for `v_act = e_(s, i)`, which is a single
   !> column of every buffer. Same fill requirements: `dv3` reads `f4`.
   !>
   !> @param[in]  work     Scratch buffers, filled at this point
   !> @param[in]  n_active Active slots of the prepared point
   !> @param[in]  s        Cartesian axis of the direction
   !> @param[in]  i        Active slot of the direction
   !> @param[out] dv0      Directional derivative of the value
   !> @param[out] dv1      Directional derivative of the spatial gradient [3]
   !> @param[out] dv2      Directional derivative of the spatial Hessian [3, 3]
   !> @param[out] dv3      Directional derivative of the third derivative [3, 3, 3]
   pure subroutine drop_field_jet_column(work, n_active, s, i, dv0, dv1, dv2, dv3)
      !> Scratch buffers, filled at this point
      type(drop_field_tangent_work_type), intent(in) :: work
      !> Active slots of the prepared point
      integer, intent(in) :: n_active
      !> Cartesian axis and active slot of the direction
      integer, intent(in) :: s, i
      !> Directional derivative of the value
      real(wp), intent(out) :: dv0
      !> Directional derivative of the spatial gradient
      real(wp), intent(out) :: dv1(3)
      !> Directional derivative of the spatial Hessian
      real(wp), intent(out) :: dv2(3, 3)
      !> Directional derivative of the third spatial derivative
      real(wp), intent(out), optional :: dv3(3, 3, 3)

      call assert_jet_filled(work, n_active, "drop_field_jet_column")
      if (present(dv3)) call assert_f4_filled(work, n_active, "drop_field_jet_column")

      dv0 = work%t0(s, i)
      dv1 = work%t1(:, s, i)
      dv2 = work%t2(:, :, s, i)
      if (present(dv3)) dv3 = work%f4(:, :, :, s, i)
   end subroutine drop_field_jet_column

   !* ================================================================================= *!
   !*                            Field-contraction tangent                              *!
   !* ================================================================================= *!

   !> Directional derivative of the jet-contracted nuclear-gradient row
   !>
   !> Returns `d_v` of the weighted row for the nuclear direction `v`, given the
   !> tangents `(dw0, dw1, dw2)` of the adjoint weights and the induced motion
   !> `dr` of the evaluation point. See the module header for the derivation
   !> and for the folding identity that turns the weight-tangent block and two
   !> of the three point-motion terms into a single contraction.
   !>
   !> `dr` is the caller's business. It is the tangent of the projected point
   !> `r*(p)` along `v` for the field row of the cavity Hessian, but nothing
   !> here assumes that; a caller wanting the partial derivative at a frozen
   !> evaluation point passes `dr = 0`.
   !>
   !> `res` follows `vjp_f1_rA`'s contract exactly: columns `1 .. active_count()`
   !> are *overwritten* with the tangent, columns beyond that are left
   !> untouched, and nothing at all is written when the active list is empty.
   !> The result is not accumulated into -- like `vjp_pt` in the gradient
   !> driver, the buffer is one point's row and the caller owns the scatter-add
   !> into the per-atom accumulator via `active_atom(i)`.
   !>
   !> This is the single-direction form: [[drop_field_tangent_point]] and
   !> [[drop_field_f4_fold]] followed by [[drop_field_tangent_dir]]. A caller
   !> sweeping a direction basis at one point calls those itself and pays for
   !> the point half once.
   !>
   !> @param[in]    lsf  LSF instance, prepared at the evaluation point
   !> @param[in]    w0   Adjoint weight of the level-set value
   !> @param[in]    w1   Adjoint weights of the spatial gradient [3]
   !> @param[in]    w2   Adjoint weights of the spatial Hessian [3, 3]
   !> @param[in]    dw0  Tangent of `w0` along `v`
   !> @param[in]    dw1  Tangent of `w1` along `v` [3]
   !> @param[in]    dw2  Tangent of `w2` along `v` [3, 3]
   !> @param[in]    dr   Induced motion of the evaluation point along `v` [3]
   !> @param[in]    v    Nuclear displacement directions [3, ncenters]
   !> @param[inout] work Scratch buffers, reused across points
   !> @param[inout] res  Tangent of the nuclear-gradient row [3, >= n_active]
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
   !> The body of [[drop_field_tangent]], less the fills, and with exactly that
   !> routine's contract: see it for the derivation, the folding identity and
   !> what `res` is. [[drop_field_tangent_point]] must have run at this
   !> evaluation point first, which is checked rather than assumed -- a stale
   !> tensor would give a wrong tangent with nothing to see it.
   !>
   !> The one LSF call left in here is `hvp_jet_rA`, the explicit nuclear
   !> motion: it is the only piece of the tangent that is not a reading of a
   !> direction-free tensor of the point. A caller that holds the weighted
   !> explicit term for its whole direction set at once -- the columns of the
   !> level set's `vjp_f2_rArB` -- passes `explicit = .false.`, gets the rest of
   !> the tangent, and adds that column itself; `v` is then not read.
   !>
   !> @param[in]    lsf      LSF instance, prepared at the evaluation point
   !> @param[in]    w0       Adjoint weight of the level-set value
   !> @param[in]    w1       Adjoint weights of the spatial gradient [3]
   !> @param[in]    w2       Adjoint weights of the spatial Hessian [3, 3]
   !> @param[in]    dw0      Tangent of `w0` along `v`
   !> @param[in]    dw1      Tangent of `w1` along `v` [3]
   !> @param[in]    dw2      Tangent of `w2` along `v` [3, 3]
   !> @param[in]    dr       Induced motion of the evaluation point along `v` [3]
   !> @param[in]    v        Nuclear displacement directions [3, ncenters]
   !> @param[in]    explicit Include the explicit nuclear-motion term, `hvp_jet_rA` along `v`
   !> @param[inout] work     Scratch buffers, filled and folded at this point
   !> @param[inout] res      Tangent of the nuclear-gradient row [3, >= n_active]
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
      !> Whether the explicit nuclear-motion term is included
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

      ! `w0 * T1_k dr_k` rides in the gradient slot and `w1_a * T2_ak dr_k` in
      ! the Hessian slot, so both join the explicit weight tangents in one
      ! contraction. `outer(w1, dr)` is deliberately not symmetrised; see the
      ! module header.
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
               acc = w0*work%hvp1(s, i)
               do a = 1, ndim
                  acc = acc + w1(a)*work%hvp2(a, s, i)
               end do
               do b = 1, ndim
                  do a = 1, ndim
                     acc = acc + w2(a, b)*work%hvp3(a, b, s, i)
                  end do
               end do
               res(s, i) = res(s, i) + acc
            end do
         end do
      end if

      !* ------------------------ Unfoldable point-motion term ------------------------ *!

      ! `w2_ab T3_abk dr_k` carries three spatial indices, one more than the
      ! row can absorb; the weights were folded into the mixed fourth
      ! derivative by the point half, once, and the direction reads the fold.
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
   !> @param[in] work     Scratch buffers
   !> @param[in] n_active Active slots the reader is about to assume
   !> @param[in] caller   Reader named in the diagnostic
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
   !> @param[in] work     Scratch buffers
   !> @param[in] n_active Active slots the reader is about to assume
   !> @param[in] caller   Reader named in the diagnostic
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
   !> @param[in] work     Scratch buffers
   !> @param[in] n_active Active slots the reader is about to assume
   !> @param[in] caller   Reader named in the diagnostic
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
