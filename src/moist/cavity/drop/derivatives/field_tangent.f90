!> Tangent of the DROP field contraction
!>
!> The reverse-mode nuclear gradient contracts each grid point's level-set
!> adjoints `(w0, w1, w2)` into a nuclear-gradient row with one `vjp_f1_rA`
!> call (see the field channel of `derivatives/nuclear.f90`):
!>
!>     res(s, i) = w0 * T0(s, i) + sum_a w1(a) * T1(a, s, i)
!>                               + sum_ab w2(a, b) * T2(a, b, s, i)
!>
!>     T0(s, i)       = lsf1_rA(s, i)      = dS / dR_(s,A)
!>     T1(a, s, i)    = lsf2_r_rA(a, s, i) = d^2 S / (dr_a dR_(s,A))
!>     T2(a, b, s, i) = lsf3_rr_rA(a,b,s,i) = d^3 S / (dr_a dr_b dR_(s,A))
!>
!> with the LSF's compact active index: slot `i` belongs to `active_atom(i)`.
!> This module supplies the directional derivative of that row along a nuclear
!> direction `v` -- the field row of the cavity Hessian.
!>
!> Two things move when the nuclei move: the adjoint weights, whose tangents
!> `(dw0, dw1, dw2)` the caller hands in, and the level-set jet itself, which
!> the LSF differentiates. The jet moves for two reasons at once, because the
!> row is evaluated at the *projected* point `r*(p)` and that point rides along:
!>
!>     d_v[T0]    = d^2S/(dR_A dR_B) v_B          + sum_k T1(k, ..) dr_k
!>     d_v[T1_a]  = d^3S/(dr_a dR_A dR_B) v_B     + sum_k T2(a, k, ..) dr_k
!>     d_v[T2_ab] = d^4S/(dr_a dr_b dR_A dR_B) v_B + sum_k T3(a, b, k, ..) dr_k
!>
!> The explicit halves are exactly `hvp_f1_rA`, `hvp_f2_r_rA` and
!> `hvp_f3_rr_rA`; the fourth-order spatial tensor is `f4_rrr_rA`.
!>
!> The folding identity
!> --------------------
!> Raising the spatial order of a term by one is precisely what `vjp_f1_rA`'s
!> `w1` and `w2` slots already do, so two of the three point-motion terms are
!> not separate contractions at all:
!>
!>     w0 * sum_k T1(k, s, i) dr_k
!>         = sum_a (w0 * dr_a) * lsf2_r_rA(a, s, i)
!>         = vjp_f1_rA(0, w0*dr, 0)
!>
!>     sum_a w1(a) sum_k T2(a, k, s, i) dr_k
!>         = sum_ab (w1(a) * dr(b)) * lsf3_rr_rA(a, b, s, i)
!>         = vjp_f1_rA(0, 0, outer(w1, dr))
!>
!> Both share the weight slots of the explicit weight-tangent term, so the
!> whole block collapses into a *single* call with shifted weights,
!>
!>     vjp_f1_rA(dw0, dw1 + w0*dr, dw2 + outer(w1, dr))
!>
!> where `outer(w1, dr)(a, b) = w1(a) * dr(b)`. That matrix is not symmetric,
!> and must not be symmetrised: `vjp_f1_rA` contracts all nine entries of `w2`
!> with no symmetry assumption and no folded factor of two, so any
!> redistribution across the diagonal would silently change the answer for a
!> caller whose own `dw2` is asymmetric. (`lsf3_rr_rA` happens to be symmetric
!> in `(a, b)`, so the ordering of the outer product is numerically inert
!> today; it is written in the order the contraction derives it in anyway.)
!>
!> The third point-motion term does *not* fold. `vjp_f1_rA` stops at two
!> spatial indices, and `w2_ab T3(a,b,k) dr_k` needs three, so it goes through
!> `f4_rrr_rA` directly. That is the one contraction here without a contracted
!> accessor, and the only reason this routine has to materialize a
!> `(3, 3, 3, 3, n_active)` tensor.
!>
!> Applicability
!> -------------
!> This is a primitive, not a driver: it neither reads cavity state nor
!> reconstructs the adjoints, and it does not know where `dr` came from. It
!> needs the full optional derivative set (`hvp_*` and `f4_rrr_rA`), which SvdW
!> and CFC implement and the isodensity LSFs do not -- passing one of those
!> aborts inside the accessor rather than returning a zero row, which is the
!> intended behaviour.
module moist_cavity_drop_derivatives_field_tangent
   use mctc_env_accuracy, only: wp
   use moist_cavity_drop_lsf_base, only: moist_cavity_drop_lsf_type

   implicit none(type, external)
   private

   public :: drop_field_tangent
   public :: drop_field_tangent_work_type

   !> Spatial dimension
   integer, parameter :: ndim = 3

   !> Scratch buffers of [[drop_field_tangent]]
   !>
   !> The three `hvp_*` accessors and `f4_rrr_rA` all return active-indexed
   !> tensors, the largest of them `(3, 3, 3, 3, n_active)`. Holding those as
   !> automatic arrays would put 120 doubles per active atom on the stack of
   !> every call, inside what is ultimately an OpenMP grid loop; holding them
   !> here instead makes the primitive allocation-free in steady state -- the
   !> buffers grow once and are then reused for every point of the grid.
   type :: drop_field_tangent_work_type
      !> Active slots the buffers are currently sized for
      integer :: capacity = 0
      !> `sum_B v_B . d^2S/(dR_A dR_B)` [3, capacity]
      real(wp), allocatable :: hvp1(:, :)
      !> `sum_B v_B . d^3S/(dr dR_A dR_B)` [3, 3, capacity]
      real(wp), allocatable :: hvp2(:, :, :)
      !> `sum_B v_B . d^4S/(dr^2 dR_A dR_B)` [3, 3, 3, capacity]
      real(wp), allocatable :: hvp3(:, :, :, :)
      !> `d^4S/(dr^3 dR_A)` [3, 3, 3, 3, capacity]
      real(wp), allocatable :: f4(:, :, :, :, :)
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

      if (allocated(self%hvp1)) deallocate (self%hvp1)
      if (allocated(self%hvp2)) deallocate (self%hvp2)
      if (allocated(self%hvp3)) deallocate (self%hvp3)
      if (allocated(self%f4)) deallocate (self%f4)

      allocate (self%hvp1(ndim, n))
      allocate (self%hvp2(ndim, ndim, n))
      allocate (self%hvp3(ndim, ndim, ndim, n))
      allocate (self%f4(ndim, ndim, ndim, ndim, n))
      self%capacity = n
   end subroutine field_tangent_work_ensure

   !* ================================================================================= *!
   !*                            Field-contraction tangent                              *!
   !* ================================================================================= *!

   !> Directional derivative of the jet-contracted nuclear-gradient row
   !>
   !> Returns `d_v` of `vjp_f1_rA(w0, w1, w2)` for the nuclear direction `v`,
   !> given the tangents `(dw0, dw1, dw2)` of the adjoint weights and the
   !> induced motion `dr` of the evaluation point. See the module header for
   !> the derivation and for the folding identity that turns the weight-tangent
   !> block and two of the three point-motion terms into a single `vjp_f1_rA`.
   !>
   !> `dr` is the caller's business. It is the tangent of the projected point
   !> `r*(p)` along `v` for the field row of the cavity Hessian, but nothing
   !> here assumes that; a caller wanting the partial derivative at a frozen
   !> evaluation point passes `dr = 0`.
   !>
   !> `res` follows `vjp_f1_rA`'s contract exactly, so the two are
   !> interchangeable at a call site: columns `1 .. active_count()` are
   !> *overwritten* with the tangent, columns beyond that are left untouched,
   !> and nothing at all is written when the active list is empty. The result
   !> is not accumulated into -- like `vjp_pt` in the gradient driver, the
   !> buffer is one point's row and the caller owns the scatter-add into the
   !> per-atom accumulator via `active_atom(i)`.
   !>
   !> The LSF must have been prepared at the evaluation point with a derivative
   !> order high enough for `f4_rrr_rA` (3 for SvdW, 4 for CFC); every accessor
   !> called here checks that itself and aborts if not.
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

      !> Shifted weights of the folded `vjp_f1_rA` call
      real(wp) :: w1_fold(ndim), w2_fold(ndim, ndim)
      !> Row accumulator and hoisted point-motion component
      real(wp) :: acc, drk
      !> Active slot, nuclear axis and spatial axes
      integer :: i, s, a, b, k
      !> Active slots of the prepared point
      integer :: n_active

      n_active = lsf%active_count()
      if (n_active == 0) return
      call work%ensure(n_active)

      !* ---------------------- Weight tangents and folded motion --------------------- *!

      ! `w0 * T1_k dr_k` rides in the gradient slot and `w1_a * T2_ak dr_k` in
      ! the Hessian slot, so both join the explicit weight tangents in one call.
      ! `outer(w1, dr)` is deliberately not symmetrised; see the module header.
      w1_fold = dw1 + w0*dr
      do b = 1, ndim
         do a = 1, ndim
            w2_fold(a, b) = dw2(a, b) + w1(a)*dr(b)
         end do
      end do
      call lsf%vjp_f1_rA(dw0, w1_fold, w2_fold, res)

      !* --------------------------- Explicit nuclear motion -------------------------- *!

      call lsf%hvp_f1_rA(v, work%hvp1)
      call lsf%hvp_f2_r_rA(v, work%hvp2)
      call lsf%hvp_f3_rr_rA(v, work%hvp3)

      !* ------------------------ Unfoldable point-motion term ------------------------ *!

      ! `w2_ab T3_abk dr_k` carries three spatial indices, one more than
      ! `vjp_f1_rA` can absorb, so the full mixed fourth derivative is formed
      ! and contracted here.
      call lsf%f4_rrr_rA(work%f4)

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
            do k = 1, ndim
               drk = dr(k)
               do b = 1, ndim
                  do a = 1, ndim
                     acc = acc + w2(a, b)*drk*work%f4(a, b, k, s, i)
                  end do
               end do
            end do
            res(s, i) = res(s, i) + acc
         end do
      end do
   end subroutine drop_field_tangent

end module moist_cavity_drop_derivatives_field_tangent
