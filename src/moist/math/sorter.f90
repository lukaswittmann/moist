!> Sorting routines
!>
!> - [[ord_sort]]: stable natural-run merge sort
!> - [[sort]]: unstable in-place introsort, nothing allocated
!> - [[sort_index]]: stable argsort; [[sort_adjoint]] moves a companion array
!> - [[counting_argsort]]: stable argsort of small integer keys, caller's counters
!>
!> The first three are ports of Fortran stdlib v0.8.1 for default-integer and
!> `real(wp)` keys, in `sorter/{merge,introsort,argsort}.f90` (MIT, see
!> `sorter/LICENSE`); `sorter/counting.f90` is moist's own
module moist_math_sorter
   use moist_math_sorter_merge, only: ord_sort
   use moist_math_sorter_introsort, only: sort
   use moist_math_sorter_argsort, only: sort_index, sort_adjoint
   use moist_math_sorter_counting, only: counting_argsort
   implicit none(type, external)
   private

   public :: ord_sort, sort, sort_index, sort_adjoint, counting_argsort
end module moist_math_sorter
