module moist_math_sorter
   use moist_math_sorter_quicksort, only: qsort
   use moist_math_sorter_counting, only: counting_argsort
   use moist_math_sorter_heapsort, only: heapsort
   implicit none(type, external)
   private

   public :: qsort, counting_argsort, heapsort
end module moist_math_sorter
