module moist_math_sorter
   use moist_math_sorter_quicksort, only: qsort
   use moist_math_sorter_counting_sort, only: counting_argsort
   implicit none
   public :: qsort, counting_argsort
end module moist_math_sorter
