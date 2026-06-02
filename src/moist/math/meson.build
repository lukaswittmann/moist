subdir('lapack')
subdir('blas')
subdir('linalg')
subdir('solver')
subdir('sorter')
subdir('grid')

srcs += files(
  'adjacency_list.f90',
  'cell_grid.f90',
  'lapack.f90',
  'blas.f90',
  'linalg.f90',
  'boys.f90',
  'solver.f90',
  'sorter.f90',
  'smoothing_kernels.f90',
  'packed_sym.f90',
  'trigonometry.f90',
  'grid.f90',
)

if use_fftw
  srcs += files('fftw3.F90')
endif
