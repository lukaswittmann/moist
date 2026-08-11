test_srcs += files('helper.f90')

foreach leaf : ['gostshyp', 'pv']
  tests += 'model_component_' + leaf
  test_srcs += files(leaf + '.f90')
endforeach

subdir('pcm')
