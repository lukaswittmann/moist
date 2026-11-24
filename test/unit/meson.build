# One suite per file, flat: suite id <id> lives in test_<id>.f90 and defines
# module test_<id>. Keeping the three names in lockstep is what lets this list
# be the single source of truth -- the file name is derived, never spelled out.
#
#   tests             suite ids, passed to the tester as its argument
#   test_timeouts     suite id -> timeout, for suites needing more than the default
#   test_priority     suite id -> start priority, higher first; the slowest suites
#                      start at t=0 instead of queueing behind the short ones
#   test_selected     suite id -> test names, for a suite whose tests must each be
#                      invoked as `tester <suite> <test>` (test-drive's run_selected, no
#                      enclosing OpenMP team) rather than `tester <suite>`
#                      (run_testsuite, whose tests run inside test-drive's own
#                      `!$omp parallel do`; see CONTRIBUTING.md:47); also splits a
#                      long suite into separately scheduled tests. Every test of the
#                      suite must be listed, an omitted one is not run
#   test_chunks       suite id -> n, for a long suite run as `tester <suite> k/n`
#                      (every n-th test from the k-th, still in test-drive's
#                      OpenMP team); coverage runs single-threaded and needs the
#                      chunks to spread such a suite over the runner's cores
tests = [
  'api',
  'parameters',
  'utils_timer',
  'utils_context',
  'utils_mem',
  'utils_prettylistprint',
  'utils_prettylistprint_threaded',
  'utils_format',
  'utils_format_threaded',
  'channels',
  'data',
  'radii',
  'math_linalg',
  'math_smoothing',
  'math_adjacency',
  'math_solvers',
  'math_sorters',
  'math_trig',
  'math_packed',
  'math_smoothstep',
  'math_symmetry',
  'math_grid_radial_rule',
  'math_grid_radial',
  'math_grid_radial_trafo',
  'math_grid_angular',
  'math_grid_atomic',
  'math_fft_dst4',
  'math_fft_3d',
  'math_fft_threads',
  'finufft',
  'math_grid_3d_fft',
  'math_grid_3d_nufft',
  'math_grid_3d_threaded',
  'math_grid_3d_partition',
  'math_grid_3d_guards',
  'math_grid_3d_derivatives',
  'math_grid_3d_integration',
  'math_grid_3d_invariance',
  'moz_1d',
  'moz_3d',
  'moz_potential_kernel',
  'moz_potential_electrostatic',
  'moz_potential_ec',
  'moz_potential_lj',
  'moz_potential_oplsaa',
  'moz_solvent',
  'cavity_drop_phi',
  'cavity_drop_cfc',
  'cavity_drop_lsf',
  'cavity_drop_lsf_reference',
  'cavity_drop_lsf_screening',
  'cavity_drop_isodensity',
  'cavity_drop_gradient',
  'cavity_drop_branching',
  'cavity_drop_nuclear_adjoint',
  'cavity_drop_potential',
  'cavity_drop_integration',
  'cavity_drop_anchor',
  'cavity_drop_projector',
  'cavity_drop_projection',
  'cavity_drop_properties',
  'cavity_iswig',
  'cavity_numsa',
  'cavity_symmetry',
  'cavity_grid',
  'cavity_mc',
  'cavity_mt_mesh',
  'cavity_mt_patches',
  'cavity_mt',
  'cavity_reference',
  'cavity_adjoint',
  'model_component_pcm_amat_kernel',
  'model_component_pcm_amat_assembly',
  'model_component_pcm_amat_adjoint',
  'model_component_pcm_electrostatics',
  'model_component_pcm_solvers',
  'model_component_pcm_interface',
  'model_component_pcm_derivatives',
  'model_component_pcm_references',
  'model_component_pcm_regression',
  'model_component_pcm_vacuum',
  'model_component_pcm_kernel_variants',
  'model_component_gostshyp',
  'model_component_pv',
  'model_continuum',
  'model_coupling',
]
test_timeouts = {}

# Ranked by CI runtime
test_priority = {
  'cavity_drop_branching': 10,
  'cavity_drop_projection': 9,
  'model_component_pcm_derivatives': 8,
  'math_grid_3d_invariance': 7,
  'math_grid_3d_derivatives': 6,
  'cavity_drop_gradient': 5,
  'cavity_drop_lsf': 4,
  'cavity_drop_isodensity': 3,
  'model_component_pcm_references': 1,
}

test_chunks = {
  'cavity_drop_branching': 4,
  'model_component_pcm_derivatives': 4,
}

# threaded tests (single testdrive thread)
test_selected = {
  'math_grid_3d_integration': [
    'single_atom', 'asymmetric_trimer', 'partition_integrals', 'multicenter_integrals',
    'rotational_convergence', 'integrate_constant_is_box_volume', 'gaussian_width_law',
  ],
  'math_grid_3d_invariance': ['single_atom', 'asymmetric_trimer'],
  'math_grid_3d_derivatives': [
    'single_atom', 'dimer', 'amino20x4_gly_xab', 'mb16_43_01', 'but14diol_1', 'il16_008',
    'partition_schemes', 'gaussian_tiny_weights', 'volume_hessian_columns',
    'refused_calls_keep_accumulators', 'width_gradient_accepted',
    'bad_volume_adjoint_unbuilt',
    'bad_volume_adjoint_uninitialized',
    'bad_volume_adjoint_size',
    'bad_volume_adjoint_nonfinite',
    'bad_volume_width_adjoint',
    'bad_volume_gradient_shape',
    'bad_volume_hvp_accumulator_shape',
    'bad_volume_hvp_direction_shape',
    'bad_volume_direction_nonfinite',
    'bad_volume_hessian_shape',
    'bad_volume_hvp_width_curvature',
    'bad_volume_hessian_curvature',
    'bad_volume_phase_reference',
  ],
  'utils_prettylistprint_threaded': ['independent_rows'],
  'utils_format_threaded': ['concurrent_rendering'],
  'cavity_drop_projection': [
    'octree_connected_minima', 'octree_shallow_minimum', 'octree_search_window',
    'octree_rigid_motion', 'octree_zero_window_on_surface', 'octree_cluster_dimer',
    'octree_cluster_cross', 'octree_cluster_octahedron', 'octree_cluster_cube',
    'octree_cluster_pentagon', 'octree_cluster_tetrahedron', 'brute_force_broad_window',
    'stationary_box_bounds', 'batch_flush_accumulator', 'newton_sphere_point',
    'newton_stage_switch', 'newton_far_seed_unreached', 'newton_far_point',
    'polish_escapes_far_pole', 'multistart_capture_stops_duplicates',
    'multistart_capture_keeps_twin_points', 'finalize_canonical_order',
    'finalize_canonical_order_chain', 'finalize_sphere_anchor',
    'finalize_branch_weight_fold', 'octree_two_sphere_branches',
    'octree_single_sphere_one_branch', 'octree_sign_change_bound_is_tight',
    'octree_certifies_empty_ball', 'octree_surface_patch_coverage',
    'octree_optional_region_certificate', 'octree_seed_anchor_within_trace_tolerance',
    'octree_far_verified_point_radius', 'octree_far_traced_radius',
    'octree_seeds_on_sign_change_segments', 'octree_sign_free_region_no_seeds',
    'octree_refinement_thin_feature', 'bad_octree_seed_size', 'bad_octree_budgets',
    'bad_octree_max_depth', 'bad_octree_pdist_max', 'bad_octree_negative_slack',
    'bad_octree_depth_cap', 'bad_octree_survivors', 'bad_octree_box_budget',
    'bad_octree_seed_before_init', 'bad_octree_run_before_init',
    'bad_octree_refine_levels', 'octree_neighbours_outside_box_budget',
    'octree_refinement_outside_box_budget', 'worker_octree_refinement_by_level',
    'bad_batch_target_size', 'bad_newton_unreached_anchor',
  ],
  'math_grid_3d_threaded': [
    'batched_ifft_uses_threaded_backend',
    'molecular_nufft_type12_threaded',
    'molecular_nufft_type3_threaded',
    'radial_trafo_one_per_thread',
    'partition_grid_threaded',
  ],
  'cavity_drop_nuclear_adjoint': [
    'all_channels', 'single_channels', 'branching_channels', 'branching_xi_fd',
    'branching_f_fd', 'reach_gate_channels', 'forward_pair_channels',
    'sibling_gate_channels', 'forward_screened_neighbour', 'forward_screened_owner',
    'forward_isodensity_owner', 'forward_frozen_field',
    'model_forward_reverse', 'cosmo_model_forward_reverse', 'shape_guard',
    'thread_count_invariance', 'isodensity_thread_count_invariance',
  ],
}

# Driver plus the shared fixture modules, which register no suite of their own.
test_srcs = files(
  'main.F90',
  'test_helpers.f90',
  'test_cavity_drop_fd.f90',
  'test_pcm_conductor_probe.f90',
  'test_fixtures.f90',
  'test_model_component_helper.f90',
)

foreach t : tests
  test_srcs += files('test_' + t + '.f90')
endforeach

tester = executable(
  'tester',
  sources: test_srcs,
  dependencies: [moist_dep, mstore_dep, testdrive_dep],
  link_language: 'fortran',
)

foreach t : tests
  runs = {t: [t]}
  if test_chunks.has_key(t)
    runs = {}
    foreach k : range(1, test_chunks[t] + 1)
      runs += {'@0@_@1@of@2@'.format(t, k, test_chunks[t]): [t, '@0@/@1@'.format(k, test_chunks[t])]}
    endforeach
  elif test_selected.has_key(t)
    runs = {}
    foreach sel : test_selected[t]
      runs += {t + '_' + sel: [t, sel]}
    endforeach
  endif
  foreach name, test_args : runs
    test(
      name,
      tester,
      args: test_args,
      env: test_env,
      is_parallel : true,
      priority: test_priority.get(t, 0),
      timeout: test_timeouts.get(t, 900),
    )
  endforeach
endforeach
