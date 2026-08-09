!> Main PCM module - re-exports all PCM functionality
!> This is the top-level module users should import to access PCM functionality.
module moist_model_component_pcm
   use moist_model_component_pcm_type, only: pcm_base, &
      & pcm_solver_type, pcm_potential_source, &
      & solver_type, potential_source
   use moist_model_component_pcm_cpcm, only: cpcm, new_cpcm
   use moist_model_component_pcm_cosmo, only: cosmo, new_cosmo
   use moist_model_component_pcm_solvers, only: solve_pcm_lu, &
      & solve_pcm_cholesky, solve_pcm_iterative, solve_pcm_inversion
   use moist_model_component_pcm_amat, only: assemble_pcm_amat, &
      & assemble_pcm_amat_with_gradient, pcm_amat_surface_weights, &
      & pcm_amat_nuclear_gradient
   use moist_model_component_pcm_electrostatics, only: &
      & pcm_electrostatic_nuclear_gradient
   implicit none (type, external)
   public

end module moist_model_component_pcm
