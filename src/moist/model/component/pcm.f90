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
   implicit none
   public

end module moist_model_component_pcm
