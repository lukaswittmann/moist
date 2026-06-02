#pragma once

#ifdef __cplusplus
#define moist_API_ENTRY extern "C"
#else
#define moist_API_ENTRY extern
#ifndef moist_CFFI
#include <stdbool.h>
#endif
#endif
#define moist_API_CALL
#define moist_API_SUFFIX__V_0_5

/// Error handle class
typedef struct _moist_error* moist_error;

/// Molecular structure data class
typedef struct _moist_structure* moist_structure;

/// Solvation model class
typedef struct _moist_model* moist_model;

/// Damping parameter class
typedef struct _moist_param* moist_param;

/// DROP cavity class
typedef struct _moist_cavity* moist_cavity;

/// Radii model class
typedef struct _moist_radii* moist_radii;

/// Callback ABI for external isodensity DROP level-set functions.
/// The callback receives a point in Bohr and must return an LSF value,
/// spatial gradient, Hessian, and third derivative using the DROP sign convention
/// (interior negative, exterior positive).
typedef void (*moist_isodensity_lsf_callback)(void* /* context */,
                                              const double* /* point[3] */,
                                              double* /* value */,
                                              double* /* grad[3] */,
                                              double* /* hess[3][3] */,
                                              double* /* third[3][3][3] */);

/*
 * Type generic macro for convenience
**/

#define moist_delete(ptr) _Generic( \
        (ptr), \
        moist_error: moist_delete_error, \
        moist_structure: moist_delete_structure, \
        moist_model: moist_delete_solvation_model, \
        moist_param: moist_delete_param, \
        moist_cavity: moist_delete_drop_cavity, \
        moist_radii: moist_delete_radii \
    )(&ptr)

/*
 * Global API queries
**/

/// Obtain library version as major * 10000 + minor * 100 + patch
moist_API_ENTRY int moist_API_CALL
moist_get_version(void) moist_API_SUFFIX__V_0_5;

/// Get version string (e.g., "0.1.0"). Pass buffersize for bounded writes.
moist_API_ENTRY void moist_API_CALL
moist_get_version_string(char* /* buffer */,
                         const int* /* buffersize */) moist_API_SUFFIX__V_0_5;

/// Print MOIST header banner to file descriptor (use 6 for stdout, 0 for stderr)
moist_API_ENTRY void moist_API_CALL
moist_print_header(int /* unit */) moist_API_SUFFIX__V_0_5;

/// Print MOIST short header to file descriptor (use 6 for stdout, 0 for stderr)
moist_API_ENTRY void moist_API_CALL
moist_print_header_short(int /* unit */) moist_API_SUFFIX__V_0_5;

/// Print MOIST ASCII banner to file descriptor (use 6 for stdout, 0 for stderr)
moist_API_ENTRY void moist_API_CALL
moist_print_header_ascii(int /* unit */) moist_API_SUFFIX__V_0_5;

/// Print MOIST version to file descriptor
moist_API_ENTRY void moist_API_CALL
moist_print_version(int /* unit */) moist_API_SUFFIX__V_0_5;

/// Print build info (version, git commit, compiler, host) to file descriptor
moist_API_ENTRY void moist_API_CALL
moist_print_build_header(int /* unit */) moist_API_SUFFIX__V_0_5;

/// Print GEMS header banner to file descriptor (use 6 for stdout, 0 for stderr)
moist_API_ENTRY void moist_API_CALL
moist_print_gems_header(int /* unit */) moist_API_SUFFIX__V_0_5;

/*
 * Error handle class
**/

/// Create new error handle object
moist_API_ENTRY moist_error moist_API_CALL
moist_new_error(void) moist_API_SUFFIX__V_0_5;

/// Check error handle status
moist_API_ENTRY int moist_API_CALL
moist_check_error(moist_error /* error */) moist_API_SUFFIX__V_0_5;

/// Check error and exit with message if error is set
/// Does nothing if no error. Prints "[MOIST Error] context: message" to stderr and exits.
moist_API_ENTRY void moist_API_CALL
moist_check_error_exit(moist_error /* error */,
                       const char* /* context */,
                       const int* /* buffersize */) moist_API_SUFFIX__V_0_5;

/// Get error message from error handle. Pass buffersize for bounded writes.
moist_API_ENTRY void moist_API_CALL
moist_get_error(moist_error /* error */,
                char* /* buffer */,
                const int* /* buffersize */) moist_API_SUFFIX__V_0_5;

/// Delete error handle object
moist_API_ENTRY void moist_API_CALL
moist_delete_error(moist_error* /* error */) moist_API_SUFFIX__V_0_5;

/*
 * Molecular structure data class
**/

/// Create new molecular structure data (quantities in Bohr)
moist_API_ENTRY moist_structure moist_API_CALL
moist_new_structure(moist_error /* error */,
                    const int /* natoms */,
                    const int* /* numbers [natoms] */,
                    const double* /* positions [natoms][3] */,
                    const double* /* lattice [3][3] */,
                    const bool* /* periodic [3] */) moist_API_SUFFIX__V_0_5;

/// Delete molecular structure data
moist_API_ENTRY void moist_API_CALL
moist_delete_structure(moist_structure* /* mol */) moist_API_SUFFIX__V_0_5;

/// Update coordinates and lattice parameters (quantities in Bohr)
moist_API_ENTRY void moist_API_CALL
moist_update_structure(moist_error /* error */,
                       moist_structure /* mol */,
                       const double* /* positions [natoms][3] */,
                       const double* /* lattice [3][3] */) moist_API_SUFFIX__V_0_5;

/*
 * Radii model class
**/

/// Create CPCM radii model
moist_API_ENTRY moist_radii moist_API_CALL
moist_new_cpcm_radii(moist_error /* error */) moist_API_SUFFIX__V_0_5;

/// Create SMD radii model
moist_API_ENTRY moist_radii moist_API_CALL
moist_new_smd_radii(moist_error /* error */) moist_API_SUFFIX__V_0_5;

/// Create D3 radii model
moist_API_ENTRY moist_radii moist_API_CALL
moist_new_d3_radii(moist_error /* error */) moist_API_SUFFIX__V_0_5;

/// Create COSMO radii model
moist_API_ENTRY moist_radii moist_API_CALL
moist_new_cosmo_radii(moist_error /* error */) moist_API_SUFFIX__V_0_5;

/// Create Bondi radii model
moist_API_ENTRY moist_radii moist_API_CALL
moist_new_bondi_radii(moist_error /* error */) moist_API_SUFFIX__V_0_5;

/// Create custom radii model (must be populated before use)
moist_API_ENTRY moist_radii moist_API_CALL
moist_new_custom_radii(moist_error /* error */) moist_API_SUFFIX__V_0_5;

/// Set custom radii from per-atom values (bohr)
moist_API_ENTRY void moist_API_CALL
moist_set_custom_radii_atoms(moist_error /* error */,
                             moist_radii /* radii */,
                             const int /* natoms */,
                             const double* /* atom_radii[natoms] */) moist_API_SUFFIX__V_0_5;

/// Set custom radii from per-element values (bohr)
moist_API_ENTRY void moist_API_CALL
moist_set_custom_radii_elements(moist_error /* error */,
                                moist_radii /* radii */,
                                const int /* nentries */,
                                const int* /* atomic_numbers[nentries] */,
                                const double* /* element_radii[nentries] */) moist_API_SUFFIX__V_0_5;

/// Delete radii model
moist_API_ENTRY void moist_API_CALL
moist_delete_radii(moist_radii* /* radii */) moist_API_SUFFIX__V_0_5;

/*
 * Solvation model class
**/

// /// Create GEMS solvation model handle from a solvent name
// /// `parameter_file` may be NULL when `read_parameters` is false.
// moist_API_ENTRY moist_model moist_API_CALL
// moist_new_gems_solvation_model(moist_error /* error */,
//                                const char* /* solvent */,
//                                const bool* /* debug */,
//                                const int* /* verbosity */,
//                                const bool* /* read_parameters */,
//                                const char* /* parameter_file */) moist_API_SUFFIX__V_0_5;

// /// Update a solvation model with a molecular structure
// moist_API_ENTRY void moist_API_CALL
// moist_update_solvation_model(moist_error /* error */,
//                              moist_model /* model */,
//                              moist_structure /* mol */) moist_API_SUFFIX__V_0_5;

// /// Get total solvation energy from a solvation model
// moist_API_ENTRY void moist_API_CALL
// moist_get_solvation_model_energy(moist_error /* error */,
//                                  moist_model /* model */,
//                                  double* /* energy */) moist_API_SUFFIX__V_0_5;

// /// Get a borrowed cavity handle from a solvation model.
// /// The returned handle is valid as long as the parent model exists.
// /// Must call moist_update_solvation_model before extracting the cavity.
// /// Use moist_delete_cavity to release the handle (does NOT destroy the model's cavity).
// moist_API_ENTRY moist_cavity moist_API_CALL
// moist_get_solvation_model_cavity(moist_error /* error */,
//                                  moist_model /* model */) moist_API_SUFFIX__V_0_5;

// /// Delete solvation model handle
// moist_API_ENTRY void moist_API_CALL
// moist_delete_solvation_model(moist_model* /* model */) moist_API_SUFFIX__V_0_5;

/*
 * DROP cavity class
**/

/*
 * Type-specific constructors
**/

/// Create DROP cavity handle (does NOT build cavity - call moist_update_cavity after)
/// Optional: nleb (Lebedev grid size), debug (enable debug output), verbose (verbosity level 0-2),
///          blendk (blending k, default 2.0), blend1b (1-body weight, default 1.0),
///          blend2b (2-body weight, default 1.0), blend3b (3-body weight, default 1.0),
///          do_fine (enable all optional properties, default false)
/// Pass NULL for any optional parameter to use the default.
moist_API_ENTRY moist_cavity moist_API_CALL
moist_new_drop_cavity(moist_error /* error */,
                     const int* /* nleb */,
                     const bool* /* debug */,
                     const int* /* verbose */,
                     const double* /* blendk */,
                     const double* /* blend1b */,
                     const double* /* blend2b */,
                     const double* /* blend3b */,
                     const bool* /* do_fine */) moist_API_SUFFIX__V_0_5;

/// Create DROP cavity handle with explicit radii model (does NOT build cavity - call moist_update_cavity after)
/// Optional parameters same as moist_new_drop_cavity. Pass NULL for any to use the default.
moist_API_ENTRY moist_cavity moist_API_CALL
moist_new_drop_cavity_with_radii(moist_error /* error */,
                                moist_radii /* radii */,
                                const int* /* nleb */,
                                const bool* /* debug */,
                                const int* /* verbose */,
                                const double* /* blendk */,
                                const double* /* blend1b */,
                                const double* /* blend2b */,
                                const double* /* blend3b */,
                                const bool* /* do_fine */) moist_API_SUFFIX__V_0_5;

/// Create DROP cavity handle backed by an external isodensity LSF callback.
/// The callback must remain valid until the cavity is deleted. The context
/// pointer is passed back unchanged on every callback invocation.
moist_API_ENTRY moist_cavity moist_API_CALL
moist_new_drop_cavity_isodensity_callback(moist_error /* error */,
                                         moist_isodensity_lsf_callback /* callback */,
                                         void* /* context */,
                                         const double* /* scale */,
                                         const int* /* nleb */,
                                         const bool* /* debug */,
                                         const int* /* verbose */,
                                         const bool* /* do_fine */,
                                         const int* /* wleb_prune_level */) moist_API_SUFFIX__V_0_5;

/*
 * Generic cavity operations (Tier 1 - work on all cavity types)
**/

/// Generic update cavity - works for all cavity types
moist_API_ENTRY void moist_API_CALL
moist_update_cavity(moist_error /* error */,
                    moist_cavity /* cavity */,
                    moist_structure /* mol */) moist_API_SUFFIX__V_0_5;

/// Get generic cavity sizes - works for all cavity types
/// Returns ngrid (number of grid points) and nsph (number of spheres)
moist_API_ENTRY void moist_API_CALL
moist_get_cavity_sizes(moist_error /* error */,
                       moist_cavity /* cavity */,
                       int* /* ngrid */,
                       int* /* nsph */) moist_API_SUFFIX__V_0_5;

/// Get generic cavity results - works for all cavity types
/// Returns only fields from base cavity_type (area, volume, xyz, a, owner, converged, radii, asph)
/// Call moist_get_cavity_sizes first to get ngrid, nsph, then allocate arrays
moist_API_ENTRY void moist_API_CALL
moist_get_cavity_results(moist_error /* error */,
                         moist_cavity /* cavity */,
                         double* /* area */,
                         double* /* volume */,
                         int* /* ngrid */,
                         int* /* nsph */,
                         double* /* xyz[3][ngrid] */,
                         double* /* a[ngrid] */,
                         int* /* owner[ngrid] - 0-based atom indices (0 to nsph-1) */,
                         bool* /* converged[ngrid] */,
                         double* /* radii[nsph] */,
                         double* /* asph[nsph] */) moist_API_SUFFIX__V_0_5;

/// Generic delete cavity - works for all cavity types
moist_API_ENTRY void moist_API_CALL
moist_delete_cavity(moist_cavity* /* cavity */) moist_API_SUFFIX__V_0_5;

/*
 * Type-specific getters (Tier 2 - DROP-specific fields)
**/

/// Get DROP-specific cavity data (only works for DROP cavities)
/// Returns DROP-only fields (nmax, normal, wleb, r_iI0, f, rho)
/// Call moist_get_cavity_sizes first to get ngrid for array allocation
moist_API_ENTRY void moist_API_CALL
moist_get_drop_specific(moist_error /* error */,
                       moist_cavity /* cavity */,
                       int* /* nmax */,
                       double* /* normal[3][ngrid] */,
                       double* /* wleb[ngrid] */,
                       double* /* r_iI0[ngrid] */,
                       double* /* f[ngrid] */,
                       double* /* rho[ngrid] */) moist_API_SUFFIX__V_0_5;

/// Assemble A-matrix and compute xi values (DROP-specific)
/// Must be called before accessing xi or using the A-matrix
/// ngrid should be obtained from moist_get_cavity_sizes first
moist_API_ENTRY void moist_API_CALL
moist_assemble_amat(moist_error /* error */,
                    moist_cavity /* cavity */,
                    const int* /* ngrid */,
                    double* /* amat0[ngrid][ngrid] */,
                    double* /* xi[ngrid] */) moist_API_SUFFIX__V_0_5;

/*
 * Gradient API (Tier 3 - Cavity and A-matrix gradients)
**/

/// Compute cavity gradient w.r.t. nuclear coordinates (DROP-specific)
/// Must be called after moist_update_cavity and before moist_get_cavity_gradient
moist_API_ENTRY void moist_API_CALL
moist_compute_cavity_gradient(moist_error /* error */,
                              moist_cavity /* cavity */) moist_API_SUFFIX__V_0_5;

/// Get cavity gradient arrays (DROP-specific)
/// Must call moist_compute_cavity_gradient first
/// Call moist_get_cavity_sizes first to get ngrid, nsph for array allocation
/// Arrays:
///   A_tot1_rA[3][nsph]           - gradient of total area w.r.t. nuclear coords
///   V_tot1_rA[3][nsph]           - gradient of total volume w.r.t. nuclear coords
///   asph1_rA[3][nsph][nsph]      - gradient of per-sphere areas (owner, perturbed atom)
///   vsph1_rA[3][nsph][nsph]      - gradient of per-sphere volumes (owner, perturbed atom)
///   xyz1_rA[3][3][nsph][ngrid]   - grid point position derivatives (xyz, perturbed_xyz, atom, grid)
///   r_iI1_rA[3][nsph][ngrid]     - gradient of grid-owner distances
///   rho1_rA[3][nsph][ngrid]      - gradient of rho (anchor-to-surface distance)
moist_API_ENTRY void moist_API_CALL
moist_get_cavity_gradient(moist_error /* error */,
                          moist_cavity /* cavity */,
                          const int* /* nsph */,
                          const int* /* ngrid */,
                          double* /* A_tot1_rA[3][nsph] */,
                          double* /* V_tot1_rA[3][nsph] */,
                          double* /* asph1_rA[3][nsph][nsph] */,
                          double* /* vsph1_rA[3][nsph][nsph] */,
                          double* /* xyz1_rA[3][3][nsph][ngrid] */,
                          double* /* r_iI1_rA[3][nsph][ngrid] */,
                          double* /* rho1_rA[3][nsph][ngrid] */) moist_API_SUFFIX__V_0_5;

/// Assemble A-matrix with gradient (DROP-specific)
/// Must call moist_compute_cavity_gradient first
/// Call moist_get_cavity_sizes first to get ngrid, nsph for array allocation
/// Arrays:
///   Amat0[ngrid][ngrid]                - CPCM A-matrix
///   Amat1_rA[3][nsph][ngrid][ngrid]    - gradient of A-matrix w.r.t. nuclear coords
///   xi[ngrid]                          - xi values (screening factors)
moist_API_ENTRY void moist_API_CALL
moist_get_amat_gradient(moist_error /* error */,
                        moist_cavity /* cavity */,
                        const int* /* nsph */,
                        const int* /* ngrid */,
                        double* /* Amat0[ngrid][ngrid] */,
                        double* /* Amat1_rA[3][nsph][ngrid][ngrid] */,
                        double* /* xi[ngrid] */) moist_API_SUFFIX__V_0_5;

/// Contract A-matrix gradient with two grid vectors (DROP-specific)
/// Computes grad_rA = sum_ij q1_i * (dA_ij/dR_A) * q2_j
/// Must call moist_compute_cavity_gradient first
/// Uses cavity-internal ngrid and nsph for array extents
/// Arrays:
///   q1[ngrid], q2[ngrid]               - contraction vectors
///   grad_rA[3][nsph]                   - contracted nuclear gradient contribution
moist_API_ENTRY void moist_API_CALL
moist_contract_amat1_q1q2_rA(moist_error /* error */,
                             moist_cavity /* cavity */,
                             const double* /* q1[ngrid] */,
                             const double* /* q2[ngrid] */,
                             double* /* grad_rA[3][nsph] */) moist_API_SUFFIX__V_0_5;

/// Contract A-matrix derivatives to per-grid surface weights (DROP-specific)
/// Computes weights satisfying:
///   q1^T dA q2 = sum_i w_xi[i] dxi[i] + w_f[i] df[i] + w_xyz[:,i].dxyz[:,i]
/// Uses cavity-internal ngrid for array extents.
/// Arrays:
///   q1[ngrid], q2[ngrid]               - contraction vectors
///   w_xi[ngrid]                        - contracted xi derivative weights
///   w_f[ngrid]                         - contracted switch-function weights
///   w_xyz[3][ngrid]                    - contracted coordinate derivative weights
moist_API_ENTRY void moist_API_CALL
moist_contract_amat1_q1q2_surface_weights(moist_error /* error */,
                                          moist_cavity /* cavity */,
                                          const double* /* q1[ngrid] */,
                                          const double* /* q2[ngrid] */,
                                          double* /* w_xi[ngrid] */,
                                          double* /* w_f[ngrid] */,
                                          double* /* w_xyz[3][ngrid] */) moist_API_SUFFIX__V_0_5;

/// Contract DROP surface weights to per-grid LSF adjoint weights (DROP-specific)
/// Includes the projected-coordinate response from w_xyz and the electronic
/// xi path through wleb, cpjac, and the critical-gradient/focus switches.
/// The exported w_f is the anchor-only iSwiG overlap and has no electronic
/// LSF response for fixed nuclei.
/// Uses cavity-internal ngrid for array extents.
/// Arrays:
///   w_xi[ngrid]                        - xi derivative weights
///   w_f[ngrid]                         - anchor switch derivative weights
///   w_xyz[3][ngrid]                    - coordinate derivative weights
///   w_lsf0[ngrid]                      - LSF value adjoint weights
///   w_lsf1[3][ngrid]                   - LSF gradient adjoint weights
///   w_lsf2[3][3][ngrid]                - LSF Hessian adjoint weights
moist_API_ENTRY void moist_API_CALL
moist_contract_surface_lsf_weights(moist_error /* error */,
                                   moist_cavity /* cavity */,
                                   const double* /* w_xi[ngrid] */,
                                   const double* /* w_f[ngrid] */,
                                   const double* /* w_xyz[3][ngrid] */,
                                   double* /* w_lsf0[ngrid] */,
                                   double* /* w_lsf1[3][ngrid] */,
                                   double* /* w_lsf2[3][3][ngrid] */) moist_API_SUFFIX__V_0_5;

/// Contract nuclear + electronic CPCM terms (DROP-specific)
/// Must call moist_compute_cavity_gradient first
/// Uses cavity-internal ngrid and nsph for array extents
/// Arrays:
///   surface_q[ngrid]                   - surface charges q_i
///   qefield[3][ngrid]                  - electronic contribution Q_i * E_elec(i)
///   za[nsph]                           - nuclear charges Z_A
///   grad_rA[3][nsph]                   - contracted nuclear gradient contribution
moist_API_ENTRY void moist_API_CALL
moist_contract_nuc_elec_qefield_rA(moist_error /* error */,
                                   moist_cavity /* cavity */,
                                   const double* /* surface_q[ngrid] */,
                                   const double* /* qefield[3][ngrid] */,
                                   const double* /* za[nsph] */,
                                   double* /* grad_rA[3][nsph] */) moist_API_SUFFIX__V_0_5;

/*
 * Legacy API (deprecated - use generic versions instead)
**/

/// Update DROP cavity (legacy - use moist_update_cavity instead)
moist_API_ENTRY void moist_API_CALL
moist_update_drop_cavity(moist_error /* error */,
                        moist_cavity /* cavity */,
                        moist_structure /* mol */,
                        const int* /* nleb */) moist_API_SUFFIX__V_0_5;

/// Get DROP sizes (legacy - use moist_get_cavity_sizes and moist_get_drop_specific instead)
moist_API_ENTRY void moist_API_CALL
moist_get_drop_sizes(moist_error /* error */,
                    moist_cavity /* cavity */,
                    int* /* ngrid */,
                    int* /* nmax */,
                    int* /* nsph */) moist_API_SUFFIX__V_0_5;

/// Get DROP results (legacy - use moist_get_cavity_results and moist_get_drop_specific instead)
moist_API_ENTRY void moist_API_CALL
moist_get_drop_results(moist_error /* error */,
                      moist_cavity /* cavity */,
                      double* /* area */,
                      double* /* volume */,
                      int* /* ngrid */,
                      int* /* nmax */,
                      int* /* nsph */,
                      double* /* xyz[3][ngrid] */,
                      double* /* normal[3][ngrid] */,
                      double* /* wleb[ngrid] */,
                      double* /* a[ngrid] */,
                      double* /* r_iI0[ngrid] */,
                      double* /* f[ngrid] */,
                      double* /* rho[ngrid] */,
                      int* /* owner[ngrid] - 0-based atom indices (0 to nsph-1) */,
                      bool* /* converged[ngrid] */,
                      double* /* radii[nsph] */,
                      double* /* asph[nsph] */) moist_API_SUFFIX__V_0_5;

/// Delete DROP cavity handle (legacy - use moist_delete_cavity instead)
moist_API_ENTRY void moist_API_CALL
moist_delete_drop_cavity(moist_cavity* /* cavity */) moist_API_SUFFIX__V_0_5;
