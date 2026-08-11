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

/// API version tags. Each declaration carries the tag of the release in which
/// its *current* form was introduced, so the list grows and older tags stay
/// defined. Because the tags are defined (to nothing) rather than omitted, a
/// consumer can feature-test with `#ifdef moist_API_SUFFIX__V_0_6` to tell
/// which contract this header describes -- useful where a declaration's C
/// signature is unchanged but the meaning of its arguments is not.
#define moist_API_SUFFIX__V_0_5
#define moist_API_SUFFIX__V_0_6

/*
 * ARRAY LAYOUT CONVENTION -- read this before allocating
 *
 * Every array crossing this API is a FLAT buffer in Fortran (column-major)
 * order. A parameter annotated `Fortran (d1,d2,...,dn)` holds d1*d2*...*dn
 * doubles (or ints/bools) and stores element (i1,i2,...,in) -- all indices
 * 0-based as seen from C -- at flat offset
 *
 *     i1 + d1*(i2 + d2*(i3 + ... + d_{n-1}*in))
 *
 * i.e. the FIRST dimension listed is the contiguous (fastest-varying) one.
 * Concretely, `Fortran (3,ngrid)` is ngrid points of 3 contiguous components:
 * component k of point i lives at `xyz[3*i + k]`. `Fortran (3,nsph,ngrid)`
 * puts element (a,A,i) at `[a + 3*(A + nsph*i)]`.
 *
 * WARNING: this is the transpose of the C declaration `double xyz[3][ngrid]`.
 * Do NOT declare these buffers as C multidimensional arrays with the
 * dimensions in the order written here -- declare them flat (e.g.
 * `double xyz[3 * ngrid];`) and index with the formula above, or declare them
 * with the dimensions REVERSED (`double xyz[ngrid][3];`).
 *
 * Annotations of the form `[n]` are plain 1-D arrays of n elements and mean
 * exactly what they say.
**/

/*
 * ARRAY CAPACITY CONVENTION -- how the array-writing getters are bounded
 *
 * Every entry point that writes into caller-allocated arrays takes the
 * capacities the caller allocated them with as leading `const int` arguments
 * (`ngrid_cap`, `nsph_cap`, ...), passed BY VALUE. A capacity is checked
 * before a single element is written: if it is smaller than the value the
 * cavity would report through moist_get_cavity_sizes, the call sets an API
 * error and returns without touching any buffer.
 *
 * A capacity LARGER than required is accepted, so a host may allocate one
 * max-size buffer once and reuse it across geometry steps. In that case the
 * capacity -- not the current ngrid/nsph -- defines the array's dimensions for
 * the layout convention above: `xyz: Fortran (3,ngrid_cap)` is indexed with
 * stride 3 for ngrid_cap points, and only the leading ngrid of them are
 * written. The same holds for the higher-rank gradient arrays, where a
 * capacity dimension sets the stride of every dimension to its left.
**/

/// Error handle class
typedef struct _moist_error* moist_error;

/// Molecular structure data class
typedef struct _moist_structure* moist_structure;

/// Solvation model class
typedef struct _moist_model* moist_model;

/// Solvation-model component class
typedef struct _moist_component* moist_component;

/// Damping parameter class.
/// Reserved: no constructor or destructor is exported yet, so it is not part
/// of the moist_delete dispatch below.
typedef struct _moist_param* moist_param;

/// PCM linear-solver selection, as accepted by moist_new_cpcm_component.
/// Mirrors the Fortran `solver_type` enumerator in
/// src/moist/model/component/pcm/type.f90; the numeric values are ABI.
typedef enum {
    /// Explicit matrix inversion
    moist_pcm_solver_inversion = 1,
    /// LU factorization (LAPACK GETRF+GETRS)
    moist_pcm_solver_lu = 2,
    /// Cholesky factorization (SPD matrices)
    moist_pcm_solver_cholesky = 3,
    /// Iterative solve (preconditioned CG)
    moist_pcm_solver_iterative = 4
} moist_pcm_solver;

/// DROP cavity class
typedef struct _moist_cavity* moist_cavity;

/// Radii model class
typedef struct _moist_radii* moist_radii;

/// Callback ABI for external isodensity DROP level set functions.
/// The callback receives a point in Bohr and must write an LSF value and
/// spatial gradient using the DROP sign convention (interior negative, exterior
/// positive).  The `hess` and `third` pointers are optional: they are NULL when
/// the cavity does not need that derivative order (e.g. the value+gradient-only
/// projection phase), and a NULL pointer means the callback should skip
/// computing that derivative, not merely skip writing it.  Dereferencing a NULL
/// `hess`/`third` is a host bug, so branch on them before writing.
/// `hess` and `third` follow the Fortran layout convention above; both are
/// symmetric under any index permutation, so the storage order is immaterial
/// in practice.
///
/// RETURN VALUE -- the failure channel.  Return 0 after writing every requested
/// buffer.  Return any nonzero value to report that the evaluation failed; the
/// output buffers are then ignored.  moist stops calling the callback for the
/// rest of that cavity build, unwinds its parallel evaluation loops, and fails
/// the enclosing moist_update_cavity() with a moist_error naming the returned
/// status.  There is no other way to abort a build from inside a callback: a
/// host that cannot signal failure leaves moist building a cavity out of
/// whatever happened to be in the buffers.
///
/// The status value itself is opaque to moist and is echoed verbatim in the
/// error message, so hosts may use it to carry their own error codes.
///
/// V_0_6 CONTRACT.  Two things changed relative to V_0_5: `hess`/`third` became
/// nullable (see above), and the return type became `int`.  The `void` -> `int`
/// change is deliberately compiler-visible: a callback still written against the
/// V_0_5 signature no longer type-checks against this typedef, so every stale
/// implementation is a compile error at the point where it is registered rather
/// than a silent misbehaviour at run time.  Fixing that compile error is the
/// prompt to also audit the nullable `hess`/`third` handling.
///
/// See test_isodensity_callback_cavity() in test/api/example.c for a callback
/// written against this contract.
typedef int (*moist_isodensity_lsf_callback)(void* /* context */,
                                             const double* /* point[3] */,
                                             double* /* value */,
                                             double* /* grad[3] */,
                                             double* /* hess: Fortran (3,3), or NULL */,
                                             double* /* third: Fortran (3,3,3), or NULL */);

/*
 * Type generic macro for convenience
**/

#define moist_delete(ptr) _Generic( \
        (ptr), \
        moist_error: moist_delete_error, \
        moist_structure: moist_delete_structure, \
        moist_model: moist_delete_solvation_model, \
        moist_component: moist_delete_solvation_component, \
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

/*
 * Error handle class
**/

/// Create new error handle object
moist_API_ENTRY moist_error moist_API_CALL
moist_new_error(void) moist_API_SUFFIX__V_0_5;

/// Check error handle status
moist_API_ENTRY int moist_API_CALL
moist_check_error(moist_error /* error */) moist_API_SUFFIX__V_0_5;

/* NOTE: there is deliberately no "check error and abort" entry point here.
 * A library must never terminate its host process: moist is linked into
 * quantum-chemistry programs, where an abort would kill an SCF mid-run, skip
 * the host's cleanup, and leave MPI/OpenMP state and scratch files behind --
 * with no way for the host to intercept it. Failure policy belongs to the
 * caller, so query the handle with moist_check_error() and, if it is set,
 * retrieve the message with moist_get_error(). A driver that does want to die
 * on the first failure can express that in a few lines of its own:
 *
 *     static void die_on_error(moist_error error, const char* context)
 *     {
 *         if (!moist_check_error(error)) return;
 *         char message[512];
 *         const int message_len = (int)sizeof(message);
 *         moist_get_error(error, message, &message_len);
 *         fprintf(stderr, "[moist Error] %s: %s\n", context, message);
 *         exit(EXIT_FAILURE);
 *     }
 *
 * A long-lived host should instead unwind, report, and keep running.
**/

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
                    const int* /* numbers[natoms] */,
                    const double* /* positions: Fortran (3,natoms) */,
                    const double* /* lattice: Fortran (3,3), vectors in columns */,
                    const bool* /* periodic[3] */) moist_API_SUFFIX__V_0_5;

/// Delete molecular structure data
moist_API_ENTRY void moist_API_CALL
moist_delete_structure(moist_structure* /* mol */) moist_API_SUFFIX__V_0_5;

/// Update coordinates and lattice parameters (quantities in Bohr)
moist_API_ENTRY void moist_API_CALL
moist_update_structure(moist_error /* error */,
                       moist_structure /* mol */,
                       const double* /* positions: Fortran (3,natoms) */,
                       const double* /* lattice: Fortran (3,3), vectors in columns */) moist_API_SUFFIX__V_0_5;

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

/// Create a CPCM component for a general solvation model.
/// `solver` uses the moist PCM solver enumeration (moist_pcm_solver).
moist_API_ENTRY moist_component moist_API_CALL
moist_new_cpcm_component(moist_error /* error */,
                         double /* epsilon */,
                         int /* solver */) moist_API_SUFFIX__V_0_6;

/// Create a pressure-volume energy component equal to `pressure * cavity volume`.
moist_API_ENTRY moist_component moist_API_CALL
moist_new_pv_component(moist_error /* error */,
                       double /* pressure */) moist_API_SUFFIX__V_0_6;

/// Create a GOSTSHYP hydrostatic-pressure component at `pressure` in
/// Hartree/bohr^3. The component cannot form its own density traces: supply the
/// Gaussian moments with moist_general_model_supply_gostshyp after every cavity
/// update, and read the amplitudes back with
/// moist_general_model_get_potential_extended.
moist_API_ENTRY moist_component moist_API_CALL
moist_new_gostshyp_component(moist_error /* error */,
                             double /* pressure */) moist_API_SUFFIX__V_0_6;

/// Delete a standalone solvation-model component handle.
moist_API_ENTRY void moist_API_CALL
moist_delete_solvation_component(moist_component* /* component */) moist_API_SUFFIX__V_0_6;

/// Create a general solvation model owning a copy of `cavity`.
/// Components can be appended until the model is first updated.
moist_API_ENTRY moist_model moist_API_CALL
moist_new_general_solvation_model(moist_error /* error */,
                                  moist_cavity /* cavity */,
                                  bool /* debug */,
                                  int /* verbosity */) moist_API_SUFFIX__V_0_6;

/// Append a copied component to a general model before its first update.
moist_API_ENTRY void moist_API_CALL
moist_general_model_add_component(moist_error /* error */,
                                  moist_model /* model */,
                                  moist_component /* component */) moist_API_SUFFIX__V_0_6;

/// Supply the external electrostatic potential and optional response arrays.
/// All arrays use the model cavity's native grid order. `phi` is required;
/// the other pointers may be NULL. Vector arrays are Fortran (3,ngrid), i.e.
/// component k of point i at `[3*i + k]`.
moist_API_ENTRY void moist_API_CALL
moist_general_model_supply_electrostatics(
    moist_error /* error */, moist_model /* model */, int /* ngrid */,
    const double* /* phi[ngrid] */, const double* /* w_xi[ngrid] or NULL */,
    const double* /* w_f[ngrid] or NULL */,
    const double* /* w_xyz: Fortran (3,ngrid), or NULL */,
    const double* /* w_n: Fortran (3,ngrid), or NULL */,
    const double* /* qefield: Fortran (3,ngrid), or NULL */) moist_API_SUFFIX__V_0_6;

/// Supply the Gaussian density moments the GOSTSHYP component consumes.
/// All arrays use the model cavity's native grid order and are required.
/// They are the moments of the solute density against the unnormalized Gaussian
/// `exp(-w_i |r - r_i|^2)` centered on each grid point: `gt = <G>`,
/// `pt = <(r-r_i) G>`, `mt = <(r-r_i)(r-r_i) G>`, `rt = <(r-r_i) |r-r_i|^2 G>`.
/// The width `w_i` is the model's own; read the grid-point areas back from the
/// live cavity before forming them.
///
/// CALLER OBLIGATION: rebuild and re-supply after EVERY cavity update.
moist_API_ENTRY void moist_API_CALL
moist_general_model_supply_gostshyp(
    moist_error /* error */, moist_model /* model */, int /* ngrid */,
    const double* /* gt[ngrid] */,
    const double* /* pt: Fortran (3,ngrid) */,
    const double* /* mt: Fortran (3,3,ngrid) */,
    const double* /* rt: Fortran (3,ngrid) */) moist_API_SUFFIX__V_0_6;

/// Return direct molecular-trace adjoints from all model components.
moist_API_ENTRY void moist_API_CALL
moist_general_model_get_trace_potential(
    moist_error /* error */, moist_model /* model */, int /* ngrid */,
    double* /* w_umol[ngrid] */, double* /* w_qmol[ngrid] */) moist_API_SUFFIX__V_0_6;

/// Return direct trace and level-set adjoints from all model components.
moist_API_ENTRY void moist_API_CALL
moist_general_model_get_potential(
    moist_error /* error */, moist_model /* model */, int /* ngrid */,
    double* /* w_umol[ngrid] */, double* /* w_qmol[ngrid] */,
    double* /* w_lsf0[ngrid] */, double* /* w_lsf1: Fortran (3,ngrid) */,
    double* /* w_lsf2: Fortran (3,3,ngrid) */) moist_API_SUFFIX__V_0_6;

/// Return every potential channel, including the GOSTSHYP host amplitudes.
/// Identical to moist_general_model_get_potential but also reports the
/// amplitudes conjugate to the host's Gaussian integral blocks, which the host
/// contracts as `F_uv += sum_i [w_gauss_g[i] g_uv,i + w_gauss_f[i] f_uv,i]`.
/// Both amplitude pointers may be NULL, and both are zero-filled when no
/// component supplies them, so a model without GOSTSHYP is not an error.
/// Prefer this over two calls: assembling a potential contracts the cavity
/// surface adjoints once, and splitting the read would pay that cost twice.
moist_API_ENTRY void moist_API_CALL
moist_general_model_get_potential_extended(
    moist_error /* error */, moist_model /* model */, int /* ngrid */,
    double* /* w_umol[ngrid] */, double* /* w_qmol[ngrid] */,
    double* /* w_lsf0[ngrid] */, double* /* w_lsf1: Fortran (3,ngrid) */,
    double* /* w_lsf2: Fortran (3,3,ngrid) */,
    double* /* w_gauss_g[ngrid] or NULL */,
    double* /* w_gauss_f[ngrid] or NULL */) moist_API_SUFFIX__V_0_6;

/// Return the accumulated nuclear gradient from all model components.
moist_API_ENTRY void moist_API_CALL
moist_general_model_get_gradient(moist_error /* error */,
                                 moist_model /* model */,
                                 int /* natoms */,
                                 double* /* gradient: Fortran (3,natoms) */) moist_API_SUFFIX__V_0_6;

/// Update a solvation model with a molecular structure
moist_API_ENTRY void moist_API_CALL
moist_update_solvation_model(moist_error /* error */,
                             moist_model /* model */,
                             moist_structure /* mol */) moist_API_SUFFIX__V_0_5;

/// Get total solvation energy from a solvation model
moist_API_ENTRY void moist_API_CALL
moist_get_solvation_model_energy(moist_error /* error */,
                                 moist_model /* model */,
                                 double* /* energy */) moist_API_SUFFIX__V_0_5;

/// Get a borrowed cavity handle from a solvation model.
/// The returned handle is valid as long as the parent model exists, but cannot
/// be rebuilt independently: moist_update_cavity and moist_update_drop_cavity
/// reject borrowed handles.
/// Use moist_delete_cavity to release the handle (does NOT destroy the model's cavity).
moist_API_ENTRY moist_cavity moist_API_CALL
moist_get_solvation_model_cavity(moist_error /* error */,
                                 moist_model /* model */) moist_API_SUFFIX__V_0_5;

/// Delete solvation model handle
moist_API_ENTRY void moist_API_CALL
moist_delete_solvation_model(moist_model* /* model */) moist_API_SUFFIX__V_0_5;

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
///          do_fine (enable all optional properties, default false),
///          tolerance (master numerical tolerance, default from the DROP parameters)
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
                     const bool* /* do_fine */,
                     const double* /* tolerance */) moist_API_SUFFIX__V_0_6;

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
                                const bool* /* do_fine */,
                                const double* /* tolerance */) moist_API_SUFFIX__V_0_6;

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
                                         const int* /* wleb_prune_level */,
                                         const double* /* tolerance */) moist_API_SUFFIX__V_0_6;

/// Create a DROP cavity backed by the internal isodensity LSF evaluator.
/// The atom-centered cartesian-monomial Gaussian basis is supplied here; the
/// (transformed) density matrix is installed each SCF step via
/// moist_set_isodensity_density. `shell_atom` uses 0-based atom indices; `exps`
/// and `coeffs` are concatenated per shell (length sum(shell_nprim)). `coeffs`
/// are the host's primitive-normalized contraction coefficients, so that
/// monomial * sum_p coeff_p exp(-a_p r^2) reproduces the host radial part; the
/// host maps its density into moist's cartesian-monomial basis (see
/// moist_get_isodensity_cart_layout). scale/nleb/debug/verbose/do_fine/
/// wleb_prune_level/tolerance are optional (NULL -> default).
moist_API_ENTRY moist_cavity moist_API_CALL
moist_new_drop_cavity_isodensity_internal(moist_error /* error */,
                                          const int /* nshell */,
                                          const int* /* shell_atom[nshell] (0-based) */,
                                          const int* /* shell_l[nshell] */,
                                          const int* /* shell_nprim[nshell] */,
                                          const double* /* exps[nprim_tot] */,
                                          const double* /* coeffs[nprim_tot] */,
                                          const double /* rho_iso */,
                                          const double* /* scale */,
                                          const int* /* nleb */,
                                          const bool* /* debug */,
                                          const int* /* verbose */,
                                          const bool* /* do_fine */,
                                          const int* /* wleb_prune_level */,
                                          const double* /* tolerance */) moist_API_SUFFIX__V_0_6;

/// Query the internal isodensity cartesian-component layout (DROP-specific).
/// Two-pass: call with the array pointers NULL to read ncart/nshell, then
/// allocate and call again. `shell_cart_offset` is 0-based (length nshell+1);
/// `comp_lx/ly/lz` are the per-component monomial powers (length ncart) defining
/// the ordering the host must match when building its density transform. Any of
/// the output pointers may be NULL to skip it.
moist_API_ENTRY void moist_API_CALL
moist_get_isodensity_cart_layout(moist_error /* error */,
                                 moist_cavity /* cavity */,
                                 int* /* ncart */,
                                 int* /* nshell */,
                                 int* /* shell_cart_offset[nshell+1] */,
                                 int* /* comp_lx[ncart] */,
                                 int* /* comp_ly[ncart] */,
                                 int* /* comp_lz[ncart] */) moist_API_SUFFIX__V_0_6;

/// Install the cartesian-monomial density matrix for the internal isodensity LSF.
/// `dcart` is the ncart-by-ncart density matrix, Fortran (ncart,ncart): element
/// (p,q) at `[p + ncart*q]`. Call before each moist_update_cavity.
moist_API_ENTRY void moist_API_CALL
moist_set_isodensity_density(moist_error /* error */,
                             moist_cavity /* cavity */,
                             const int /* ncart */,
                             const double* /* dcart: Fortran (ncart,ncart) */) moist_API_SUFFIX__V_0_6;

/// Return the master numerical tolerance currently configured on a DROP cavity.
/// Reports an error for non-DROP cavity types.
moist_API_ENTRY void moist_API_CALL
moist_get_drop_cavity_tolerance(moist_error /* error */,
                                moist_cavity /* cavity */,
                                double* /* tolerance */) moist_API_SUFFIX__V_0_6;

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
/// Returns fields shared by all cavity types. `converged` contains DROP's
/// per-point projection status; cavity types without a projection solve return
/// true for every retained point.
/// Call moist_get_cavity_sizes first to get ngrid, nsph, then allocate arrays
/// and pass the capacities you allocated them with (see the capacity
/// convention at the top of this header). Nothing is written when either
/// capacity is too small; an API error is set instead.
moist_API_ENTRY void moist_API_CALL
moist_get_cavity_results(moist_error /* error */,
                         moist_cavity /* cavity */,
                         const int /* ngrid_cap: allocated grid capacity,
                                      must be >= the ngrid reported by
                                      moist_get_cavity_sizes */,
                         const int /* nsph_cap: allocated sphere capacity,
                                      must be >= the nsph reported by
                                      moist_get_cavity_sizes */,
                         double* /* area */,
                         double* /* volume */,
                         int* /* ngrid */,
                         int* /* nsph */,
                         double* /* xyz: Fortran (3,ngrid_cap) */,
                         double* /* a[ngrid_cap] */,
                         int* /* owner[ngrid_cap] - 0-based atom indices (0 to nsph-1) */,
                         bool* /* converged[ngrid_cap] */,
                         double* /* radii[nsph_cap] */,
                         double* /* asph[nsph_cap] */) moist_API_SUFFIX__V_0_6;

/// Generic delete cavity - works for all cavity types
moist_API_ENTRY void moist_API_CALL
moist_delete_cavity(moist_cavity* /* cavity */) moist_API_SUFFIX__V_0_5;

/*
 * Type-specific getters (Tier 2 - DROP-specific fields)
**/

/// Get DROP-specific cavity data (only works for DROP cavities)
/// Returns DROP-only fields (nmax, normal, wleb, r_iI0, f, rho)
/// Call moist_get_cavity_sizes first to get ngrid for array allocation, then
/// pass the capacity you allocated with (see the capacity convention at the
/// top of this header). Nothing is written when the capacity is too small.
moist_API_ENTRY void moist_API_CALL
moist_get_drop_specific(moist_error /* error */,
                       moist_cavity /* cavity */,
                       const int /* ngrid_cap: allocated grid capacity, must be
                                    >= the ngrid reported by
                                    moist_get_cavity_sizes */,
                       int* /* nmax */,
                       double* /* normal: Fortran (3,ngrid_cap) */,
                       double* /* wleb[ngrid_cap] */,
                       double* /* r_iI0[ngrid_cap] */,
                       double* /* f[ngrid_cap] */,
                       double* /* rho[ngrid_cap] */) moist_API_SUFFIX__V_0_6;

/// Get the stable per-grid point numbering of a DROP cavity (DROP-specific)
/// The numbering uniquely identifies the same physical surface point across
/// cavity rebuilds, so callers can match points between successive cavities.
/// Call moist_get_cavity_sizes first to get ngrid for array allocation, then
/// pass the capacity you allocated with. Nothing is written when it is too small.
moist_API_ENTRY void moist_API_CALL
moist_get_drop_numbering(moist_error /* error */,
                         moist_cavity /* cavity */,
                         const int /* ngrid_cap: allocated grid capacity, must
                                      be >= the ngrid reported by
                                      moist_get_cavity_sizes */,
                         int* /* numbering[ngrid_cap] */) moist_API_SUFFIX__V_0_6;

/// Assemble A-matrix and compute xi values
/// Must be called before accessing xi or using the A-matrix
/// Call moist_get_cavity_sizes first to get ngrid for array allocation, then
/// pass the capacity you allocated with. Nothing is written when it is too small.
/// Works for every Gaussian-discretized cavity (no longer DROP-specific)
moist_API_ENTRY void moist_API_CALL
moist_assemble_amat(moist_error /* error */,
                    moist_cavity /* cavity */,
                    const int /* ngrid_cap: allocated grid capacity, must be >=
                                 the ngrid reported by moist_get_cavity_sizes */,
                    double* /* amat0: Fortran (ngrid_cap,ngrid_cap); symmetric,
                               so the transpose reading gives the same buffer */,
                    double* /* xi[ngrid_cap] */) moist_API_SUFFIX__V_0_6;

/// Return Gaussian PCM widths and switching factors.
/// Call moist_get_cavity_sizes first and pass the capacity allocated for both
/// arrays. Nothing is written when it is too small.
moist_API_ENTRY void moist_API_CALL
moist_get_cavity_gaussian(moist_error /* error */,
                          moist_cavity /* cavity */,
                          const int /* ngrid_cap: allocated grid capacity, must
                                       be >= the ngrid reported by
                                       moist_get_cavity_sizes */,
                          double* /* xi[ngrid_cap] */,
                          double* /* f[ngrid_cap] */) moist_API_SUFFIX__V_0_6;

/*
 * Gradient API (Tier 3 - Cavity and A-matrix gradients)
**/

/// Compute cavity gradient w.r.t. nuclear coordinates
/// Must be called after moist_update_cavity and before moist_get_cavity_gradient
moist_API_ENTRY void moist_API_CALL
moist_compute_cavity_gradient(moist_error /* error */,
                              moist_cavity /* cavity */) moist_API_SUFFIX__V_0_5;

/// Compute anchor-only nuclear derivatives (callback/isodensity LSF, DROP-specific)
/// Must be called after moist_update_cavity and before the *_rA contractions.
/// Restricts each grid point's nuclear coupling to its owner atom's rigid anchor
/// motion (the level set field's nuclear derivatives are zero for callback LSFs).
moist_API_ENTRY void moist_API_CALL
moist_compute_anchor_gradient(moist_error /* error */,
                              moist_cavity /* cavity */) moist_API_SUFFIX__V_0_6;

/// Get anchor-channel nuclear derivatives from the anchor pass
/// Must call moist_compute_anchor_gradient (or moist_compute_cavity_gradient) first
/// Call moist_get_cavity_sizes first and pass the capacities you allocated the
/// arrays with (see the capacity convention at the top of this header).
/// Arrays (Fortran order; see the layout convention at the top of this header):
///   xyz1_rA  Fortran (3,3,nsph_cap,ngrid_cap) - d(r_i)_j / d(R_A)_alpha  (j, alpha, A, i)
///                                       element at [j + 3*(alpha + 3*(A + nsph_cap*i))]
///   xi1_rA   Fortran (3,nsph_cap,ngrid_cap)   - d(xi_i)  / d(R_A)_alpha  (alpha, A, i)
///                                       element at [alpha + 3*(A + nsph_cap*i)]
///   a_i1_rA  Fortran (3,nsph_cap,ngrid_cap)   - d(a_i)   / d(R_A)_alpha  (alpha, A, i)
///   v_i1_rA  Fortran (3,nsph_cap,ngrid_cap)   - d(v_i)   / d(R_A)_alpha  (alpha, A, i)
///   A_tot1_rA Fortran (3,nsph_cap)    - d(total area)   / d(R_A)_alpha  (alpha, A)
///   V_tot1_rA Fortran (3,nsph_cap)    - d(total volume) / d(R_A)_alpha  (alpha, A)
/// a_i1_rA/v_i1_rA are the un-summed per-point counterparts of A_tot1_rA/V_tot1_rA.
/// The grid point area carries a switching-function dependence (a_i ~ f_i/xi_i^2),
/// so a_i1_rA is not recoverable from xi1_rA alone; it is the area route a generic
/// geometric surface functional (e.g. GOSTSHYP) needs.
moist_API_ENTRY void moist_API_CALL
moist_get_anchor_gradient(moist_error /* error */,
                          moist_cavity /* cavity */,
                          const int /* nsph_cap: allocated sphere capacity, must
                                       be >= the nsph reported by
                                       moist_get_cavity_sizes */,
                          const int /* ngrid_cap: allocated grid capacity, must
                                       be >= the ngrid reported by
                                       moist_get_cavity_sizes */,
                          double* /* xyz1_rA: Fortran (3,3,nsph_cap,ngrid_cap) */,
                          double* /* xi1_rA: Fortran (3,nsph_cap,ngrid_cap) */,
                          double* /* a_i1_rA: Fortran (3,nsph_cap,ngrid_cap) */,
                          double* /* v_i1_rA: Fortran (3,nsph_cap,ngrid_cap) */,
                          double* /* A_tot1_rA: Fortran (3,nsph_cap) */,
                          double* /* V_tot1_rA: Fortran (3,nsph_cap) */) moist_API_SUFFIX__V_0_6;

/// Get cavity gradient arrays (DROP-specific)
/// Must call moist_compute_cavity_gradient first
/// Call moist_get_cavity_sizes first to get ngrid, nsph for array allocation,
/// then pass the capacities you allocated with (see the capacity convention at
/// the top of this header).
/// Arrays (Fortran order; see the layout convention at the top of this header):
///   A_tot1_rA Fortran (3,nsph_cap)    - gradient of total area w.r.t. nuclear coords
///   V_tot1_rA Fortran (3,nsph_cap)    - gradient of total volume w.r.t. nuclear coords
///   asph1_rA  Fortran (3,nsph_cap,nsph_cap) - gradient of per-sphere areas (xyz, owner, perturbed atom)
///   vsph1_rA  Fortran (3,nsph_cap,nsph_cap) - gradient of per-sphere volumes (xyz, owner, perturbed atom)
///   xyz1_rA   Fortran (3,3,nsph_cap,ngrid_cap) - grid point position derivatives (xyz, perturbed_xyz, atom, grid)
///   r_iI1_rA  Fortran (3,nsph_cap,ngrid_cap) - gradient of grid-owner distances
///   rho1_rA   Fortran (3,nsph_cap,ngrid_cap) - gradient of rho (anchor-to-surface distance)
moist_API_ENTRY void moist_API_CALL
moist_get_cavity_gradient(moist_error /* error */,
                          moist_cavity /* cavity */,
                          const int /* nsph_cap: allocated sphere capacity, must
                                       be >= the nsph reported by
                                       moist_get_cavity_sizes */,
                          const int /* ngrid_cap: allocated grid capacity, must
                                       be >= the ngrid reported by
                                       moist_get_cavity_sizes */,
                          double* /* A_tot1_rA: Fortran (3,nsph_cap) */,
                          double* /* V_tot1_rA: Fortran (3,nsph_cap) */,
                          double* /* asph1_rA: Fortran (3,nsph_cap,nsph_cap) */,
                          double* /* vsph1_rA: Fortran (3,nsph_cap,nsph_cap) */,
                          double* /* xyz1_rA: Fortran (3,3,nsph_cap,ngrid_cap) */,
                          double* /* r_iI1_rA: Fortran (3,nsph_cap,ngrid_cap) */,
                          double* /* rho1_rA: Fortran (3,nsph_cap,ngrid_cap) */) moist_API_SUFFIX__V_0_6;

/// Assemble a Gaussian PCM A-matrix with its nuclear derivatives.
/// Must call moist_compute_cavity_gradient first
/// Call moist_get_cavity_sizes first to get ngrid, nsph for array allocation,
/// then pass the capacities you allocated with (see the capacity convention at
/// the top of this header).
/// Arrays (Fortran order; see the layout convention at the top of this header):
///   Amat0     Fortran (ngrid_cap,ngrid_cap)  - CPCM A-matrix (symmetric, so the
///                                            transpose reading is identical)
///   Amat1_rA  Fortran (3,nsph_cap,ngrid_cap,ngrid_cap) - gradient of A-matrix
///                                            w.r.t. nuclear coords; element
///                                            (alpha,A,i,j) at
///                                            [alpha + 3*(A + nsph_cap*(i + ngrid_cap*j))]
///   xi        [ngrid_cap]                    - xi values (screening factors)
moist_API_ENTRY void moist_API_CALL
moist_get_amat_gradient(moist_error /* error */,
                        moist_cavity /* cavity */,
                        const int /* nsph_cap: allocated sphere capacity, must be
                                     >= the nsph reported by
                                     moist_get_cavity_sizes */,
                        const int /* ngrid_cap: allocated grid capacity, must be
                                     >= the ngrid reported by
                                     moist_get_cavity_sizes */,
                        double* /* Amat0: Fortran (ngrid_cap,ngrid_cap) */,
                        double* /* Amat1_rA: Fortran (3,nsph_cap,ngrid_cap,ngrid_cap) */,
                        double* /* xi[ngrid_cap] */) moist_API_SUFFIX__V_0_6;

/// Contract Gaussian PCM A-matrix derivatives with two grid vectors.
/// Computes grad_rA = sum_ij q1_i * (dA_ij/dR_A) * q2_j
/// Must call moist_compute_cavity_gradient first
/// Uses cavity-internal ngrid and nsph for array extents
/// Arrays (Fortran order; see the layout convention at the top of this header):
///   q1[ngrid], q2[ngrid]               - contraction vectors
///   grad_rA   Fortran (3,nsph)         - contracted nuclear gradient contribution
moist_API_ENTRY void moist_API_CALL
moist_contract_amat1_q1q2_rA(moist_error /* error */,
                             moist_cavity /* cavity */,
                             const double* /* q1[ngrid] */,
                             const double* /* q2[ngrid] */,
                             double* /* grad_rA: Fortran (3,nsph) */) moist_API_SUFFIX__V_0_5;

/// Contract Gaussian PCM A-matrix derivatives to per-grid surface weights.
/// Computes weights satisfying:
///   q1^T dA q2 = sum_i w_xi[i] dxi[i] + w_f[i] df[i] + w_xyz[:,i].dxyz[:,i]
/// Uses cavity-internal ngrid for array extents.
/// Arrays (Fortran order; see the layout convention at the top of this header):
///   q1[ngrid], q2[ngrid]               - contraction vectors
///   w_xi[ngrid]                        - contracted xi derivative weights
///   w_f[ngrid]                         - contracted switch-function weights
///   w_xyz     Fortran (3,ngrid)        - contracted coordinate derivative weights
moist_API_ENTRY void moist_API_CALL
moist_contract_amat1_q1q2_surface_weights(moist_error /* error */,
                                          moist_cavity /* cavity */,
                                          const double* /* q1[ngrid] */,
                                          const double* /* q2[ngrid] */,
                                          double* /* w_xi[ngrid] */,
                                          double* /* w_f[ngrid] */,
                                          double* /* w_xyz: Fortran (3,ngrid) */) moist_API_SUFFIX__V_0_5;

/// Contract DROP surface weights to per-grid LSF adjoint weights (DROP-specific)
/// Includes the projected-coordinate response from w_xyz and the electronic
/// xi path through wleb, cpjac, and the critical-gradient/focus switches.
/// The exported w_f is the anchor-only iSwiG overlap and has no electronic
/// LSF response for fixed nuclei.
/// Uses cavity-internal ngrid for array extents.
/// Arrays (Fortran order; see the layout convention at the top of this header):
///   w_xi[ngrid]                        - xi derivative weights
///   w_f[ngrid]                         - anchor switch derivative weights
///   w_xyz     Fortran (3,ngrid)        - coordinate derivative weights
///   w_lsf0[ngrid]                      - LSF value adjoint weights
///   w_lsf1    Fortran (3,ngrid)        - LSF gradient adjoint weights
///   w_lsf2    Fortran (3,3,ngrid)      - LSF Hessian adjoint weights
/// This entry point covers the xi/f/xyz channels only; it is a thin wrapper
/// that forwards NULL for the extended channels below.
moist_API_ENTRY void moist_API_CALL
moist_contract_surface_lsf_weights(moist_error /* error */,
                                   moist_cavity /* cavity */,
                                   const double* /* w_xi[ngrid] */,
                                   const double* /* w_f[ngrid] */,
                                   const double* /* w_xyz: Fortran (3,ngrid) */,
                                   double* /* w_lsf0[ngrid] */,
                                   double* /* w_lsf1: Fortran (3,ngrid) */,
                                   double* /* w_lsf2: Fortran (3,3,ngrid) */) moist_API_SUFFIX__V_0_5;

/// Extended DROP surface-to-LSF contraction with optional normal/curvature channels.
/// Same contract as moist_contract_surface_lsf_weights plus:
///   w_n       Fortran (3,ngrid)        - outward-normal derivative weights (or NULL)
///   w_k1[ngrid]                        - first principal curvature weights (or NULL)
///   w_k2[ngrid]                        - second principal curvature weights (or NULL)
/// Pass NULL for any optional channel that should be skipped.
moist_API_ENTRY void moist_API_CALL
moist_contract_surface_lsf_weights_extended(
                                   moist_error /* error */,
                                   moist_cavity /* cavity */,
                                   const double* /* w_xi[ngrid] */,
                                   const double* /* w_f[ngrid] */,
                                   const double* /* w_xyz: Fortran (3,ngrid) */,
                                   double* /* w_lsf0[ngrid] */,
                                   double* /* w_lsf1: Fortran (3,ngrid) */,
                                   double* /* w_lsf2: Fortran (3,3,ngrid) */,
                                   const double* /* w_n: Fortran (3,ngrid), or NULL */,
                                   const double* /* w_k1[ngrid] or NULL */,
                                   const double* /* w_k2[ngrid] or NULL */) moist_API_SUFFIX__V_0_6;

/// Contract direct nuclear and electronic PCM terms.
/// Must call moist_compute_cavity_gradient first
/// Uses cavity-internal ngrid and nsph for array extents
/// Arrays (Fortran order; see the layout convention at the top of this header):
///   surface_q[ngrid]                   - surface charges q_i
///   qefield   Fortran (3,ngrid)        - electronic contribution Q_i * E_elec(i)
///   za[nsph]                           - nuclear charges Z_A
///   grad_rA   Fortran (3,nsph)         - contracted nuclear gradient contribution
moist_API_ENTRY void moist_API_CALL
moist_contract_nuc_elec_qefield_rA(moist_error /* error */,
                                   moist_cavity /* cavity */,
                                   const double* /* surface_q[ngrid] */,
                                   const double* /* qefield: Fortran (3,ngrid) */,
                                   const double* /* za[nsph] */,
                                   double* /* grad_rA: Fortran (3,nsph) */) moist_API_SUFFIX__V_0_5;

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
/// Call moist_get_cavity_sizes first and pass the capacities you allocated the
/// arrays with (see the capacity convention at the top of this header).
moist_API_ENTRY void moist_API_CALL
moist_get_drop_results(moist_error /* error */,
                      moist_cavity /* cavity */,
                      const int /* ngrid_cap: allocated grid capacity, must be
                                   >= the ngrid reported by
                                   moist_get_cavity_sizes */,
                      const int /* nsph_cap: allocated sphere capacity, must be
                                   >= the nsph reported by
                                   moist_get_cavity_sizes */,
                      double* /* area */,
                      double* /* volume */,
                      int* /* ngrid */,
                      int* /* nmax */,
                      int* /* nsph */,
                      double* /* xyz: Fortran (3,ngrid_cap) */,
                      double* /* normal: Fortran (3,ngrid_cap) */,
                      double* /* wleb[ngrid_cap] */,
                      double* /* a[ngrid_cap] */,
                      double* /* r_iI0[ngrid_cap] */,
                      double* /* f[ngrid_cap] */,
                      double* /* rho[ngrid_cap] */,
                      int* /* owner[ngrid_cap] - 0-based atom indices (0 to nsph-1) */,
                      bool* /* converged[ngrid_cap] */,
                      double* /* radii[nsph_cap] */,
                      double* /* asph[nsph_cap] */) moist_API_SUFFIX__V_0_6;

/// Delete DROP cavity handle (legacy - use moist_delete_cavity instead)
moist_API_ENTRY void moist_API_CALL
moist_delete_drop_cavity(moist_cavity* /* cavity */) moist_API_SUFFIX__V_0_5;
