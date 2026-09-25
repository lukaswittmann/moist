#pragma once
#ifndef moist_CFFI
#include <stdint.h>
#include <stddef.h>
#endif

#ifdef __cplusplus
#define moist_API_ENTRY extern "C"
#else
#define moist_API_ENTRY extern
#ifndef moist_CFFI
#include <stdbool.h>
#endif
#endif
#define moist_API_CALL

/// Current API contract. Feature-test with #ifdef moist_API_SUFFIX__V_1_0.
/// Older contract tags are no longer defined.
#define moist_API_SUFFIX__V_1_0

/*
 * ARRAY LAYOUT CONVENTION -- read this before allocating
 *
 * Every array crossing this API is a FLAT buffer in C row-major order.
 * The LAST dimension is contiguous. C dimensions reverse the native Fortran
 * dimensions without changing the flat buffer. Thus xyz[ngrid][3] corresponds
 * to Fortran xyz(3,ngrid), with component k of point i at xyz[3*i+k].
 * Derivative axes follow that same exact reversal, not a universal
 * value-then-derivative convention: see each tensor's named indices below.
 * A row-major (d1,d2,d3) buffer stores (i,j,k) at (i*d2+j)*d3+k.
 *
 * Annotations of the form `[n]` are plain 1-D arrays of n elements and mean
 * exactly what they say.
**/

/*
 * ARRAY SIZES -- who states how large a buffer is
 *
 * The entries of the host loop take no size: moist_answer_coupling_request,
 * moist_get_coupling_request_width, moist_get_response_array, the two name
 * getters moist_get_coupling_request_name and moist_get_response_item_name,
 * and the moist_get_cavity_field_* entries. moist knows the size already --
 * the grid size of the coupling, response or cavity and the extents of the
 * named output, array or field -- and reads or writes exactly that many
 * elements: a row-major (ngrid, dims...) array, or the `count`
 * moist_get_cavity_field_info reports. A pointer cannot say how large the
 * host's buffer is, so nothing checks it: the host allocates the documented
 * shape. Names are written with their terminator into MOIST_NAME_MAX + 1
 * (field names MOIST_FIELD_NAME_MAX + 1) characters.
 *
 * The other numerical array getters take the
 * capacities the caller allocated them with as `int` arguments
 * (`ngrid_cap`, `nsph_cap`, ...), passed BY VALUE. A capacity is checked
 * before a single element is written: if it is smaller than the value the
 * corresponding size query reports, the call sets an API
 * error and returns without touching any buffer.
 *
 * A capacity LARGER than required is accepted, so a host may allocate one
 * max-size buffer once and reuse it across geometry steps. In that case the
 * capacity defines the array's declared dimensions for
 * the layout convention above: `xyz: row-major (ngrid_cap, 3)` is indexed with
 * stride 3 for ngrid_cap points, and only the leading ngrid of them are
 * written; padding is untouched and must not be contracted. For higher-rank
 * arrays, the stride of each dimension is the product of dimensions to its
 * right, including capacities.
**/

/// Error handle class
typedef struct moist_error_s* moist_error;

/// Status returned by moist_check_error; details are read with moist_get_error.
typedef enum {
    moist_success = 0,
    moist_failure = 1,
    moist_invalid_error = 2
} moist_status;

/// Molecular structure data class
typedef struct moist_structure_s* moist_structure;

/// Solvation model class
typedef struct moist_model_s* moist_model;

/// Solvation-model component class
typedef struct moist_component_s* moist_component;

/// PCM linear-solver selection, as accepted by the CPCM and COSMO constructors
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
typedef struct moist_cavity_s* moist_cavity;

/// Radii model class
typedef struct moist_radii_s* moist_radii;

/// Host coupling class: the request list of one general model, minted by
/// moist_new_coupling and staged per phase by the moist_prepare_model_*
/// entry points (see HOST COUPLING PROTOCOL below).
typedef struct moist_coupling_s* moist_coupling;

/// Response class: the host part of one phase, filled by the
/// moist_get_model_* entry points and read back by walking its items.
typedef struct moist_response_s* moist_response;

/// Requests, response items, outputs and arrays are identified by scientific
/// NAME, never by a numeric tag or handle. Request names: "point_potential",
/// "gaussian_potential", "gaussian_moments"; their outputs: "phi", "dphi_dr",
/// "dphi_dxi", "gt", "pt", "mt", "rt". Response items: "potential_adjoint"
/// ("w_phi"), "density" ("w_rho", "w_grad_rho", "w_hess_rho"),
/// "gostshyp_amplitude" ("w_overlap", "w_normal_deriv"). Names are at most
/// MOIST_NAME_MAX characters.
#define MOIST_NAME_MAX 32


/* Handles and their borrowed aliases must not be used concurrently. Serialize
 * calls sharing a model, coupling, cavity, response, or error handle. Independent
 * object graphs may be used on separate threads; callbacks follow the rule below.
 */

/// Callback ABI for external isodensity DROP level set functions.
/// The callback receives a point in Bohr and must write the ELECTRON DENSITY
/// rho(r) there and its spatial gradient -- the bare density, not a level set:
/// do not subtract an isovalue and do not flip the sign.  moist forms the level
/// set itself as S = scale * (rho_iso - rho), with rho_iso and scale supplied to
/// moist_new_isodensity_callback_lsf(), which is what puts S in the DROP
/// sign convention (interior negative, exterior positive).
/// Calls may run concurrently on multiple threads with the same context pointer.
/// The callback must be reentrant; synchronize shared state or use thread-local scratch.
/// The `d2rho` and `d3rho` pointers are optional: they are NULL when
/// the cavity does not need that derivative order (e.g. the value+gradient-only
/// projection phase), and a NULL pointer means the callback should skip
/// computing that derivative, not merely skip writing it.  Dereferencing a NULL
/// `d2rho`/`d3rho` is a host bug, so branch on them before writing.
/// `d2rho` and `d3rho` follow the C row-major layout convention above; both are
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
/// See test_isodensity_callback_cavity() in test/api/example.c for a callback
/// written against this contract.
#ifdef __cplusplus
extern "C" {
#endif
typedef int (*moist_isodensity_lsf_callback)(void* context,
                                             const double* point /* [3] */,
                                             double* rho,
                                             double* drho /* [3] */,
                                             double* d2rho /* : row-major (3, 3), or NULL */,
                                             double* d3rho /* : row-major (3, 3, 3), or NULL */);


#ifdef __cplusplus
}
#endif

/// Level-set function, copied into a DROP cavity.
typedef struct moist_lsf_s* moist_lsf;

/* Options: pass NULL to a constructor for defaults, or initialize with
 * moist_init_*_options(error, &options, sizeof options) before overriding fields.
 * struct_size is maintained by the initializer. Fields may only be appended;
 * existing field offsets and meanings are stable within this API version.
 * New fields must start at or beyond the original sizeof, including its tail padding.
 * The full 1.0 layout is the minimum accepted size. Future libraries accept
 * older prefixes and default appended fields. Larger structs are accepted;
 * unknown trailing fields are ignored on input and untouched by initializers.
 * Initializers zero reserved fields and padding within the supported layout.
 */
typedef struct {
    size_t struct_size;
    int nleb;
    bool debug;
    int verbosity;
    bool do_fine;
    double tolerance;
    int proj_maxiter;
    int proj_level;
    double branch_weight_s;
    double rho_grid_h;
    int wleb_prune_level;
    int reserved0;  /* Initialized to zero, ignored on input. Do not reuse. */
    /* --- end of 1.0 layout --- */
} moist_drop_options;

moist_API_ENTRY void moist_API_CALL
moist_init_drop_options(moist_error error, moist_drop_options *options,
                        size_t struct_size) moist_API_SUFFIX__V_1_0;

typedef struct {
    size_t struct_size;
    int nleb;
    bool debug;
    int verbosity;
    double cut_a;
    double cut_f;
    /* --- end of 1.0 layout --- */
} moist_iswig_options;

moist_API_ENTRY void moist_API_CALL
moist_init_iswig_options(moist_error error, moist_iswig_options *options,
                        size_t struct_size) moist_API_SUFFIX__V_1_0;

typedef struct {
    size_t struct_size;
    double blend_k;
    double blend_1b;
    double blend_2b;
    double blend_3b;
    /* --- end of 1.0 layout --- */
} moist_svdw_options;

moist_API_ENTRY void moist_API_CALL
moist_init_svdw_options(moist_error error, moist_svdw_options *options,
                        size_t struct_size) moist_API_SUFFIX__V_1_0;

typedef struct {
    size_t struct_size;
    double a1;
    double a2;
    double c;
    int m;
    int reserved0;  /* Initialized to zero, ignored on input. Do not reuse. */
    /* --- end of 1.0 layout --- */
} moist_cfc_options;

moist_API_ENTRY void moist_API_CALL
moist_init_cfc_options(moist_error error, moist_cfc_options *options,
                        size_t struct_size) moist_API_SUFFIX__V_1_0;

typedef struct {
    size_t struct_size;
    double rho_iso;
    double scale;
    /* --- end of 1.0 layout --- */
} moist_isodensity_options;

moist_API_ENTRY void moist_API_CALL
moist_init_isodensity_options(moist_error error, moist_isodensity_options *options,
                        size_t struct_size) moist_API_SUFFIX__V_1_0;

typedef struct {
    size_t struct_size;
    bool debug;
    int verbosity;
    /* --- end of 1.0 layout --- */
} moist_model_options;

moist_API_ENTRY void moist_API_CALL
moist_init_model_options(moist_error error, moist_model_options *options,
                        size_t struct_size) moist_API_SUFFIX__V_1_0;

typedef struct {
    size_t struct_size;
    moist_pcm_solver solver;
    int solver_maxiter;  /* Iteration cap; only read by moist_pcm_solver_iterative. */
    double solver_tol;   /* Residual threshold; only read by moist_pcm_solver_iterative. */
    /* --- end of 1.0 layout --- */
} moist_pcm_options;

moist_API_ENTRY void moist_API_CALL
moist_init_pcm_options(moist_error error, moist_pcm_options *options,
                        size_t struct_size) moist_API_SUFFIX__V_1_0;

/// Banner text returned by moist_get_banner.
typedef enum {
    moist_banner_full = 0, moist_banner_short = 1,
    moist_banner_ascii = 2, moist_banner_build = 3
} moist_banner;

/// Copy banner text, or query its length with buffer=NULL and capacity=0.
/// length is required and receives the byte count excluding the NUL terminator.
/// With capacity>0, buffer is required and always NUL-terminated, even on error.
/// Capacities above SIZE_MAX/2 are rejected, but a supplied buffer is still cleared.
/// Truncation is successful: length >= capacity means a larger buffer is needed.
/// On error length is unchanged. The host chooses where to print.
moist_API_ENTRY void moist_API_CALL
moist_get_banner(moist_error error, moist_banner style, char *buffer,
                 size_t capacity, size_t *length) moist_API_SUFFIX__V_1_0;

/// Create an LSF. NULL options selects defaults.
moist_API_ENTRY moist_lsf moist_API_CALL
moist_new_svdw_lsf(moist_error error, const moist_svdw_options *options) moist_API_SUFFIX__V_1_0;
moist_API_ENTRY moist_lsf moist_API_CALL
moist_new_cfc_lsf(moist_error error, const moist_cfc_options *options) moist_API_SUFFIX__V_1_0;

/// Create a callback LSF. Callback and context must outlive every cavity/model copy.
moist_API_ENTRY moist_lsf moist_API_CALL
moist_new_isodensity_callback_lsf(moist_error error, moist_isodensity_lsf_callback callback,
                                  void *context, const moist_isodensity_options *options) moist_API_SUFFIX__V_1_0;

/// Create an LSF owning its Gaussian basis. shell_atom indices are zero-based.
/// exps and coeffs contain sum(shell_nprim) values; coefficients include primitive
/// normalization. Set the density on the constructed cavity, not the source LSF.
moist_API_ENTRY moist_lsf moist_API_CALL
moist_new_isodensity_lsf(moist_error error, int nshell, const int *shell_atom,
                        const int *shell_l, const int *shell_nprim,
                        const double *exps, const double *coeffs,
                        const moist_isodensity_options *options) moist_API_SUFFIX__V_1_0;
moist_API_ENTRY void moist_API_CALL
moist_delete_lsf(moist_lsf *lsf) moist_API_SUFFIX__V_1_0;

/// Construct DROP from a required LSF. LSF and radii are copied.
/// NULL radii selects CPCM radii; NULL options selects defaults.
moist_API_ENTRY moist_cavity moist_API_CALL
moist_new_drop_cavity(moist_error error, moist_lsf lsf, moist_radii radii,
                      const moist_drop_options *options) moist_API_SUFFIX__V_1_0;

/// Construct iSwiG, copying radii. NULL radii selects CPCM; NULL options defaults.
moist_API_ENTRY moist_cavity moist_API_CALL
moist_new_iswig_cavity(moist_error error, moist_radii radii,
                       const moist_iswig_options *options) moist_API_SUFFIX__V_1_0;

/// Construct a model owning a cavity copy. Add components before the first update.
moist_API_ENTRY moist_model moist_API_CALL
moist_new_model(moist_error error, moist_cavity cavity,
                const moist_model_options *options) moist_API_SUFFIX__V_1_0;

/// Construct PCM components; NULL options selects the default solver.
moist_API_ENTRY moist_component moist_API_CALL
moist_new_cpcm_component(moist_error error, double epsilon,
                         const moist_pcm_options *options) moist_API_SUFFIX__V_1_0;
moist_API_ENTRY moist_component moist_API_CALL
moist_new_cosmo_component(moist_error error, double epsilon,
                          const moist_pcm_options *options) moist_API_SUFFIX__V_1_0;

/*
 * Global API queries
**/

/// Obtain library version as major * 10000 + minor * 100 + patch
moist_API_ENTRY int moist_API_CALL
moist_get_version(void) moist_API_SUFFIX__V_1_0;

/// Copy the full release version, including prerelease suffix.
/// Buffer, length-query, truncation, and failure rules follow moist_get_banner.
moist_API_ENTRY void moist_API_CALL
moist_get_version_string(moist_error error, char *buffer, size_t capacity,
                         size_t *length) moist_API_SUFFIX__V_1_0;

/*
 * Error handle class
**/

/// Create new error handle object
moist_API_ENTRY moist_error moist_API_CALL
moist_new_error(void) moist_API_SUFFIX__V_1_0;

/// Return moist_success, moist_failure, or moist_invalid_error for a NULL handle.
/// Returns int so C and C++ hosts can use it directly as a status.
moist_API_ENTRY int moist_API_CALL
moist_check_error(moist_error error) moist_API_SUFFIX__V_1_0;

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
 *         int message_len = (int)sizeof(message);
 *         moist_get_error(error, message, &message_len);
 *         fprintf(stderr, "[moist Error] %s: %s\n", context, message);
 *         exit(EXIT_FAILURE);
 *     }
 *
 * A long-lived host should instead unwind, report, and keep running.
**/

/// Copy the diagnostic; no error produces an empty string, NULL error a diagnostic.
/// buffer and buffersize are required; *buffersize must be positive.
/// Writes at most *buffersize bytes, always NUL-terminating; long messages truncate.
moist_API_ENTRY void moist_API_CALL
moist_get_error(moist_error error,
                char* buffer,
                const int* buffersize) moist_API_SUFFIX__V_1_0;

/// Delete error handle object
moist_API_ENTRY void moist_API_CALL
moist_delete_error(moist_error* error) moist_API_SUFFIX__V_1_0;

/* ERROR CONVENTION
 *
 * Fallible operations report through a required, non-NULL moist_error handle.
 * Each fallible call replaces the previous diagnostic. Check moist_check_error
 * after each call, before using any result. Queries
 * write values through output parameters; values never encode error status.
 * Required pointers are checked; numeric outputs are unchanged on error.
 * Banner, version, and field-description getters clear writable buffers on error.
 * Constructors return a handle and report failure through the error handle;
 * release owned handles with their matching delete function.
 */

/*
 * Molecular structure data class
**/

/// Create new molecular structure data (quantities in Bohr)
moist_API_ENTRY moist_structure moist_API_CALL
moist_new_structure(moist_error error,
                    int natoms,
                    const int* numbers /* [natoms] */,
                    const double* positions /* : row-major (natoms, 3) */,
                    const double* lattice /* : row-major (3, 3), vectors in rows */,
                    const bool* periodic /* [3] */) moist_API_SUFFIX__V_1_0;

/// Delete molecular structure data
moist_API_ENTRY void moist_API_CALL
moist_delete_structure(moist_structure* mol) moist_API_SUFFIX__V_1_0;

/// Update coordinates and lattice parameters (quantities in Bohr)
moist_API_ENTRY void moist_API_CALL
moist_update_structure(moist_error error,
                       moist_structure mol,
                       const double* positions /* : row-major (natoms, 3) */,
                       const double* lattice /* : row-major (3, 3), vectors in rows */) moist_API_SUFFIX__V_1_0;

/*
 * Radii model class
**/

/// Create CPCM radii model
moist_API_ENTRY moist_radii moist_API_CALL
moist_new_cpcm_radii(moist_error error) moist_API_SUFFIX__V_1_0;

/// Create SMD radii model
moist_API_ENTRY moist_radii moist_API_CALL
moist_new_smd_radii(moist_error error) moist_API_SUFFIX__V_1_0;

/// Create D3 radii model
moist_API_ENTRY moist_radii moist_API_CALL
moist_new_d3_radii(moist_error error) moist_API_SUFFIX__V_1_0;

/// Create COSMO radii model
moist_API_ENTRY moist_radii moist_API_CALL
moist_new_cosmo_radii(moist_error error) moist_API_SUFFIX__V_1_0;

/// Create Bondi radii model
moist_API_ENTRY moist_radii moist_API_CALL
moist_new_bondi_radii(moist_error error) moist_API_SUFFIX__V_1_0;

/// Create custom radii model (must be populated before use)
moist_API_ENTRY moist_radii moist_API_CALL
moist_new_custom_radii(moist_error error) moist_API_SUFFIX__V_1_0;

/// Set custom radii from per-atom values (bohr)
moist_API_ENTRY void moist_API_CALL
moist_set_custom_radii_atoms(moist_error error,
                             moist_radii radii,
                             int natoms,
                             const double* atom_radii /* [natoms] */) moist_API_SUFFIX__V_1_0;

/// Set custom radii from per-element values (bohr)
moist_API_ENTRY void moist_API_CALL
moist_set_custom_radii_elements(moist_error error,
                                moist_radii radii,
                                int nentries,
                                const int* atomic_numbers /* [nentries] */,
                                const double* element_radii /* [nentries] */) moist_API_SUFFIX__V_1_0;

/// Delete radii model
moist_API_ENTRY void moist_API_CALL
moist_delete_radii(moist_radii* radii) moist_API_SUFFIX__V_1_0;

/*
 * Solvation model class
**/



/// Create a pressure-volume energy component equal to `pressure * cavity volume`.
moist_API_ENTRY moist_component moist_API_CALL
moist_new_pv_component(moist_error error,
                       double pressure) moist_API_SUFFIX__V_1_0;

/// Create a GOSTSHYP hydrostatic-pressure component at `pressure` in
/// Hartree/bohr^3. The component cannot form its own density traces: answer the
/// Gaussian-moment request of the coupling ("gt", "pt", "mt", "rt") with
/// moist_answer_coupling_request in every phase, and contract the
/// "gostshyp_amplitude" item of the response walk.
moist_API_ENTRY moist_component moist_API_CALL
moist_new_gostshyp_component(moist_error error,
                             double pressure) moist_API_SUFFIX__V_1_0;

/// Delete a standalone solvation-model component handle.
moist_API_ENTRY void moist_API_CALL
moist_delete_component(moist_component* component) moist_API_SUFFIX__V_1_0;


/// Append a copied component to a general model before its first update.
moist_API_ENTRY void moist_API_CALL
moist_add_model_component(moist_error error, moist_model model, moist_component component) moist_API_SUFFIX__V_1_0;

/* HOST COUPLING PROTOCOL
 *
 * The same six steps as the Fortran interface:
 *
 * 1. Update the model, then mint a coupling with moist_new_coupling.
 * 2. Stage a phase with moist_prepare_model_energy/response/gradient. Energy
 *    clears all answers; response and gradient keep valid answers. Every
 *    prepare starts a new walk over the requests.
 * 3. Walk the requests that still miss an output with
 *    `while (moist_next_coupling_request(err, cpl))`. Each such request is
 *    visited once per pass; false ends the pass and rewinds, so a second loop
 *    retries what is still missing. False also reports a failure: check err
 *    after the loop. Leaving the loop early resumes the pass at the next call.
 * 4. For the current request, read its name, ask which outputs are missing,
 *    compute them on the cavity grid -- read the grid from the model's cavity
 *    via moist_get_model_cavity + moist_get_cavity_field_real ("xyz", "xi0",
 *    ...) -- and submit each with moist_answer_coupling_request. A rejected
 *    output stays missing until a valid retry; the others survive.
 * 5. Call the moist_get_model_* entry matching the staged phase. It fails by
 *    name on a wrong staging or any missing output; evaluation does not
 *    consume answers.
 * 6. Walk the response with `while (moist_next_response_item(err, resp))`:
 *    every item present is visited once per pass. Read the current item's
 *    name, copy its arrays with moist_get_response_array and contract them;
 *    stop on an item the host cannot contract. An absent item is physics, not
 *    an error.
 *
 *     moist_prepare_model_energy(err, model, cpl);
 *     while (moist_next_coupling_request(err, cpl)) {
 *         char name[MOIST_NAME_MAX + 1];
 *         bool missing;
 *         moist_get_coupling_request_name(err, cpl, name);
 *         if (strcmp(name, "gaussian_potential") == 0) {
 *             moist_get_coupling_request_missing(err, cpl, "phi", &missing);
 *             if (missing) moist_answer_coupling_request(err, cpl, "phi", phi);
 *         }
 *     }
 *     moist_get_model_energy(err, model, cpl, &energy);
 *
 *     moist_prepare_model_response(err, model, cpl);
 *     ... walk the requests again as above ...
 *     moist_get_model_response(err, model, cpl, resp);
 *     while (moist_next_response_item(err, resp)) {
 *         char name[MOIST_NAME_MAX + 1];
 *         moist_get_response_item_name(err, resp, name);
 *         if (strcmp(name, "potential_adjoint") == 0) {
 *             moist_get_response_array(err, resp, "w_phi", w_phi);
 *             host_fock_potential(ngrid, w_phi, fock);
 *         } else {
 *             host_abort("unsupported response item", name);
 *         }
 *     }
 *
 * The "density" item is produced by the response phase only. Its weights are
 * dE/drho at fixed nuclei -- the same object in both phases -- so moist forms
 * them once and leaves the contraction out of the gradient phase, which a host
 * that does not need them would otherwise pay for. On a cavity whose surface
 * follows the density, run the response phase before the gradient and reuse
 * those weights: contracted with d rho/dP they complete the Fock matrix, and
 * with the basis-centre derivative d rho/dR they complete the nuclear gradient.
 * Going straight to the gradient phase there drops that term silently.
 *
 * A model update invalidates all its couplings' answers, even if the grid size
 * is unchanged; prepare again to refresh a coupling. The host must prepare
 * energy again for a new density, and also update a density-dependent cavity
 * before preparing. Release couplings before their model. Response arrays are
 * independent copies; get_response/get_gradient clear the response after input
 * validation, which also starts a new walk over its items.
 */
/// Mint a borrowed coupling after the model's first successful update.
/// Several couplings may coexist. Every component declares its requests.
/// The coupling survives later model updates; the prepare_* entry points
/// refresh it.
moist_API_ENTRY moist_coupling moist_API_CALL
moist_new_coupling(moist_error error,
                   moist_model model) moist_API_SUFFIX__V_1_0;

/// Delete a coupling handle and release its model-owned collection.
moist_API_ENTRY void moist_API_CALL
moist_delete_coupling(moist_coupling* coupling) moist_API_SUFFIX__V_1_0;

/// Create an empty response handle
moist_API_ENTRY moist_response moist_API_CALL
moist_new_response(moist_error error) moist_API_SUFFIX__V_1_0;

/// Delete a response handle
moist_API_ENTRY void moist_API_CALL
moist_delete_response(moist_response* response) moist_API_SUFFIX__V_1_0;

/// Stage the energy phase: mark every answer stale and start a new walk
moist_API_ENTRY void moist_API_CALL
moist_prepare_model_energy(moist_error error,
                           moist_model model,
                           moist_coupling coupling) moist_API_SUFFIX__V_1_0;

/// Stage the response phase, retaining valid answers; starts a new walk
moist_API_ENTRY void moist_API_CALL
moist_prepare_model_response(moist_error error,
                             moist_model model,
                             moist_coupling coupling) moist_API_SUFFIX__V_1_0;

/// Stage the gradient phase, retaining valid answers; starts a new walk
moist_API_ENTRY void moist_API_CALL
moist_prepare_model_gradient(moist_error error,
                             moist_model model,
                             moist_coupling coupling) moist_API_SUFFIX__V_1_0;

/// Solvation energy from a staged coupling. Fails by name on a coupling not
/// staged for the energy phase and on an unanswered required output, including
/// one whose answer was rejected. Adds the energy contribution to *energy.
/// Initialize it before the first call. On failure the accumulator is unchanged.
moist_API_ENTRY void moist_API_CALL
moist_get_model_energy(moist_error error,
                       moist_model model,
                       moist_coupling coupling,
                       double* energy) moist_API_SUFFIX__V_1_0;

/// Host part of the response phase from a staged coupling. `response` is
/// cleared after input validation and returns every item the model produces.
moist_API_ENTRY void moist_API_CALL
moist_get_model_response(moist_error error,
                         moist_model model,
                         moist_coupling coupling,
                         moist_response response) moist_API_SUFFIX__V_1_0;

/// Nuclear gradient from a staged coupling, plus the host part of the gradient
/// phase in `response` (cleared after input validation): the potential adjoint
/// and the GOSTSHYP amplitudes. Not the "density" item -- see HOST COUPLING
/// PROTOCOL above for the weights a density-backed cavity carries over from the
/// response phase. `gradient` is row-major (nat_cap, 3);
/// the capacity is checked against the atom count before anything is written.
/// Adds to gradient; initialize logical entries before the first call.
/// On failure the gradient accumulator is unchanged; padding is untouched.
moist_API_ENTRY void moist_API_CALL
moist_get_model_gradient(moist_error error,
                         moist_model model,
                         moist_coupling coupling,
                         moist_response response,
                         int nat_cap,
                         double* gradient /* : row-major (nat_cap, 3) */) moist_API_SUFFIX__V_1_0;

/// Advance to the next request with a missing output of the staged phase.
/// Returns false after the last one, rewinding so that the next call starts a
/// new pass, and on any failure, so it can drive a while loop: check error
/// afterwards. Leaving the loop early resumes the pass at the next call.
moist_API_ENTRY bool moist_API_CALL
moist_next_coupling_request(moist_error error, moist_coupling coupling) moist_API_SUFFIX__V_1_0;
/// NUL-terminated scientific name of the current request, written into
/// name[MOIST_NAME_MAX + 1]. An error when no request is current.
moist_API_ENTRY void moist_API_CALL
moist_get_coupling_request_name(moist_error error, moist_coupling coupling, char* name) moist_API_SUFFIX__V_1_0;
/// Whether the named output of the current request is required by the staged
/// phase and still missing, evaluated at the call rather than frozen when the
/// cursor arrived -- so answering an output and asking again reports it as no
/// longer missing. A name the request does not declare is not missing. An
/// error when no request is current.
moist_API_ENTRY void moist_API_CALL
moist_get_coupling_request_missing(moist_error error, moist_coupling coupling, const char* output,
    bool* missing) moist_API_SUFFIX__V_1_0;

/// Submit one named output of the current request. values is row-major
/// (ngrid, dims...) on the cavity grid, with the output's leading extents:
/// phi, dphi_dxi, gt: [ngrid]; dphi_dr, pt, rt: [ngrid][3]; mt: [ngrid][3][3].
/// moist reads exactly that many values. An error when no request is current,
/// for an unknown output, a NULL buffer or a non-finite value. A rejected
/// output stays missing until a valid retry; other outputs are unaffected.
moist_API_ENTRY void moist_API_CALL
moist_answer_coupling_request(moist_error error, moist_coupling coupling, const char* output,
    const double* values) moist_API_SUFFIX__V_1_0;
/// Copy the Gaussian moment exponents, width[ngrid] in bohr**-2, of the
/// current "gaussian_moments" request. They are the component's choice, not a
/// cavity field: read them, never recompute them.
moist_API_ENTRY void moist_API_CALL
moist_get_coupling_request_width(moist_error error, moist_coupling coupling, double* width) moist_API_SUFFIX__V_1_0;

/// Advance to the next item of the response. Every item present is visited
/// once per pass. Returns false after the last one, rewinding so that the next
/// call starts a new pass, at once for an empty response, and on any failure,
/// so it can drive a while loop: check error afterwards.
moist_API_ENTRY bool moist_API_CALL
moist_next_response_item(moist_error error, moist_response response) moist_API_SUFFIX__V_1_0;
/// NUL-terminated scientific name of the current item, written into
/// name[MOIST_NAME_MAX + 1]. An error when no item is current.
moist_API_ENTRY void moist_API_CALL
moist_get_response_item_name(moist_error error, moist_response response,
    char* name) moist_API_SUFFIX__V_1_0;
/// Copy one named array of the current item, row-major with the grid axis
/// first; moist writes exactly the shape listed here.
///   "potential_adjoint": "w_phi" [ngrid], w_phi_i = dE/dphi_i, contracted by
///     the host as `F_uv += sum_i w_phi_i V_uv(r_i)`; for a stationary PCM
///     this is the surface charge q_i.
///   "density": "w_rho" [ngrid], "w_grad_rho" [ngrid][3], "w_hess_rho"
///     [ngrid][3][3], the weights conjugate to the density value, gradient and
///     Hessian. Hessian weights use C indices [point][b][a] for native
///     weights(a,b,point); they need not be symmetric, so preserve this axis
///     order when contracting.
///   "gostshyp_amplitude": "w_overlap" [ngrid], "w_normal_deriv" [ngrid],
///     contracted as `F_uv += sum_i [w_overlap[i] g_uv,i + w_normal_deriv[i] f_uv,i]`.
/// An error when no item is current, for an array the current item does not
/// have, or a NULL buffer.
moist_API_ENTRY void moist_API_CALL
moist_get_response_array(moist_error error,
                         moist_response response,
                         const char* array,
                         double* values) moist_API_SUFFIX__V_1_0;

/// Update a solvation model with a molecular structure
moist_API_ENTRY void moist_API_CALL
moist_update_model(moist_error error,
                             moist_model model,
                             moist_structure mol) moist_API_SUFFIX__V_1_0;

/// Get a borrowed cavity handle from a solvation model.
/// The returned handle is valid as long as the parent model exists, but cannot
/// be rebuilt independently: moist_update_cavity
/// reject borrowed handles.
/// Use moist_delete_cavity to release the handle (does NOT destroy the model's cavity).
moist_API_ENTRY moist_cavity moist_API_CALL
moist_get_model_cavity(moist_error error,
                                 moist_model model) moist_API_SUFFIX__V_1_0;

/// Delete solvation model handle
moist_API_ENTRY void moist_API_CALL
moist_delete_model(moist_model* model) moist_API_SUFFIX__V_1_0;

/*
 * DROP cavity class
**/

/*
 * Type-specific constructors
**/







/// Query the internal isodensity cartesian-component layout (DROP-specific).
/// Two-pass: call with the array pointers NULL to read ncart/nshell, then
/// allocate and call again. `shell_cart_offset` is 0-based (length nshell+1);
/// `comp_lx/ly/lz` are the per-component monomial powers (length ncart) defining
/// the ordering the host must match when building its density transform. Any of
/// the output pointers may be NULL to skip it.
moist_API_ENTRY void moist_API_CALL
moist_get_isodensity_cart_layout(moist_error error,
                                 moist_cavity cavity,
                                 int* ncart,
                                 int* nshell,
                                 int* shell_cart_offset /* [nshell+1] */,
                                 int* comp_lx /* [ncart] */,
                                 int* comp_ly /* [ncart] */,
                                 int* comp_lz /* [ncart] */) moist_API_SUFFIX__V_1_0;

/// Install the cartesian-monomial density matrix for the internal isodensity LSF.
/// `dcart` is the ncart-by-ncart density matrix, row-major (ncart, ncart): element
/// (p,q) at `[p*ncart + q]`. Call before each moist_update_cavity.
moist_API_ENTRY void moist_API_CALL
moist_set_isodensity_density(moist_error error,
                             moist_cavity cavity,
                             int ncart,
                             const double* dcart /* : row-major (ncart, ncart) */) moist_API_SUFFIX__V_1_0;

/// Set density on a model owning an internal isodensity cavity. Copies dcart.
/// Invalidates model results and every coupling of the model; call moist_update_model next.
moist_API_ENTRY void moist_API_CALL
moist_set_model_isodensity_density(moist_error error, moist_model model, int ncart,
                                   const double *dcart /* [ncart][ncart] */) moist_API_SUFFIX__V_1_0;

/// Return the master numerical tolerance currently configured on a DROP cavity.
/// Reports an error for non-DROP cavity types.
moist_API_ENTRY void moist_API_CALL
moist_get_drop_cavity_tolerance(moist_error error,
                                moist_cavity cavity,
                                double* tolerance) moist_API_SUFFIX__V_1_0;

/*
 * Generic cavity operations (Tier 1 - work on all cavity types)
**/

/// Generic update cavity - works for all cavity types
moist_API_ENTRY void moist_API_CALL
moist_update_cavity(moist_error error,
                    moist_cavity cavity,
                    moist_structure mol) moist_API_SUFFIX__V_1_0;

/// Get generic cavity sizes - works for all cavity types
/// Returns ngrid (number of grid points) and nsph (number of spheres)
moist_API_ENTRY void moist_API_CALL
moist_get_cavity_sizes(moist_error error,
                       moist_cavity cavity,
                       int* ngrid,
                       int* nsph) moist_API_SUFFIX__V_1_0;

/// Get generic cavity results - works for all cavity types
/// Returns fields shared by all cavity types. `converged` contains DROP's
/// per-point projection status; cavity types without a projection solve return
/// true for every retained point.
/// Call moist_get_cavity_sizes first to get ngrid, nsph, then allocate arrays
/// and pass the capacities you allocated them with (see the capacity
/// convention at the top of this header). Nothing is written when either
/// capacity is smaller than the current size; an API error is set instead.
moist_API_ENTRY void moist_API_CALL
moist_get_cavity_results(moist_error error,
                         moist_cavity cavity,
                         int ngrid_cap /* : allocated grid capacity,
                                      at least the ngrid reported by
                                      moist_get_cavity_sizes */,
                         int nsph_cap /* : allocated sphere capacity,
                                      at least the nsph reported by
                                      moist_get_cavity_sizes */,
                         double* area,
                         double* volume,
                         int* ngrid,
                         int* nsph,
                         double* xyz /* : row-major (ngrid_cap, 3) */,
                         double* a /* [ngrid_cap] */,
                         int* owner /* [ngrid_cap] - 0-based atom indices (0 to nsph-1) */,
                         bool* converged /* [ngrid_cap] */,
                         double* radii /* [nsph_cap] */,
                         double* asph /* [nsph_cap] */) moist_API_SUFFIX__V_1_0;

/// Generic delete cavity - works for all cavity types
moist_API_ENTRY void moist_API_CALL
moist_delete_cavity(moist_cavity* cavity) moist_API_SUFFIX__V_1_0;

/*
 * Type-specific getters (Tier 2 - DROP-specific fields)
**/

/// Named result fields (Tier 2 - everything a cavity holds, by name)
///
/// A cavity declares the per-point and per-sphere arrays it currently holds,
/// each under the name it uses internally. Enumerate them with
/// moist_get_cavity_field_count + moist_get_cavity_field_info, then read one
/// with the accessor matching its type tag. Extending a cavity with a new
/// result therefore needs no new entry point here.
///
/// A field that was not computed is NOT declared: the optional DROP properties
/// (curvature, grid-point density, ...) appear only once the matching property
/// request was set, and asking for one that is absent is an error rather than a
/// buffer of zeros. Nothing is written to a rejected buffer.
///
/// Names are the cavity's own, so `xyz`, `a`, `owner`, `radii` and `asph` carry
/// exactly the values moist_get_cavity_results reports, `owner` 0-based
/// included. DROP adds the projection results -- `numbering`, `anchor_id`,
/// `branch`, `branch_count`, `wbranch`, `wleb`, `rho`, `r_iI0`, `normal0`,
/// `converged` and the diagnostics -- and iSwiG adds its own `numbering`.

/// Element type tags reported by moist_get_cavity_field_info. The values are
/// part of the contract; read a field with the accessor matching its tag.
#define MOIST_FIELD_REAL 1
#define MOIST_FIELD_INT 2
#define MOIST_FIELD_BOOL 3

/// Highest rank a field can have, i.e. the length of the `dims` buffer
/// moist_get_cavity_field_info writes.
#define MOIST_FIELD_MAX_RANK 2
/// Maximum field-name length excluding the NUL terminator. Mirrored by
/// `max_field_name_len` (api.f90).
#define MOIST_FIELD_NAME_MAX 64

/// Number of named result fields the cavity currently holds.
/// The count changes when the cavity is rebuilt or its property request
/// changes, so read it again after moist_update_cavity.
moist_API_ENTRY void moist_API_CALL
moist_get_cavity_field_count(moist_error error,
                             moist_cavity cavity,
                             int* nfield) moist_API_SUFFIX__V_1_0;

/// Describe one readable field by position.
/// `index` is 0-based and must be below the reported count. `dtype` is one of
/// the MOIST_FIELD_* tags, `rank` is 0 for a scalar, `dims` receives
/// MOIST_FIELD_MAX_RANK extents with the slowest-varying one first (unused
/// entries are 1), and `count` is the number of elements a read writes --
/// allocate that many for it. The name is written with its terminator into
/// name[MOIST_FIELD_NAME_MAX + 1].
moist_API_ENTRY void moist_API_CALL
moist_get_cavity_field_info(moist_error error,
                            moist_cavity cavity,
                            int index /* : 0-based, below the count from
                                         moist_get_cavity_field_count */,
                            char* name /* [MOIST_FIELD_NAME_MAX + 1] */,
                            int* dtype /* : one of MOIST_FIELD_* */,
                            int* rank /* : 0 for a scalar */,
                            int* dims /* [MOIST_FIELD_MAX_RANK] */,
                            int* count /* : elements a read writes */) moist_API_SUFFIX__V_1_0;

/// Copy a field description or query its length, following moist_get_banner.
/// Pass about=NULL and capacity=0 to query length (excluding the terminator).
moist_API_ENTRY void moist_API_CALL
moist_get_cavity_field_about(moist_error error, moist_cavity cavity,
                             const char* name, char* about,
                             size_t capacity, size_t* length) moist_API_SUFFIX__V_1_0;

/// Read a MOIST_FIELD_REAL field by name: the `count` elements
/// moist_get_cavity_field_info reports; a rank-2 field is written flat in
/// row-major order per the layout convention at the top of this header.
moist_API_ENTRY void moist_API_CALL
moist_get_cavity_field_real(moist_error error,
                            moist_cavity cavity,
                            const char* name,
                            double* values /* [count] */) moist_API_SUFFIX__V_1_0;

/// Read a MOIST_FIELD_INT field by name.
/// Sphere indices are handed out 0-based, matching moist_get_cavity_results;
/// the DROP point ids `numbering` and `anchor_id` and the 1-based `branch` are
/// passed through as the cavity stores them.
moist_API_ENTRY void moist_API_CALL
moist_get_cavity_field_int(moist_error error,
                           moist_cavity cavity,
                           const char* name,
                           int* values /* [count] */) moist_API_SUFFIX__V_1_0;

/// Read a MOIST_FIELD_BOOL field by name.
moist_API_ENTRY void moist_API_CALL
moist_get_cavity_field_bool(moist_error error,
                            moist_cavity cavity,
                            const char* name,
                            bool* values /* [count] */) moist_API_SUFFIX__V_1_0;

/// Assemble A-matrix and compute xi values
/// Must be called before accessing xi or using the A-matrix
/// Call moist_get_cavity_sizes first to get ngrid for array allocation, then
/// pass the capacity you allocated with. Nothing is written when it is smaller than the current size.
/// Works for every Gaussian-discretized cavity (no longer DROP-specific)
moist_API_ENTRY void moist_API_CALL
moist_assemble_amat(moist_error error,
                    moist_cavity cavity,
                    int ngrid_cap /* : allocated grid capacity, at least
                                 the ngrid reported by moist_get_cavity_sizes */,
                    double* amat0 /* : row-major (ngrid_cap, ngrid_cap), symmetric,
                               so the transpose reading gives the same buffer */,
                    double* xi /* [ngrid_cap] */) moist_API_SUFFIX__V_1_0;

/// Return Gaussian PCM widths and switching factors.
/// Call moist_get_cavity_sizes first and pass the capacity allocated for both
/// arrays. Nothing is written when it is smaller than the current size.
moist_API_ENTRY void moist_API_CALL
moist_get_cavity_gaussian(moist_error error,
                          moist_cavity cavity,
                          int ngrid_cap /* : allocated grid capacity, at
                                       least the ngrid reported by
                                       moist_get_cavity_sizes */,
                          double* xi /* [ngrid_cap] */,
                          double* f /* [ngrid_cap] */) moist_API_SUFFIX__V_1_0;

/*
 * Gradient API (Tier 3 - Cavity and A-matrix gradients)
**/

/// Compute cavity gradient w.r.t. nuclear coordinates
/// Must be called after moist_update_cavity and before moist_get_cavity_gradient
moist_API_ENTRY void moist_API_CALL
moist_compute_cavity_gradient(moist_error error,
                              moist_cavity cavity) moist_API_SUFFIX__V_1_0;

/// Compute anchor-only nuclear derivatives (callback/isodensity LSF, DROP-specific)
/// Must be called after moist_update_cavity and before the *_rA contractions.
/// Restricts each grid point's nuclear coupling to its owner atom's rigid anchor
/// motion (the level set field's nuclear derivatives are zero for callback LSFs).
moist_API_ENTRY void moist_API_CALL
moist_compute_anchor_gradient(moist_error error,
                              moist_cavity cavity) moist_API_SUFFIX__V_1_0;

/// Get anchor-channel nuclear derivatives from the anchor pass
/// Must call moist_compute_anchor_gradient (or moist_compute_cavity_gradient) first
/// Call moist_get_cavity_sizes first and pass the capacities you allocated the
/// arrays with (see the capacity convention at the top of this header).
/// Arrays (row-major order; see the layout convention at the top of this header):
///   xyz1_rA  row-major (ngrid_cap, nsph_cap, 3, 3) - d(r_i)_j / d(R_A)_alpha  (i, A, alpha, j)
///                                       element at [j + 3*(alpha + 3*(A + nsph_cap*i))]
///   xi1_rA   row-major (ngrid_cap, nsph_cap, 3)   - d(xi_i)  / d(R_A)_alpha  (i, A, alpha)
///                                       element at [alpha + 3*(A + nsph_cap*i)]
///   a_i1_rA  row-major (ngrid_cap, nsph_cap, 3)   - d(a_i)   / d(R_A)_alpha  (i, A, alpha)
///   v_i1_rA  row-major (ngrid_cap, nsph_cap, 3)   - d(v_i)   / d(R_A)_alpha  (i, A, alpha)
///   A_tot1_rA row-major (nsph_cap, 3)    - d(total area)   / d(R_A)_alpha  (A, alpha)
///   V_tot1_rA row-major (nsph_cap, 3)    - d(total volume) / d(R_A)_alpha  (A, alpha)
/// a_i1_rA/v_i1_rA are the un-summed per-point counterparts of A_tot1_rA/V_tot1_rA.
/// The grid point area carries a switching-function dependence (a_i ~ f_i/xi_i^2),
/// so a_i1_rA is not recoverable from xi1_rA alone; it is the area route a generic
/// geometric surface functional (e.g. GOSTSHYP) needs.
moist_API_ENTRY void moist_API_CALL
moist_get_anchor_gradient(moist_error error,
                          moist_cavity cavity,
                          int nsph_cap /* : allocated sphere capacity, at
                                       least the nsph reported by
                                       moist_get_cavity_sizes */,
                          int ngrid_cap /* : allocated grid capacity, at
                                       least the ngrid reported by
                                       moist_get_cavity_sizes */,
                          double* xyz1_rA /* : row-major (ngrid_cap, nsph_cap, 3, 3) */,
                          double* xi1_rA /* : row-major (ngrid_cap, nsph_cap, 3) */,
                          double* a_i1_rA /* : row-major (ngrid_cap, nsph_cap, 3) */,
                          double* v_i1_rA /* : row-major (ngrid_cap, nsph_cap, 3) */,
                          double* A_tot1_rA /* : row-major (nsph_cap, 3) */,
                          double* V_tot1_rA /* : row-major (nsph_cap, 3) */) moist_API_SUFFIX__V_1_0;

/// Get cavity gradient arrays (DROP-specific)
/// Must call moist_compute_cavity_gradient first
/// Call moist_get_cavity_sizes first to get ngrid, nsph for array allocation,
/// then pass the capacities you allocated with (see the capacity convention at
/// the top of this header).
/// Arrays (row-major order; see the layout convention at the top of this header):
///   A_tot1_rA row-major (nsph_cap, 3)    - gradient of total area w.r.t. nuclear coords
///   V_tot1_rA row-major (nsph_cap, 3)    - gradient of total volume w.r.t. nuclear coords
///   asph1_rA  row-major (nsph_cap, nsph_cap, 3) - gradient of per-sphere areas (perturbed atom, owner, xyz)
///   vsph1_rA  row-major (nsph_cap, nsph_cap, 3) - gradient of per-sphere volumes (perturbed atom, owner, xyz)
///   xyz1_rA   row-major (ngrid_cap, nsph_cap, 3, 3) - grid point position derivatives (grid, atom, perturbed_xyz, xyz)
///   r_iI1_rA  row-major (ngrid_cap, nsph_cap, 3) - gradient of grid-owner distances
///   rho1_rA   row-major (ngrid_cap, nsph_cap, 3) - gradient of rho (anchor-to-surface distance)
moist_API_ENTRY void moist_API_CALL
moist_get_cavity_gradient(moist_error error,
                          moist_cavity cavity,
                          int nsph_cap /* : allocated sphere capacity, at
                                       least the nsph reported by
                                       moist_get_cavity_sizes */,
                          int ngrid_cap /* : allocated grid capacity, at
                                       least the ngrid reported by
                                       moist_get_cavity_sizes */,
                          double* A_tot1_rA /* : row-major (nsph_cap, 3) */,
                          double* V_tot1_rA /* : row-major (nsph_cap, 3) */,
                          double* asph1_rA /* : row-major (nsph_cap, nsph_cap, 3) */,
                          double* vsph1_rA /* : row-major (nsph_cap, nsph_cap, 3) */,
                          double* xyz1_rA /* : row-major (ngrid_cap, nsph_cap, 3, 3) */,
                          double* r_iI1_rA /* : row-major (ngrid_cap, nsph_cap, 3) */,
                          double* rho1_rA /* : row-major (ngrid_cap, nsph_cap, 3) */) moist_API_SUFFIX__V_1_0;

/// Assemble a Gaussian PCM A-matrix with its nuclear derivatives.
/// Must call moist_compute_cavity_gradient first
/// Call moist_get_cavity_sizes first to get ngrid, nsph for array allocation,
/// then pass the capacities you allocated with (see the capacity convention at
/// the top of this header).
/// Arrays (row-major order; see the layout convention at the top of this header):
///   Amat0     row-major (ngrid_cap, ngrid_cap)  - CPCM A-matrix (symmetric, so the
///                                            transpose reading is identical)
///   Amat1_rA  row-major (ngrid_cap, ngrid_cap, nsph_cap, 3) - gradient of A-matrix
///                                            w.r.t. nuclear coords; element
///                                            (j,i,A,alpha) at
///                                            [alpha + 3*(A + nsph_cap*(i + ngrid_cap*j))]
///   xi        [ngrid_cap]                    - xi values (screening factors)
moist_API_ENTRY void moist_API_CALL
moist_get_amat_gradient(moist_error error,
                        moist_cavity cavity,
                        int nsph_cap /* : allocated sphere capacity, at
                                     least the nsph reported by
                                     moist_get_cavity_sizes */,
                        int ngrid_cap /* : allocated grid capacity, at
                                     least the ngrid reported by
                                     moist_get_cavity_sizes */,
                        double* Amat0 /* : row-major (ngrid_cap, ngrid_cap) */,
                        double* Amat1_rA /* : row-major (ngrid_cap, ngrid_cap, nsph_cap, 3) */,
                        double* xi /* [ngrid_cap] */) moist_API_SUFFIX__V_1_0;

/// Contract Gaussian PCM A-matrix derivatives with two grid vectors.
/// Computes grad_rA = sum_ij q1_i * (dA_ij/dR_A) * q2_j
/// Must call moist_compute_cavity_gradient first
/// Uses cavity-internal ngrid and nsph for array extents
/// Arrays (row-major order; see the layout convention at the top of this header):
///   q1[ngrid], q2[ngrid]               - contraction vectors
///   grad_rA   row-major (nsph, 3)         - contracted nuclear gradient contribution
moist_API_ENTRY void moist_API_CALL
moist_contract_amat1_q1q2_rA(moist_error error,
                             moist_cavity cavity,
                             const double* q1 /* [ngrid] */,
                             const double* q2 /* [ngrid] */,
                             double* grad_rA /* : row-major (nsph, 3) */) moist_API_SUFFIX__V_1_0;

/// Contract Gaussian PCM A-matrix derivatives to per-grid surface weights.
/// Computes weights satisfying:
///   q1^T dA q2 = sum_i w_xi[i] dxi[i] + w_f[i] df[i] + w_xyz[i,:].dxyz[i,:]
/// Uses cavity-internal ngrid for array extents.
/// Arrays (row-major order; see the layout convention at the top of this header):
///   q1[ngrid], q2[ngrid]               - contraction vectors
///   w_xi[ngrid]                        - contracted xi derivative weights
///   w_f[ngrid]                         - contracted switch-function weights
///   w_xyz     row-major (ngrid, 3)        - contracted coordinate derivative weights
moist_API_ENTRY void moist_API_CALL
moist_contract_amat1_q1q2_surface_weights(moist_error error,
                                          moist_cavity cavity,
                                          const double* q1 /* [ngrid] */,
                                          const double* q2 /* [ngrid] */,
                                          double* w_xi /* [ngrid] */,
                                          double* w_f /* [ngrid] */,
                                          double* w_xyz /* : row-major (ngrid, 3) */) moist_API_SUFFIX__V_1_0;

/// Contract DROP surface weights to per-grid LSF adjoint weights (DROP-specific)
/// Includes the projected-coordinate response from w_xyz and the electronic
/// xi path through wleb, cpjac, and the critical-gradient/focus switches.
/// The exported w_f is the anchor-only iSwiG overlap and has no electronic
/// LSF response for fixed nuclei.
/// Uses cavity-internal ngrid for array extents.
/// Arrays (row-major order; see the layout convention at the top of this header):
///   w_xi[ngrid]                        - xi derivative weights
///   w_f[ngrid]                         - anchor switch derivative weights
///   w_xyz     row-major (ngrid, 3)        - coordinate derivative weights
///   w_lsf0[ngrid]                      - LSF value adjoint weights
///   w_lsf1    row-major (ngrid, 3)        - LSF gradient adjoint weights
///   w_lsf2    row-major (ngrid, 3, 3)      - LSF Hessian adjoint weights
/// w_lsf2[i][b][a] is conjugate to the native Hessian(a,b,i).
/// Note: w_lsf2 is NOT symmetric
/// This entry point covers the xi/f/xyz channels only; it is a thin wrapper
/// that forwards NULL for the extended channels below.
moist_API_ENTRY void moist_API_CALL
moist_contract_surface_lsf_weights(moist_error error,
                                   moist_cavity cavity,
                                   const double* w_xi /* [ngrid] */,
                                   const double* w_f /* [ngrid] */,
                                   const double* w_xyz /* : row-major (ngrid, 3) */,
                                   double* w_lsf0 /* [ngrid] */,
                                   double* w_lsf1 /* : row-major (ngrid, 3) */,
                                   double* w_lsf2 /* : row-major (ngrid, 3, 3) */) moist_API_SUFFIX__V_1_0;

/// Extended DROP surface-to-LSF contraction with optional normal/curvature channels.
/// Same contract as moist_contract_surface_lsf_weights plus:
///   w_n       row-major (ngrid, 3)        - outward-normal derivative weights (or NULL)
///   w_k1[ngrid]                        - first principal curvature weights (or NULL)
///   w_k2[ngrid]                        - second principal curvature weights (or NULL)
/// Pass NULL for any optional channel that should be skipped.
moist_API_ENTRY void moist_API_CALL
moist_contract_surface_lsf_weights_extended(
                                   moist_error error,
                                   moist_cavity cavity,
                                   const double* w_xi /* [ngrid] */,
                                   const double* w_f /* [ngrid] */,
                                   const double* w_xyz /* : row-major (ngrid, 3) */,
                                   double* w_lsf0 /* [ngrid] */,
                                   double* w_lsf1 /* : row-major (ngrid, 3) */,
                                   double* w_lsf2 /* : row-major (ngrid, 3, 3) */,
                                   const double* w_n /* : row-major (ngrid, 3), or NULL */,
                                   const double* w_k1 /* [ngrid] or NULL */,
                                   const double* w_k2 /* [ngrid] or NULL */) moist_API_SUFFIX__V_1_0;

/// Contract direct nuclear and electronic PCM terms.
/// Must call moist_compute_cavity_gradient first
/// Uses cavity-internal ngrid and nsph for array extents
/// Arrays (row-major order; see the layout convention at the top of this header):
///   w_phi[ngrid]                       - potential adjoint dE/dphi_i (surface charge q_i of a stationary PCM)
///   w_xyz     row-major (ngrid, 3)        - total dE_host/dr_i (nuclear + electronic)
///   za[nsph]                           - nuclear charges Z_A
///   grad_rA   row-major (nsph, 3)         - contracted nuclear gradient contribution
moist_API_ENTRY void moist_API_CALL
moist_contract_pcm_nuclear_gradient(moist_error error,
                                   moist_cavity cavity,
                                   const double* w_phi /* [ngrid] */,
                                   const double* w_xyz /* : row-major (ngrid, 3) */,
                                   const double* za /* [nsph] */,
                                   double* grad_rA /* : row-major (nsph, 3) */) moist_API_SUFFIX__V_1_0;


/* Ownership: each new_* handle has one owner. Constructors copy configuration.
 * Delete couplings before their model. Borrowed model cavities must not outlive
 * that model; deleting their wrapper does not delete the model's cavity.
 * Responses own their results. Deletion sets the supplied handle to NULL.
 */
#ifndef moist_CFFI
#ifdef __cplusplus
inline void moist_delete(moist_error& handle) { moist_delete_error(&handle); }
inline void moist_delete(moist_structure& handle) { moist_delete_structure(&handle); }
inline void moist_delete(moist_model& handle) { moist_delete_model(&handle); }
inline void moist_delete(moist_component& handle) { moist_delete_component(&handle); }
inline void moist_delete(moist_cavity& handle) { moist_delete_cavity(&handle); }
inline void moist_delete(moist_radii& handle) { moist_delete_radii(&handle); }
inline void moist_delete(moist_lsf& handle) { moist_delete_lsf(&handle); }
inline void moist_delete(moist_coupling& handle) { moist_delete_coupling(&handle); }
inline void moist_delete(moist_response& handle) { moist_delete_response(&handle); }
#else
#define moist_delete(handle) _Generic((handle), \
    moist_error: moist_delete_error, \
    moist_structure: moist_delete_structure, \
    moist_model: moist_delete_model, \
    moist_component: moist_delete_component, \
    moist_cavity: moist_delete_cavity, \
    moist_radii: moist_delete_radii, \
    moist_lsf: moist_delete_lsf, \
    moist_coupling: moist_delete_coupling, \
    moist_response: moist_delete_response)(&(handle))
#endif
#endif
