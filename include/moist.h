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

/// Current API contract, feature-test with #ifdef moist_API_SUFFIX__V_1_0
#define moist_API_SUFFIX__V_1_0

/*
 * ARRAY LAYOUT CONVENTION -- read this before allocating
 *
 * - Every array crossing this API is a flat buffer in C row-major order
 * - Last dimension contiguous; C dimensions reverse the native Fortran
 *   dimensions without changing the flat buffer
 * - xyz[ngrid][3] <-> Fortran xyz(3,ngrid); component k of point i at xyz[3*i+k]
 * - Derivative axes follow the same reversal, not a universal
 *   value-then-derivative convention -- see each tensor's named indices below
 * - Row-major (d1,d2,d3) buffer: (i,j,k) at (i*d2+j)*d3+k
 * - `[n]` annotations are plain 1-D arrays of n elements
**/

/*
 * ARRAY SIZES -- who states how large a buffer is
 *
 * - No-size entries: moist_answer_coupling_request, the two name getters
 *   moist_get_coupling_request_name/moist_get_response_item_name, and the
 *   moist_get_cavity_field_*, moist_get_model_field_*,
 *   moist_get_model_component_field_*,
 *   moist_get_coupling_request_field_* and moist_get_response_field_* entries
 * - Each reads or writes exactly the grid size of the coupling or cavity and
 *   the extents of the named output or field: a row-major (ngrid, dims...)
 *   answer, or the `count` the matching *_field_info reports
 * - Not checked against the buffer size; host allocates the documented shape
 * - Names written with terminator into MOIST_NAME_MAX + 1 characters
 *   (field names: MOIST_FIELD_NAME_MAX + 1)
 *
 * - Other numerical array getters take caller-allocated capacities as `int`
 *   arguments (`ngrid_cap`, `nsph_cap`, ...), passed by value
 * - Capacity checked before any element is written: smaller than the size
 *   query's value sets an API error and returns without touching any buffer
 *
 * - Capacity larger than required is accepted; it then defines the array's
 *   declared dimensions for the layout convention above: `xyz: row-major
 *   (ngrid_cap, 3)` indexed with stride 3 for ngrid_cap points, only the
 *   leading ngrid written, padding untouched and must not be contracted
 * - Higher-rank arrays: stride of each dimension is the product of
 *   dimensions to its right, including capacities
**/

/// Error handle class
typedef struct moist_error_s* moist_error;

/// Shared run context for cavities and models
typedef struct moist_context_s* moist_context;

/// Status returned by moist_check_error; details read with moist_get_error
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

/// PCM linear-solver selection, shared by all PCM component constructors
/// Mirrors the Fortran `solver_type` enumerator in
/// src/moist/model/continuum/component/pcm/type.f90; numeric values are ABI
typedef enum {
    /// Dense factorization
    moist_pcm_solver_dense = 1,
    /// Restarted GMRES with block-Jacobi preconditioning
    moist_pcm_solver_iterative = 2
} moist_pcm_solver;

/// Density representation of the internal isodensity LSF
/// Mirrors the Fortran `moist_iso_basis` enumerator in
/// src/moist/cavity/drop/lsf/isodensity/gto.f90; numeric values are ABI
typedef enum {
    /// Cartesian monomials
    moist_isodensity_basis_cartesian = 1,
    /// Harmonic polynomials, 2l+1 per shell
    moist_isodensity_basis_spherical = 2
} moist_isodensity_basis;

/// Matrix-storage policy for a PCM solver
typedef enum {
    moist_pcm_storage_dense = 1,
    moist_pcm_storage_matrix_free = 2
} moist_pcm_storage;

/// DROP cavity class
typedef struct moist_cavity_s* moist_cavity;

/// Radii model class
typedef struct moist_radii_s* moist_radii;

/// Host coupling class: the request list of one solvation model or standalone
/// cavity, minted by moist_new_coupling or moist_new_cavity_coupling and staged
/// per phase by the moist_prepare_* entry points (see HOST COUPLING PROTOCOL
/// below)
typedef struct moist_coupling_s* moist_coupling;

/// Response class: the host part of one phase, filled by the
/// moist_get_model_* entry points and read back by walking its items
typedef struct moist_response_s* moist_response;

/// Requests, response items, outputs and arrays are identified by canonical
/// NAME, never by a numeric tag or handle. Request names: "point_potential",
/// "gaussian_potential", "gaussian_moments", "atomic_multipoles",
/// "atomic_charges", "radial_potential", "nuclear_charges", "density_matrix";
/// their outputs: "phi", "dphi_dr", "dphi_dxi", "d2phi_dr2", "d2phi_dr_dxi",
/// "gt", "pt", "mt", "rt", "q",
/// "mu", "theta", "z", "dcart"; request inputs:
/// "width" ("gaussian_moments"), "nbasis" ("density_matrix"). Response items: "potential_adjoint"
/// ("w_phi", "w_dphi_dr" when computed), "density" ("w_rho", "w_grad_rho", "w_hess_rho"),
/// "gaussian_amplitude" ("w_overlap", "w_normal_deriv"),
/// "atomic_multipole_adjoint" ("dg_dq", "dg_dmu", "dg_dtheta"),
/// "atomic_charge_adjoint" ("dg_dq"), "radial_potential_adjoint" ("dg_dphi"),
/// "parameter_adjoint" (scalars "dg_deps", ...; response phase only, on request).
/// Names at most MOIST_NAME_MAX characters
#define MOIST_NAME_MAX 32


/* Handles and their borrowed aliases: not for concurrent use
 * - Serialize calls sharing a context, model, coupling, cavity, response, or error handle
 * - Independent object graphs usable on separate threads; callbacks follow the rule below
 */

/// Callback ABI for external isodensity DROP level set functions
/// Point in Bohr; write the ELECTRON DENSITY rho(r) there and its spatial
/// gradient -- the bare density, not a level set: do not subtract an
/// isovalue, do not flip the sign. moist forms the level set itself as
/// S = scale * (rho_iso - rho), with rho_iso and scale supplied to
/// moist_new_isodensity_callback_lsf(): interior negative, exterior positive
/// in the DROP sign convention
/// Calls may run concurrently on multiple threads with the same context
/// pointer. Callback must be reentrant; synchronize shared state or use
/// thread-local scratch
/// `d2rho`/`d3rho` optional: NULL when the cavity does not need that
/// derivative order (e.g. the value+gradient-only projection phase); NULL
/// means skip computing that derivative, not merely skip writing it
/// Dereferencing a NULL `d2rho`/`d3rho` is a host bug -- branch on them
/// before writing
/// `d2rho`/`d3rho` follow the C row-major layout convention above; both
/// symmetric under any index permutation, storage order immaterial
///
/// RETURN VALUE -- the failure channel
/// - 0 after writing every requested buffer
/// - Nonzero reports evaluation failure: output buffers ignored, moist stops
///   calling the callback for the rest of the cavity build, unwinds its
///   parallel evaluation loops, and fails the enclosing moist_update_cavity(),
///   moist_build_cavity(), moist_update_model() or moist_build_model()
///   with a moist_error naming the returned status
/// - Only abort path from inside a callback; without it moist builds a
///   cavity from whatever is in the buffers
/// - Status value opaque to moist, echoed verbatim in the error message; may
///   carry host-defined error codes
///
/// Example: test_isodensity_callback_cavity() in test/api/example.c
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

/// Level-set function, copied into a DROP cavity
typedef struct moist_lsf_s* moist_lsf;

/// Construct a shared run context
/// nthreads is fixed for the context's lifetime: > 0 as given, 0 takes the
/// calling thread's omp_get_max_threads at construction, negative is an error.
/// The count sizes MOIST's own OpenMP regions; the host's OpenMP runtime and
/// the BLAS/LAPACK threading are left alone. verbosity (0 silent ... 6 debug,
/// 2 is the library default; levels in the API reference), debug and
/// warnings_to_stderr are the only logging settings; debug raises the
/// verbosity to 6. Warnings print to standard output, warnings_to_stderr also
/// copies them to standard error. Every model constructor requires a context,
/// a cavity or component constructor takes one optionally; each retains it
moist_API_ENTRY moist_context moist_API_CALL
moist_new_context(moist_error error, int nthreads, int verbosity,
                  bool debug, bool warnings_to_stderr) moist_API_SUFFIX__V_1_0;
moist_API_ENTRY void moist_API_CALL
moist_delete_context(moist_context *context) moist_API_SUFFIX__V_1_0;
moist_API_ENTRY void moist_API_CALL
moist_get_context_num_threads(moist_error error, moist_context context,
                             int *nthreads) moist_API_SUFFIX__V_1_0;
/// Print the timings recorded on a context to its output: nothing at
/// verbosity 0, the top-level phases at 1, the whole tree from 4
moist_API_ENTRY void moist_API_CALL
moist_print_context_timer(moist_error error,
                          moist_context context) moist_API_SUFFIX__V_1_0;
/// Print the timings recorded on a context as a compact list of seconds per
/// step, at most two levels deep, at the context's verbosity and to its output
moist_API_ENTRY void moist_API_CALL
moist_print_context_timer_simple(moist_error error,
                                 moist_context context) moist_API_SUFFIX__V_1_0;
/// Host timers: time the host's own work in the context's timer report
/// - Every host timer lands under a top-level "Host" node, apart from moist's
///   own timers: neither ever nests under the other. "Host" runs while any
///   host timer is open and carries the "Host" category of the report
/// - Host timers nest among themselves; a stop names the innermost open one
///   and anything else is an error that stops nothing
/// - Time only host work: moist calls made with a host timer open are counted
///   both there and in moist's own timers
/// - A name has 1 to MOIST_NAME_MAX bytes of UTF-8 (an ASCII character takes
///   one, any other two to four), no control character and no leading or
///   trailing space; a repeated name under the same parent accumulates into
///   one node
/// - Stop on failure paths too, into a separate moist_error: every call
///   replaces the diagnostic it is given, so reusing the failed one loses it
/// - Serialized with every other call on the context; never from a callback
moist_API_ENTRY void moist_API_CALL
moist_start_context_timer(moist_error error, moist_context context,
                          const char* name) moist_API_SUFFIX__V_1_0;
/// Stop the innermost open host timer, which must carry `name`
moist_API_ENTRY void moist_API_CALL
moist_stop_context_timer(moist_error error, moist_context context,
                         const char* name) moist_API_SUFFIX__V_1_0;

/* Options: pass NULL to a constructor for defaults, or initialize with
 * moist_init_*_options(error, &options, sizeof options) before overriding fields
 * - struct_size maintained by the initializer
 * - Fields may only be appended; existing field offsets and meanings stable
 *   within this API version
 * - New fields start at or beyond the original sizeof, including tail padding
 * - Full 1.0 layout is the minimum accepted size
 * - Future libraries accept older prefixes and default appended fields
 * - Larger structs accepted; unknown trailing fields ignored on input and
 *   untouched by initializers
 * - Initializers zero padding within the supported layout
 */
typedef struct {
    size_t struct_size;
    double spacing; /* Maximum Cartesian surface-cell side in bohr */
    double target_area; /* Target planar patch area in bohr squared */
    double area_tolerance; /* Relative patch-area tolerance */
    double gaussian_scale; /* Gaussian width scale */
    double w_0ls_from; /* Lower low-gradient switch bound */
    double w_0ls_to; /* Upper low-gradient switch bound */
    double w_0ls_p; /* Low-gradient switch power */
    double w_0ls_a; /* Low-gradient switch shape */
    double weight_cutoff; /* Minimum retained switch weight */
    /* --- end of 1.0 layout --- */
} moist_mt_options;

moist_API_ENTRY void moist_API_CALL
moist_init_mt_options(moist_error error,
    moist_mt_options *options, size_t struct_size) moist_API_SUFFIX__V_1_0;

typedef struct {
    size_t struct_size;
    int nleb;
    bool do_fine;
    double tolerance;
    int proj_maxiter;
    int proj_level;
    double branch_weight_s;
    double rho_grid_h;
    int wleb_prune_level;
    /* --- end of 1.0 layout --- */
} moist_drop_options;

moist_API_ENTRY void moist_API_CALL
moist_init_drop_options(moist_error error, moist_drop_options *options,
                        size_t struct_size) moist_API_SUFFIX__V_1_0;

typedef struct {
    size_t struct_size;
    int nleb;
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
    /* --- end of 1.0 layout --- */
} moist_cfc_options;

moist_API_ENTRY void moist_API_CALL
moist_init_cfc_options(moist_error error, moist_cfc_options *options,
                        size_t struct_size) moist_API_SUFFIX__V_1_0;

typedef struct {
    size_t struct_size;
    double rho_iso;
    double scale;
    moist_isodensity_basis basis;  /* Only read by moist_new_isodensity_lsf */
    bool core_fill;  /* Fill ECP nuclei (nuclear charge z < atomic number) against density dips */
    double core_radius_fraction;  /* Fill radius as a fraction of the atom radius */
    /* --- end of 1.0 layout --- */
} moist_isodensity_options;

moist_API_ENTRY void moist_API_CALL
moist_init_isodensity_options(moist_error error, moist_isodensity_options *options,
                        size_t struct_size) moist_API_SUFFIX__V_1_0;

/// Current generated-PCM layout; rebuild callers using the old PCM options
typedef struct {
    size_t struct_size;
    moist_pcm_solver solver;
    int solver_maxiter;
    int solver_restart;
    double solver_tol;
    moist_pcm_storage storage;
    double newton_tol;
    int newton_maxiter;
    bool parameter_derivatives; /* return the "parameter_adjoint" response item */
} moist_pcm_options;

moist_API_ENTRY void moist_API_CALL
moist_init_pcm_options(moist_error error, moist_pcm_options *options,
                        size_t struct_size) moist_API_SUFFIX__V_1_0;

/// Numerical regularization of normalized GOSTSHYP gradient traces
typedef struct {
    size_t struct_size;
    double regularization_start;  /* |ftilde| at and below which a point is off, bohr^-4 */
    double regularization_end;    /* |ftilde| at and above which 1/ftilde is exact, bohr^-4 */
    bool suppress_negative_amplitudes;
    /* --- end of 1.0 layout --- */
} moist_gostshyp_options;

moist_API_ENTRY void moist_API_CALL
moist_init_gostshyp_options(moist_error error, moist_gostshyp_options *options,
                            size_t struct_size) moist_API_SUFFIX__V_1_0;

/// Banner text returned by moist_get_banner
typedef enum {
    moist_banner_full = 0, moist_banner_short = 1,
    moist_banner_ascii = 2, moist_banner_build = 3
} moist_banner;

/// Copy banner text, or query its length with buffer=NULL and capacity=0
/// length required, receives the byte count excluding the NUL terminator
/// capacity>0: buffer required, always NUL-terminated, even on error
/// Capacities above SIZE_MAX/2 rejected; supplied buffer still cleared
/// Truncation is successful: length >= capacity means a larger buffer needed
/// On error length unchanged; host chooses where to print
moist_API_ENTRY void moist_API_CALL
moist_get_banner(moist_error error, moist_banner style, char *buffer,
                 size_t capacity, size_t *length) moist_API_SUFFIX__V_1_0;

/// Create an LSF; NULL options selects defaults
moist_API_ENTRY moist_lsf moist_API_CALL
moist_new_svdw_lsf(moist_error error, const moist_svdw_options *options) moist_API_SUFFIX__V_1_0;
moist_API_ENTRY moist_lsf moist_API_CALL
moist_new_cfc_lsf(moist_error error, const moist_cfc_options *options) moist_API_SUFFIX__V_1_0;

/// Create a callback LSF; callback and context must outlive every cavity/model copy
moist_API_ENTRY moist_lsf moist_API_CALL
moist_new_isodensity_callback_lsf(moist_error error, moist_isodensity_lsf_callback callback,
                                  void *context, const moist_isodensity_options *options) moist_API_SUFFIX__V_1_0;

/// Create an LSF owning its Gaussian basis; shell_atom indices zero-based
/// exps and coeffs contain sum(shell_nprim) values; coefficients include
/// primitive normalization. Set the density on the constructed cavity, not
/// the source LSF
moist_API_ENTRY moist_lsf moist_API_CALL
moist_new_isodensity_lsf(moist_error error, int nshell, const int *shell_atom,
                        const int *shell_l, const int *shell_nprim,
                        const double *exps, const double *coeffs,
                        const moist_isodensity_options *options) moist_API_SUFFIX__V_1_0;
moist_API_ENTRY void moist_API_CALL
moist_delete_lsf(moist_lsf *lsf) moist_API_SUFFIX__V_1_0;

/* A model requires a context (NULL is an error); a cavity or component
 * takes one optionally. Every handle retains its context until deleted; the
 * context alone sets threads, verbosity and debug. A cavity or component
 * without one runs on its model's context and fails when used standalone */

/// Construct DROP from a required LSF; LSF and radii are copied
/// NULL radii selects CPCM radii; NULL options selects defaults
moist_API_ENTRY moist_cavity moist_API_CALL
moist_new_drop_cavity(moist_error error, moist_context context, moist_lsf lsf,
                      moist_radii radii,
                      const moist_drop_options *options) moist_API_SUFFIX__V_1_0;

/// Construct a frozen marching-tetrahedra surface, copying LSF and radii
/// NULL radii selects CPCM radii; NULL options selects native defaults
/// Supports every native LSF; geometry derivatives are unavailable
moist_API_ENTRY moist_cavity moist_API_CALL
moist_new_mt_cavity(moist_error error, moist_context context,
    moist_lsf lsf, moist_radii radii,
    const moist_mt_options *options) moist_API_SUFFIX__V_1_0;

/// Construct iSwiG, copying radii; NULL radii selects CPCM, NULL options defaults
moist_API_ENTRY moist_cavity moist_API_CALL
moist_new_iswig_cavity(moist_error error, moist_context context, moist_radii radii,
                       const moist_iswig_options *options) moist_API_SUFFIX__V_1_0;

/// Construct a model owning a cavity copy; add components before the first update
/// A cavity or component copy keeps its own context; one without runs on this
moist_API_ENTRY moist_model moist_API_CALL
moist_new_model(moist_error error, moist_context context,
                moist_cavity cavity) moist_API_SUFFIX__V_1_0;

/// Construct PCM components; NULL options selects the default solver
/// NULL context runs the component on its model's context
moist_API_ENTRY moist_component moist_API_CALL
moist_new_cpcm_component(moist_error error, moist_context context, double epsilon,
                         const moist_pcm_options *options) moist_API_SUFFIX__V_1_0;
moist_API_ENTRY moist_component moist_API_CALL
moist_new_cosmo_component(moist_error error, moist_context context, double epsilon,
                          const moist_pcm_options *options) moist_API_SUFFIX__V_1_0;
moist_API_ENTRY moist_component moist_API_CALL
moist_new_iefpcm_component(moist_error error, moist_context context, double epsilon,
                           const moist_pcm_options *options) moist_API_SUFFIX__V_1_0;
moist_API_ENTRY moist_component moist_API_CALL
moist_new_ssvpe_component(moist_error error, moist_context context, double epsilon,
                          const moist_pcm_options *options) moist_API_SUFFIX__V_1_0;
moist_API_ENTRY moist_component moist_API_CALL
moist_new_dpcm_component(moist_error error, moist_context context, double epsilon,
                          const moist_pcm_options *options) moist_API_SUFFIX__V_1_0;
moist_API_ENTRY moist_component moist_API_CALL
moist_new_dpcm_original_component(moist_error error, moist_context context, double epsilon,
                                   const moist_pcm_options *options) moist_API_SUFFIX__V_1_0;
moist_API_ENTRY moist_component moist_API_CALL
moist_new_lpb_component(moist_error error, moist_context context, double epsilon, double kappa,
                         const moist_pcm_options *options) moist_API_SUFFIX__V_1_0;
moist_API_ENTRY moist_component moist_API_CALL
moist_new_slic_component(moist_error error, moist_context context, double epsilon, double alpha,
                          double beta, double gamma, double mu, double phi_static,
                          const moist_pcm_options *options) moist_API_SUFFIX__V_1_0;
moist_API_ENTRY moist_component moist_API_CALL
moist_new_nonlocal_component(moist_error error, moist_context context, double eps_solvent, double eps_inf,
                              double lambda, double eps_solute,
                              const moist_pcm_options *options) moist_API_SUFFIX__V_1_0;

/*
 * Global API queries
**/

/// Library version as major * 10000 + minor * 100 + patch
moist_API_ENTRY int moist_API_CALL
moist_get_version(void) moist_API_SUFFIX__V_1_0;

/// Copy the full release version, including prerelease suffix
/// Buffer, length-query, truncation, and failure rules follow moist_get_banner
moist_API_ENTRY void moist_API_CALL
moist_get_version_string(moist_error error, char *buffer, size_t capacity,
                         size_t *length) moist_API_SUFFIX__V_1_0;

/*
 * Solvent descriptors -- bulk properties of one solvent table entry
 *
 * - Table values follow the Minnesota Solvent Descriptor Database: epsilon
 *   and surface tension at 298.15 K, refractive index at 293.15 K, density
 *   mostly at 293.15 K
 * - Derived values (molar mass, number density, molecular volume) come from
 *   the stored solvent geometry
 * - A quiet NaN marks a value the table lacks and carries into the values
 *   derived from it; 0 is a real value (e.g. abraham_alpha of aprotic
 *   solvents). Every table value is present; decalin (ID 94) has no stored
 *   geometry, so its six geometry-derived values are NaN
 * - IDs are frozen: an ID never changes meaning, new solvents are appended
 *
 * Size contract
 * - Pass sizeof the struct as struct_size; smaller than the full 1.0 layout
 *   is an error
 * - moist writes min(struct_size, its own layout) bytes and records that
 *   count in descriptors->struct_size; fields beyond it are untouched
 * - No initialization needed; on any error *descriptors is untouched
 *
 * Need no context; pure table reads, safe on any thread
 */
typedef struct {
    size_t struct_size;            /* bytes filled by moist */
    int id;                        /* frozen table ID */
    char name[MOIST_NAME_MAX + 1]; /* canonical name, lowercase */
    double epsilon;                /* static relative permittivity */
    double refractive_index;       /* sodium D line */
    double abraham_alpha;          /* Abraham H-bond acidity */
    double abraham_beta;           /* Abraham H-bond basicity */
    double surface_tension_si;     /* macroscopic, N/m */
    double surface_tension_au;     /* macroscopic, Eh/bohr^2 */
    double mass_density_si;        /* kg/m^3 */
    double mass_density_au;        /* m_e/bohr^3 */
    double molar_mass_si;          /* kg/mol */
    double molecular_mass_au;      /* mass of one molecule, m_e */
    double number_density_si;      /* 1/m^3 */
    double number_density_au;      /* 1/bohr^3 */
    double molecular_volume_si;    /* volume per molecule, m^3 */
    double molecular_volume_au;    /* volume per molecule, bohr^3 */
    /* --- end of 1.0 layout --- */
} moist_solvent_descriptors;

/// Descriptors of a solvent by name or any alias (e.g. "H2O", "DMSO")
/// Case-insensitive, surrounding whitespace ignored. Errors: NULL name or
/// descriptors, struct_size below the 1.0 layout, blank or unknown name
moist_API_ENTRY void moist_API_CALL
moist_get_solvent_descriptors(moist_error error, const char *name,
                              moist_solvent_descriptors *descriptors,
                              size_t struct_size) moist_API_SUFFIX__V_1_0;

/// Descriptors of a solvent by its table ID; errors as by name, and an
/// unknown ID
moist_API_ENTRY void moist_API_CALL
moist_get_solvent_descriptors_by_id(moist_error error, int id,
                                    moist_solvent_descriptors *descriptors,
                                    size_t struct_size) moist_API_SUFFIX__V_1_0;

/// Copy the "Solvent properties" printout of filled descriptors, or query its
/// length, following moist_get_banner
/// Formats the values passed, without a table lookup; struct_size below the
/// 1.0 layout is an error. The host chooses where to print
moist_API_ENTRY void moist_API_CALL
moist_get_solvent_descriptors_text(moist_error error,
                                   const moist_solvent_descriptors *descriptors,
                                   char *buffer, size_t capacity,
                                   size_t *length) moist_API_SUFFIX__V_1_0;

/*
 * Error handle class
**/

/// Create new error handle object
moist_API_ENTRY moist_error moist_API_CALL
moist_new_error(void) moist_API_SUFFIX__V_1_0;

/// moist_success, moist_failure, or moist_invalid_error for a NULL handle
moist_API_ENTRY int moist_API_CALL
moist_check_error(moist_error error) moist_API_SUFFIX__V_1_0;

/// Copy the diagnostic; no error produces an empty string, NULL error a diagnostic
/// buffer and buffersize required; *buffersize must be positive
/// Writes at most *buffersize bytes, always NUL-terminating; long messages truncate
moist_API_ENTRY void moist_API_CALL
moist_get_error(moist_error error,
                char* buffer,
                const int* buffersize) moist_API_SUFFIX__V_1_0;

/// Delete error handle object
moist_API_ENTRY void moist_API_CALL
moist_delete_error(moist_error* error) moist_API_SUFFIX__V_1_0;

/* ERROR CONVENTION
 *
 * - Fallible operations report through a required, non-NULL moist_error handle
 * - Each fallible call replaces the previous diagnostic
 * - Check moist_check_error after each call, before using any result
 * - Queries write values through output parameters; values never encode
 *   error status
 * - Required pointers checked; numeric outputs unchanged on error
 * - Banner, version, and field-description getters clear writable buffers
 *   on error
 * - Constructors return a handle and report failure through the error
 *   handle; release owned handles with their matching delete function
 */

/*
 * Molecular structure data class
**/

/// Create new molecular structure data (quantities in Bohr)
/// The structure holds element numbers, geometry and the net charge; the
/// nuclear charges of the host Hamiltonian arrive through the
/// "nuclear_charges" coupling request
moist_API_ENTRY moist_structure moist_API_CALL
moist_new_structure(moist_error error,
                    int natoms,
                    const int* numbers /* [natoms] */,
                    const double* positions /* : row-major (natoms, 3) */,
                    const double* lattice /* : row-major (3, 3), vectors in rows */,
                    const bool* periodic /* [3] */,
                    const double* charge /* : net charge in e, NULL for neutral */) moist_API_SUFFIX__V_1_0;

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

/// Create Gaussian charge scheme radii model (includes the 1.2 scale factor)
moist_API_ENTRY moist_radii moist_API_CALL
moist_new_gauss_radii(moist_error error) moist_API_SUFFIX__V_1_0;

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


/// Create a pressure-volume energy component equal to `pressure * cavity volume`
moist_API_ENTRY moist_component moist_API_CALL
moist_new_pv_component(moist_error error, moist_context context,
                       double pressure) moist_API_SUFFIX__V_1_0;

/// Create a GOSTSHYP hydrostatic-pressure component at `pressure` in
/// Hartree/bohr^3. Cannot form its own density traces: answer the
/// Gaussian-moment request of the coupling ("gt", "pt", "mt", "rt") with
/// unit-integral G=(omega/pi)^1.5 exp(-omega |r-C|^2). Use the same
/// normalization for every moment and for AO response/gradient integrals. Use
/// moist_answer_coupling_request in every phase, and contract the
/// "gaussian_amplitude" item of the response walk. NULL options selects
/// compiled defaults
moist_API_ENTRY moist_component moist_API_CALL
moist_new_gostshyp_component(moist_error error, moist_context context, double pressure,
                             const moist_gostshyp_options *options) moist_API_SUFFIX__V_1_0;

/// Delete a standalone solvation-model component handle
moist_API_ENTRY void moist_API_CALL
moist_delete_component(moist_component* component) moist_API_SUFFIX__V_1_0;


/// Append a copied component to a continuum model before its first update
moist_API_ENTRY void moist_API_CALL
moist_add_model_component(moist_error error, moist_model model, moist_component component) moist_API_SUFFIX__V_1_0;

/* HOST COUPLING PROTOCOL
 *
 * The same steps as the Fortran interface:
 *
 * 1. Update the model with the structure, moist_update_model(err, model, mol);
 *    it does everything that needs no host state. Mint a coupling with
 *    moist_new_coupling, before or after; it survives later updates
 * 2. Build phase, part of every host loop, after every update and for every
 *    new density: moist_prepare_model_build(err, model, cpl) clears every
 *    answer of cpl and stages the requests the build reads; walk and answer
 *    them as in 4 and 5; moist_build_model(err, model, cpl) completes the
 *    model. On a geometric cavity the walk is empty and the build changes
 *    nothing, so a host never branches on the cavity type. Build-phase
 *    requests:
 *    - "density_matrix" ("dcart", nbasis x nbasis, symmetric) of an internal
 *      isodensity cavity, in the cartesian-monomial basis of
 *      moist_get_isodensity_cart_layout
 *    - "nuclear_charges" ("z", one per atom) of an isodensity cavity with the
 *      core fill on: the nuclear charges of the host Hamiltonian, atomic
 *      numbers minus the core electrons an ECP replaces; an all-electron host
 *      answers the atomic numbers
 *    The callback isodensity level set is pulled point by point during the
 *    surface construction and stays outside the coupling. Without core fill,
 *    callbacks build during update; repeat the update for each new density.
 *    A standalone cavity runs
 *    moist_update_cavity, moist_prepare_cavity_build and moist_build_cavity
 *    on a coupling from moist_new_cavity_coupling the same way
 * 3. Stage a phase with moist_prepare_model_energy/response/gradient; energy
 *    clears all answers except those of the build phase, response and
 *    gradient keep valid answers; every prepare starts a new walk over the
 *    requests. Before the build a deferring cavity holds no surface: no
 *    component declares, grid outputs are refused for want of a point count
 *    and get_* fails on its own.
 *    PCM components ask for "nuclear_charges" in the gradient phase unless
 *    the build phase already holds them
 * 4. Walk the requests that still miss an output with
 *    `while (moist_next_coupling_request(err, cpl))`. Each such request is
 *    visited once per pass; false ends the pass and rewinds, so a second loop
 *    retries what is still missing; false also reports a failure, check err
 *    after the loop. Leaving the loop early resumes the pass at the next call
 * 5. For the current request, read its name, ask which outputs are missing,
 *    read the inputs it lists (moist_get_coupling_request_field_real, e.g.
 *    the "width" of "gaussian_moments"), compute the outputs on the model's
 *    grid -- read it with moist_get_model_field_real ("xyz", "xi0", ...),
 *    which works for every model family -- and submit each with
 *    moist_answer_coupling_request. A rejected output stays missing until a
 *    valid retry; the others survive
 * 6. Call the moist_get_model_* entry matching the staged phase; fails by
 *    name on a wrong staging or any missing output, evaluation does not
 *    consume answers
 * 7. Walk the response with `while (moist_next_response_item(err, resp))`:
 *    every item present is visited once per pass. Read the current item's
 *    name, copy its arrays with moist_get_response_field_real (their shapes
 *    from moist_get_response_field_info) and contract them; stop on an item
 *    the host cannot contract. An absent item is physics, not an error
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
 *             moist_get_response_field_real(err, resp, "w_phi", w_phi);
 *             host_fock_potential(ngrid, w_phi, fock);
 *         } else {
 *             host_abort("unsupported response item", name);
 *         }
 *     }
 *
 * The "density" item is produced by the response phase only. Its weights are
 * dE/drho at fixed nuclei -- the same object in both phases; moist forms them
 * once and does not contract them again in the gradient phase. On a cavity
 * whose surface follows the density, run the response phase before the
 * gradient and reuse those weights: contracted with d rho/dP they complete
 * the Fock matrix, and with the basis-center derivative d rho/dR they
 * complete the nuclear gradient. Going straight to the gradient phase there
 * drops that term silently
 *
 * - A model update invalidates the answers of all its couplings, even if the
 *   grid size is unchanged; a build invalidates every coupling but the one it
 *   read, which keeps its build-phase answers until the next prepare_build
 * - Host must prepare energy again for a new density, and also run the
 *   build phase of a density-dependent cavity before preparing
 * - Release couplings before their model or cavity
 * - Response arrays are independent copies; get_response/get_gradient clear
 *   the response once they accept the request, which also starts a new walk
 *   over its items; a rejected request, including an unimplemented theory,
 *   leaves it untouched
 */
/// Mint a borrowed coupling; needs an updated model
/// Several couplings may coexist; every component declares its requests
/// Coupling survives later model updates; the prepare_* entry points
/// refresh it. Grid outputs are refused until the model is built
moist_API_ENTRY moist_coupling moist_API_CALL
moist_new_coupling(moist_error error,
                   moist_model model) moist_API_SUFFIX__V_1_0;

/// Mint a borrowed coupling of a standalone cavity for its build phase
/// Release it with moist_delete_coupling before the cavity. A model-owned
/// cavity handle refuses: its model's coupling answers the build phase
moist_API_ENTRY moist_coupling moist_API_CALL
moist_new_cavity_coupling(moist_error error,
                          moist_cavity cavity) moist_API_SUFFIX__V_1_0;

/// Delete a coupling handle and release its model- or cavity-owned collection
moist_API_ENTRY void moist_API_CALL
moist_delete_coupling(moist_coupling* coupling) moist_API_SUFFIX__V_1_0;

/// Create an empty response handle
moist_API_ENTRY moist_response moist_API_CALL
moist_new_response(moist_error error) moist_API_SUFFIX__V_1_0;

/// Delete a response handle
moist_API_ENTRY void moist_API_CALL
moist_delete_response(moist_response* response) moist_API_SUFFIX__V_1_0;

/// Stage the build phase: clear every answer of the coupling and start a new
/// walk over the requests the build reads; the atom count comes from the
/// structure of the latest moist_update_model. A model whose cavity waits for
/// its build phase is not built until the next build succeeds, and every
/// coupling of the model loses its answers over it; the cavity keeps its last
/// built surface until that build replaces it
moist_API_ENTRY void moist_API_CALL
moist_prepare_model_build(moist_error error,
                          moist_model model,
                          moist_coupling coupling) moist_API_SUFFIX__V_1_0;

/// Complete a model from the build phase of a coupling staged by
/// moist_prepare_model_build. Fails by name on a wrong staging or a missing
/// required output. The coupling keeps its build-phase answers; every other
/// coupling of the model is invalidated. Only a deferred cavity and its
/// components are rebuilt
moist_API_ENTRY void moist_API_CALL
moist_build_model(moist_error error,
                  moist_model model,
                  moist_coupling coupling) moist_API_SUFFIX__V_1_0;

/// Stage the energy phase: mark every answer but those of the build phase
/// stale and start a new walk. On a cavity without a surface the snapshot
/// has no point count, so grid outputs are refused
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
/// staged for the energy phase and on an unanswered required output,
/// including one whose answer was rejected. Adds the energy contribution to
/// *energy; initialize it before the first call. On failure the accumulator
/// is unchanged
moist_API_ENTRY void moist_API_CALL
moist_get_model_energy(moist_error error,
                       moist_model model,
                       moist_coupling coupling,
                       double* energy) moist_API_SUFFIX__V_1_0;

/// Host part of the response phase from a staged coupling. `response` is
/// cleared once the request is accepted and returns every item the model
/// produces; a rejected request, including an unimplemented theory, leaves it
/// untouched
moist_API_ENTRY void moist_API_CALL
moist_get_model_response(moist_error error,
                         moist_model model,
                         moist_coupling coupling,
                         moist_response response) moist_API_SUFFIX__V_1_0;

/// Nuclear gradient from a staged coupling, plus the host part of the gradient
/// phase in `response` (cleared once the request is accepted; a rejected
/// request, including an unimplemented theory, leaves it untouched): the
/// potential adjoint and the Gaussian amplitudes. Not the "density" item -- see
/// HOST COUPLING PROTOCOL above for the weights a density-backed cavity carries
/// over from the response phase. `gradient` is row-major (nat_cap, 3); capacity
/// checked against the atom count before anything is written. Adds to gradient;
/// initialize logical entries before the first call. On failure the gradient
/// accumulator is unchanged; padding is untouched
moist_API_ENTRY void moist_API_CALL
moist_get_model_gradient(moist_error error,
                         moist_model model,
                         moist_coupling coupling,
                         moist_response response,
                         int nat_cap,
                         double* gradient /* : row-major (nat_cap, 3) */) moist_API_SUFFIX__V_1_0;

/// Advance to the next request with a missing output of the staged phase
/// Returns false after the last one, rewinding so that the next call starts a
/// new pass, and on any failure, so it can drive a while loop: check error
/// afterwards. Leaving the loop early resumes the pass at the next call
moist_API_ENTRY bool moist_API_CALL
moist_next_coupling_request(moist_error error, moist_coupling coupling) moist_API_SUFFIX__V_1_0;
/// NUL-terminated canonical name of the current request, written into
/// name[MOIST_NAME_MAX + 1]. An error when no request is current
moist_API_ENTRY void moist_API_CALL
moist_get_coupling_request_name(moist_error error, moist_coupling coupling, char* name) moist_API_SUFFIX__V_1_0;
/// Whether the named output of the current request is required by the staged
/// phase and still missing, evaluated at the call rather than frozen when the
/// cursor arrived -- so answering an output and asking again reports it as no
/// longer missing. A name the request does not declare is not missing. An
/// error when no request is current
moist_API_ENTRY void moist_API_CALL
moist_get_coupling_request_missing(moist_error error, moist_coupling coupling, const char* output,
    bool* missing) moist_API_SUFFIX__V_1_0;

/// Submit one named output of the current request. values is row-major
/// (ngrid, dims...) on the cavity grid, with the output's leading extents:
/// phi, dphi_dxi, gt: [ngrid]; dphi_dr, d2phi_dr_dxi, pt, rt: [ngrid][3];
/// d2phi_dr2, mt: [ngrid][3][3].
/// per atom: q, z: [natom], mu: [natom][3], theta: [natom][3][3]; over the
/// basis: dcart: [nbasis][nbasis], symmetric, with nbasis the request's own
/// input "nbasis". moist reads exactly that many
/// values. An error when no request is current,
/// for an unknown output, a NULL buffer or a non-finite value. A rejected
/// output stays missing until a valid retry; other outputs are unaffected
moist_API_ENTRY void moist_API_CALL
moist_answer_coupling_request(moist_error error, moist_coupling coupling, const char* output,
    const double* values) moist_API_SUFFIX__V_1_0;

/// Named inputs of the current request
///
/// A request kind may carry inputs its component chose, which the host reads
/// and never recomputes: "gaussian_moments" lists "width" [ngrid], the
/// Gaussian moment exponents in bohr**-2 -- not a cavity field; "density_matrix"
/// lists "nbasis", a scalar, both extents of "dcart". A kind without inputs
/// lists none. Enumerate them with
/// moist_get_coupling_request_field_count +
/// moist_get_coupling_request_field_info, then read one by name; every input
/// is MOIST_FIELD_REAL, and `dims` lists its extents slowest-varying first,
/// as for the cavity fields. An error when no request is current

/// Number of named inputs of the current request
moist_API_ENTRY void moist_API_CALL
moist_get_coupling_request_field_count(moist_error error,
                                       moist_coupling coupling,
                                       int* nfield) moist_API_SUFFIX__V_1_0;

/// Describe one input of the current request by position; outputs as
/// moist_get_cavity_field_info: `dims` receives MOIST_FIELD_MAX_RANK extents
/// with the slowest-varying one first (unused entries are 1), and `count` is
/// the number of values a read writes
moist_API_ENTRY void moist_API_CALL
moist_get_coupling_request_field_info(moist_error error,
                                      moist_coupling coupling,
                                      int index /* : 0-based, below the count from
                                                   moist_get_coupling_request_field_count */,
                                      char* name /* [MOIST_FIELD_NAME_MAX + 1] */,
                                      int* dtype /* : MOIST_FIELD_REAL */,
                                      int* rank /* : 0 for a scalar */,
                                      int* dims /* [MOIST_FIELD_MAX_RANK] */,
                                      int* count /* : elements a read writes */) moist_API_SUFFIX__V_1_0;

/// Copy an input description of the current request or query its length,
/// following moist_get_banner
moist_API_ENTRY void moist_API_CALL
moist_get_coupling_request_field_about(moist_error error, moist_coupling coupling,
                                       const char* name, char* about,
                                       size_t capacity, size_t* length) moist_API_SUFFIX__V_1_0;

/// Read one input of the current request by name: the `count` values
/// moist_get_coupling_request_field_info reports, row-major with its `dims`.
/// An error for a name the request does not list or a NULL buffer; nothing
/// is written then
moist_API_ENTRY void moist_API_CALL
moist_get_coupling_request_field_real(moist_error error,
                                      moist_coupling coupling,
                                      const char* name,
                                      double* values /* [count] */) moist_API_SUFFIX__V_1_0;

/// Advance to the next item of the response. Every item present is visited
/// once per pass. Returns false after the last one, rewinding so that the next
/// call starts a new pass, at once for an empty response, and on any failure,
/// so it can drive a while loop: check error afterwards
moist_API_ENTRY bool moist_API_CALL
moist_next_response_item(moist_error error, moist_response response) moist_API_SUFFIX__V_1_0;
/// NUL-terminated canonical name of the current item, written into
/// name[MOIST_NAME_MAX + 1]. An error when no item is current
moist_API_ENTRY void moist_API_CALL
moist_get_response_item_name(moist_error error, moist_response response,
    char* name) moist_API_SUFFIX__V_1_0;

/// Named arrays of the current response item
///
/// Every item describes its arrays like a cavity its fields: enumerate them
/// with moist_get_response_field_count + moist_get_response_field_info, then
/// read one by name with moist_get_response_field_real. An array the item was
/// accumulated without is not listed. Every array is MOIST_FIELD_REAL, and
/// `dims` lists its extents slowest-varying first, the reverse of the native
/// Fortran shape: "w_grad_rho" (3, ngrid) reports {ngrid, 3}, "dg_dphi"
/// (ngrid, natom) reports {natom, ngrid}. An error when no item is current.
/// The items, their arrays and the host's contraction:
///   "potential_adjoint": "w_phi" [ngrid], w_phi_i = dE/dphi_i, contracted by
///     the host as `F_uv += sum_i w_phi_i V_uv(r_i)`; for a stationary PCM
///     this is the surface charge q_i
///     "w_dphi_dr" [ngrid][3] weights for the raw spatial potential derivative;
///     present only when computed. Contract them as
///     `F_uv += sum_i,a w_dphi_dr[i][a] dV_uv(r_i)/dr_a`
///   "density": "w_rho" [ngrid], "w_grad_rho" [ngrid][3], "w_hess_rho"
///     [ngrid][3][3], the weights conjugate to the density value, gradient and
///     Hessian. Hessian weights use C indices [point][b][a] for native
///     weights(a,b,point); need not be symmetric, preserve this axis order
///     when contracting
///   "gaussian_amplitude": "w_overlap" [ngrid], "w_normal_deriv" [ngrid],
///     contracted as `F_uv += sum_i [w_overlap[i] g_uv,i + w_normal_deriv[i] f_uv,i]`
///     with g_uv,i = <u|G_i|v> and f_uv,i = n_i . grad_r g_uv,i
///   "atomic_multipole_adjoint": "dg_dq" [natom], "dg_dmu" [natom][3],
///     "dg_dtheta" [natom][3][3], dE/dq_A, dE/dmu_A and dE/dtheta_A of the
///     solute's atomic multipoles, contracted as
///     `F_uv += sum_A [dg_dq[A] dq_A/dP_uv + dg_dmu[A] . dmu_A/dP_uv + dg_dtheta[A] : dtheta_A/dP_uv]`
///     and with d/dR of the same multipoles for the gradient. Quadrupole
///     weights use C indices [atom][b][a] for native weights(a,b,atom); an
///     order the model never consumed is not listed
///   "atomic_charge_adjoint": "dg_dq" [natom], dE/dq_A, contracted as
///     `F_uv += sum_A dg_dq[A] dq_A/dP_uv` and with dq_A/dR for the gradient
///   "radial_potential_adjoint": "dg_dphi" [natom][ngrid], dE/dphi_a(r_i) of
///     the site-resolved radial potential, contracted as
///     `F_uv += sum_a sum_i dg_dphi[a][i] dphi_a(r_i)/dP_uv` and with
///     dphi_a(r_i)/dR for the gradient
///   "parameter_adjoint": scalar fields, dE/dp of the physical model parameters
///     at fixed geometry and traces: "dg_deps", "dg_deps_inf", "dg_dlambda",
///     "dg_deps_solute" (nonlocal), "dg_dkappa" (LPB), "dg_dalpha", "dg_dbeta",
///     "dg_dgamma", "dg_dmu", "dg_dphi_static" (SLIC); a parameter the model
///     lacks is not listed. Present only when the PCM setting
///     `parameter_derivatives` is true (default false), and only in the response
///     phase, never in a gradient read; a host that stops on unknown items must
///     not enable it

/// Number of named arrays of the current response item
moist_API_ENTRY void moist_API_CALL
moist_get_response_field_count(moist_error error,
                               moist_response response,
                               int* nfield) moist_API_SUFFIX__V_1_0;

/// Describe one array of the current item by position; outputs as
/// moist_get_cavity_field_info: `dims` receives MOIST_FIELD_MAX_RANK extents
/// with the slowest-varying one first (unused entries are 1), and `count` is
/// the number of values a read writes
moist_API_ENTRY void moist_API_CALL
moist_get_response_field_info(moist_error error,
                              moist_response response,
                              int index /* : 0-based, below the count from
                                           moist_get_response_field_count */,
                              char* name /* [MOIST_FIELD_NAME_MAX + 1] */,
                              int* dtype /* : MOIST_FIELD_REAL */,
                              int* rank /* : 0 for a scalar */,
                              int* dims /* [MOIST_FIELD_MAX_RANK] */,
                              int* count /* : elements a read writes */) moist_API_SUFFIX__V_1_0;

/// Copy an array description of the current item or query its length,
/// following moist_get_banner
moist_API_ENTRY void moist_API_CALL
moist_get_response_field_about(moist_error error, moist_response response,
                               const char* name, char* about,
                               size_t capacity, size_t* length) moist_API_SUFFIX__V_1_0;

/// Read one array of the current item by name: the `count` values
/// moist_get_response_field_info reports, row-major with its `dims`. An
/// error for an array the current item does not list or a NULL buffer;
/// nothing is written then
moist_API_ENTRY void moist_API_CALL
moist_get_response_field_real(moist_error error,
                              moist_response response,
                              const char* name,
                              double* values /* [count] */) moist_API_SUFFIX__V_1_0;

/// Update a solvation model with a molecular structure: everything that needs
/// no host state; a model whose cavity waits for host answers completes in
/// the build phase. Invalidates every coupling of the model
moist_API_ENTRY void moist_API_CALL
moist_update_model(moist_error error,
                   moist_model model,
                   moist_structure mol) moist_API_SUFFIX__V_1_0;

/// Get a borrowed cavity handle from a solvation model
/// Valid as long as the parent model exists, but cannot be rebuilt
/// independently: moist_update_cavity, moist_prepare_cavity_build,
/// moist_build_cavity and moist_new_cavity_coupling reject borrowed handles
/// Use moist_delete_cavity to release the handle (does NOT destroy the
/// model's cavity); the handle retains the model's context until then
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





/// Query the internal isodensity cartesian-component layout (DROP-specific)
/// Two-pass: call with the array pointers NULL to read ncart/nshell, then
/// allocate and call again. `shell_cart_offset` is 0-based (length nshell+1);
/// `comp_lx/ly/lz` are the per-component monomial powers (length ncart)
/// defining the ordering the host must match when building its density
/// transform; the "density_matrix" request of the build phase runs over it.
/// With moist_isodensity_basis_spherical, a density outside the span of the
/// spherical functions (a host with cartesian basis functions) is refused.
/// Any of the output pointers may be NULL to skip it
moist_API_ENTRY void moist_API_CALL
moist_get_isodensity_cart_layout(moist_error error,
                                 moist_cavity cavity,
                                 int* ncart,
                                 int* nshell,
                                 int* shell_cart_offset /* [nshell+1] */,
                                 int* comp_lx /* [ncart] */,
                                 int* comp_ly /* [ncart] */,
                                 int* comp_lz /* [ncart] */) moist_API_SUFFIX__V_1_0;

/// Master numerical tolerance currently configured on a DROP cavity
/// Reports an error for non-DROP cavity types
moist_API_ENTRY void moist_API_CALL
moist_get_drop_cavity_tolerance(moist_error error,
                                moist_cavity cavity,
                                double* tolerance) moist_API_SUFFIX__V_1_0;

/*
 * Generic cavity operations (Tier 1 - work on all cavity types)
**/

/// Generic update cavity - works for all cavity types
/// Everything that needs no host state; a geometric cavity is built
/// afterwards, an isodensity cavity waiting for host answers completes in
/// the build phase. Invalidates every coupling of the cavity
moist_API_ENTRY void moist_API_CALL
moist_update_cavity(moist_error error,
                    moist_cavity cavity,
                    moist_structure mol) moist_API_SUFFIX__V_1_0;

/// Stage the build phase of a standalone cavity: clear every answer of the
/// coupling (minted by moist_new_cavity_coupling) and start a new walk over
/// the requests the build reads. A cavity that waits for its build phase
/// drops its surface until the next build succeeds; its surface readers
/// refuse it
moist_API_ENTRY void moist_API_CALL
moist_prepare_cavity_build(moist_error error,
                           moist_cavity cavity,
                           moist_coupling coupling) moist_API_SUFFIX__V_1_0;

/// Complete a cavity of any type from the build phase of a coupling staged by
/// moist_prepare_cavity_build. Fails by name on a wrong staging or a missing
/// required output. The coupling keeps its build-phase answers; every other
/// coupling of the cavity is invalidated. Only a deferred surface is rebuilt
moist_API_ENTRY void moist_API_CALL
moist_build_cavity(moist_error error,
                   moist_cavity cavity,
                   moist_coupling coupling) moist_API_SUFFIX__V_1_0;

/// Get generic cavity sizes - works for all cavity types
/// Returns ngrid (number of grid points) and nsph (number of spheres)
/// Refuses a cavity that holds no surface, naming the build phase; its named
/// fields then list only nsph, radii and sphxyz
moist_API_ENTRY void moist_API_CALL
moist_get_cavity_sizes(moist_error error,
                       moist_cavity cavity,
                       int* ngrid,
                       int* nsph) moist_API_SUFFIX__V_1_0;

/// Get generic cavity results - works for all cavity types
/// Returns fields shared by all cavity types. `converged` contains DROP's
/// per-point projection status; cavity types without a projection solve return
/// true for every retained point
/// Call moist_get_cavity_sizes first to get ngrid, nsph, then allocate arrays
/// and pass the capacities you allocated them with (see the capacity
/// convention at the top of this header). Nothing is written when either
/// capacity is smaller than the current size; an API error is set instead
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
/// Enumerate with moist_get_cavity_field_count + moist_get_cavity_field_info,
/// then read one with the accessor matching its type tag. A field that was not
/// computed is not listed; nothing is written to a rejected buffer. Rules and
/// conventions: docs/reference/fields.rst

/// Element type tags reported by the moist_get_*_field_info entries. Values
/// are part of the contract; read a field with the accessor matching its tag
#define MOIST_FIELD_REAL 1
#define MOIST_FIELD_INT 2
#define MOIST_FIELD_BOOL 3

/// Highest rank a field can have, i.e. the length of the `dims` buffer
/// the *_field_info entries write
#define MOIST_FIELD_MAX_RANK 3
/// Maximum field-name length excluding the NUL terminator. Mirrored by
/// `max_field_name_len` (api.f90)
#define MOIST_FIELD_NAME_MAX 64

/// Number of named result fields the cavity currently holds
/// Count changes when the cavity is rebuilt or its property request
/// changes; read it again after moist_update_cavity
moist_API_ENTRY void moist_API_CALL
moist_get_cavity_field_count(moist_error error,
                             moist_cavity cavity,
                             int* nfield) moist_API_SUFFIX__V_1_0;

/// Describe one readable field by position
/// `index` is 0-based and must be below the reported count. `dtype` is one of
/// the MOIST_FIELD_* tags, `rank` is 0 for a scalar, `dims` receives
/// MOIST_FIELD_MAX_RANK extents with the slowest-varying one first (unused
/// entries are 1), and `count` is the number of elements a read writes --
/// allocate that many for it. Name written with its terminator into
/// name[MOIST_FIELD_NAME_MAX + 1]
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

/// Copy a field description or query its length, following moist_get_banner
/// Pass about=NULL and capacity=0 to query length (excluding the terminator)
moist_API_ENTRY void moist_API_CALL
moist_get_cavity_field_about(moist_error error, moist_cavity cavity,
                             const char* name, char* about,
                             size_t capacity, size_t* length) moist_API_SUFFIX__V_1_0;

/// Read a MOIST_FIELD_REAL field by name: the `count` elements
/// moist_get_cavity_field_info reports; a rank-2 field is written flat in
/// row-major order per the layout convention at the top of this header
moist_API_ENTRY void moist_API_CALL
moist_get_cavity_field_real(moist_error error,
                            moist_cavity cavity,
                            const char* name,
                            double* values /* [count] */) moist_API_SUFFIX__V_1_0;

/// Read a MOIST_FIELD_INT field by name
/// Sphere indices handed out 0-based, matching moist_get_cavity_results;
/// the DROP point ids `numbering` and `anchor_id` and the 1-based `branch` are
/// passed through as the cavity stores them
moist_API_ENTRY void moist_API_CALL
moist_get_cavity_field_int(moist_error error,
                           moist_cavity cavity,
                           const char* name,
                           int* values /* [count] */) moist_API_SUFFIX__V_1_0;

/// Read a MOIST_FIELD_BOOL field by name
moist_API_ENTRY void moist_API_CALL
moist_get_cavity_field_bool(moist_error error,
                            moist_cavity cavity,
                            const char* name,
                            bool* values /* [count] */) moist_API_SUFFIX__V_1_0;

/// Named fields of a model's evaluation domain (Tier 2)
///
/// The cavity field getters through a model handle, for every model family: a
/// continuum model reports its cavity's fields, a volume model its grid. No
/// update is required. See docs/reference/fields.rst

/// Number of named fields of the model's evaluation domain
/// Read it again after moist_update_model
moist_API_ENTRY void moist_API_CALL
moist_get_model_field_count(moist_error error,
                            moist_model model,
                            int* nfield) moist_API_SUFFIX__V_1_0;

/// Describe one model field by position; outputs as moist_get_cavity_field_info
moist_API_ENTRY void moist_API_CALL
moist_get_model_field_info(moist_error error,
                           moist_model model,
                           int index /* : 0-based, below the count from
                                        moist_get_model_field_count */,
                           char* name /* [MOIST_FIELD_NAME_MAX + 1] */,
                           int* dtype /* : one of MOIST_FIELD_* */,
                           int* rank /* : 0 for a scalar */,
                           int* dims /* [MOIST_FIELD_MAX_RANK] */,
                           int* count /* : elements a read writes */) moist_API_SUFFIX__V_1_0;

/// Copy a model field description or query its length, following
/// moist_get_banner
moist_API_ENTRY void moist_API_CALL
moist_get_model_field_about(moist_error error, moist_model model,
                            const char* name, char* about,
                            size_t capacity, size_t* length) moist_API_SUFFIX__V_1_0;

/// Read a MOIST_FIELD_REAL model field by name: the `count` elements
/// moist_get_model_field_info reports, rank-2 fields flat in row-major order
moist_API_ENTRY void moist_API_CALL
moist_get_model_field_real(moist_error error,
                           moist_model model,
                           const char* name,
                           double* values /* [count] */) moist_API_SUFFIX__V_1_0;

/// Read a MOIST_FIELD_INT model field by name; indices as the owner declares
/// them, e.g. a cavity's 0-based `owner`
moist_API_ENTRY void moist_API_CALL
moist_get_model_field_int(moist_error error,
                          moist_model model,
                          const char* name,
                          int* values /* [count] */) moist_API_SUFFIX__V_1_0;

/// Read a MOIST_FIELD_BOOL model field by name
moist_API_ENTRY void moist_API_CALL
moist_get_model_field_bool(moist_error error,
                           moist_model model,
                           const char* name,
                           bool* values /* [count] */) moist_API_SUFFIX__V_1_0;

/// Components of a continuum model and their named results (Tier 2)
///
/// Components are addressed by their 0-based position in the order they were
/// added; a component added twice repeats its name. Count and names need no
/// update. Each component lists "energy" (MOIST_FIELD_REAL, rank 0), its
/// contribution to the latest moist_get_model_energy in Hartree, after a
/// successful evaluation; moist_update_model clears it, and a failed
/// evaluation leaves it unavailable for the failing and later components.
/// A model of another family and an index outside the count are errors.
/// See docs/reference/fields.rst

/// Number of components of a continuum model
moist_API_ENTRY void moist_API_CALL
moist_get_model_component_count(moist_error error,
                                moist_model model,
                                int* ncomponents) moist_API_SUFFIX__V_1_0;

/// Copy a component name or query its length, following moist_get_banner
moist_API_ENTRY void moist_API_CALL
moist_get_model_component_name(moist_error error, moist_model model,
                               int component /* : 0-based, below the count from
                                                moist_get_model_component_count */,
                               char* name, size_t capacity,
                               size_t* length) moist_API_SUFFIX__V_1_0;

/// Copy a one-line component description, e.g. "Pressure-volume work,
/// pressure times cavity volume", or query its length, following
/// moist_get_banner
moist_API_ENTRY void moist_API_CALL
moist_get_model_component_description(moist_error error, moist_model model,
                                      int component /* : 0-based, below the count from
                                                       moist_get_model_component_count */,
                                      char* description, size_t capacity,
                                      size_t* length) moist_API_SUFFIX__V_1_0;

/// Number of named results the component currently holds
moist_API_ENTRY void moist_API_CALL
moist_get_model_component_field_count(moist_error error,
                                      moist_model model,
                                      int component,
                                      int* nfield) moist_API_SUFFIX__V_1_0;

/// Describe one component result by position
moist_API_ENTRY void moist_API_CALL
moist_get_model_component_field_info(moist_error error,
                                     moist_model model,
                                     int component,
                                     int index /* : 0-based, below the count from
                                                  moist_get_model_component_field_count */,
                                     char* name /* [MOIST_FIELD_NAME_MAX + 1] */,
                                     int* dtype /* : MOIST_FIELD_REAL */,
                                     int* rank /* : 0 for a scalar */,
                                     int* dims /* [MOIST_FIELD_MAX_RANK] */,
                                     int* count /* : elements a read writes */) moist_API_SUFFIX__V_1_0;

/// Copy a component result description or query its length, following
/// moist_get_banner
moist_API_ENTRY void moist_API_CALL
moist_get_model_component_field_about(moist_error error, moist_model model,
                                      int component, const char* name, char* about,
                                      size_t capacity, size_t* length) moist_API_SUFFIX__V_1_0;

/// Read a MOIST_FIELD_REAL component result by name: the `count` elements
/// moist_get_model_component_field_info reports. An error for a result the
/// component does not hold or a NULL buffer; nothing is written then
moist_API_ENTRY void moist_API_CALL
moist_get_model_component_field_real(moist_error error,
                                     moist_model model,
                                     int component,
                                     const char* name,
                                     double* values /* [count] */) moist_API_SUFFIX__V_1_0;

/// Copy the settings printout of a continuum model, or query its length,
/// following moist_get_banner
/// The cavity section first, then one section per component headed by its
/// 1-based position: name, description, solvent inputs and registered
/// settings, one "key ... value" line each. Empty at context verbosity 0, the
/// normal settings from 1, the detailed ones from 3. Needs no update. A model
/// of another family is an error; the host chooses where to print
moist_API_ENTRY void moist_API_CALL
moist_get_model_parameters_text(moist_error error, moist_model model,
                                char* buffer, size_t capacity,
                                size_t* length) moist_API_SUFFIX__V_1_0;

/// Copy the radii table of the current structure, or query its length,
/// following moist_get_banner
/// One row per element of the structure: atomic number, symbol and radius in
/// Angstrom and Bohr; for custom radii a short summary of the model. Empty at
/// context verbosity 0, printed from 1. Available right after
/// moist_update_model, also when the cavity waits for the build phase; an
/// error before it. No update prints the table on its own. A model of
/// another family is an error
moist_API_ENTRY void moist_API_CALL
moist_get_model_radii_text(moist_error error, moist_model model,
                           char* buffer, size_t capacity,
                           size_t* length) moist_API_SUFFIX__V_1_0;

/// Gaussian PCM widths and switching factors
/// Call moist_get_cavity_sizes first and pass the capacity allocated for both
/// arrays. Nothing is written when it is smaller than the current size
moist_API_ENTRY void moist_API_CALL
moist_get_cavity_gaussian(moist_error error,
                          moist_cavity cavity,
                          int ngrid_cap /* : allocated grid capacity, at
                                       least the ngrid reported by
                                       moist_get_cavity_sizes */,
                          double* xi /* [ngrid_cap] */,
                          double* f /* [ngrid_cap] */) moist_API_SUFFIX__V_1_0;

/*
 * Gradient API (Tier 3 - Cavity gradients)
**/

/// Compute cavity gradient w.r.t. nuclear coordinates
/// Must be called after moist_update_cavity and before moist_get_cavity_gradient
moist_API_ENTRY void moist_API_CALL
moist_compute_cavity_gradient(moist_error error,
                              moist_cavity cavity) moist_API_SUFFIX__V_1_0;

/// Compute anchor-only nuclear derivatives (callback/isodensity LSF, DROP-specific)
/// Must be called after moist_update_cavity and before the *_rA contractions
/// Restricts each grid point's nuclear coupling to its owner atom's rigid anchor
/// motion (the level set field's nuclear derivatives are zero for callback LSFs)
moist_API_ENTRY void moist_API_CALL
moist_compute_anchor_gradient(moist_error error,
                              moist_cavity cavity) moist_API_SUFFIX__V_1_0;

/// Get anchor-channel nuclear derivatives from the anchor pass
/// Must call moist_compute_anchor_gradient (or moist_compute_cavity_gradient) first
/// Call moist_get_cavity_sizes first and pass the capacities you allocated the
/// arrays with (see the capacity convention at the top of this header)
/// Arrays (row-major order; see the layout convention at the top of this header):
///   xyz1_rA  row-major (ngrid_cap, nsph_cap, 3, 3) - d(r_i)_j / d(R_A)_alpha  (i, A, alpha, j)
///                                       element at [j + 3*(alpha + 3*(A + nsph_cap*i))]
///   xi1_rA   row-major (ngrid_cap, nsph_cap, 3)   - d(xi_i)  / d(R_A)_alpha  (i, A, alpha)
///                                       element at [alpha + 3*(A + nsph_cap*i)]
///   a_i1_rA  row-major (ngrid_cap, nsph_cap, 3)   - d(a_i)   / d(R_A)_alpha  (i, A, alpha)
///   v_i1_rA  row-major (ngrid_cap, nsph_cap, 3)   - d(v_i)   / d(R_A)_alpha  (i, A, alpha)
///   A_tot1_rA row-major (nsph_cap, 3)    - d(total area)   / d(R_A)_alpha  (A, alpha)
///   V_tot1_rA row-major (nsph_cap, 3)    - d(total volume) / d(R_A)_alpha  (A, alpha)
/// a_i1_rA/v_i1_rA are the un-summed per-point counterparts of A_tot1_rA/V_tot1_rA
/// Grid point area carries a switching-function dependence (a_i ~ f_i/xi_i^2):
/// a_i1_rA is not recoverable from xi1_rA alone; it is the area route a generic
/// geometric surface functional (e.g. GOSTSHYP) needs
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
/// the top of this header)
/// Arrays (row-major order; see the layout convention at the top of this header):
///   A_tot1_rA row-major (nsph_cap, 3)    - gradient of total area w.r.t. nuclear coords
///   V_tot1_rA row-major (nsph_cap, 3)    - gradient of total volume w.r.t. nuclear coords
///   asph1_rA  row-major (nsph_cap, nsph_cap, 3) - gradient of per-sphere areas (perturbed atom, owner, xyz)
///   vsph1_rA  row-major (nsph_cap, nsph_cap, 3) - gradient of per-sphere volumes (perturbed atom, owner, xyz)
///   xyz1_rA   row-major (ngrid_cap, nsph_cap, 3, 3) - grid point position derivatives (grid, atom, perturbed_xyz, xyz)
///   r_iI1_rA  row-major (ngrid_cap, nsph_cap, 3) - gradient of grid-owner distances (optional)
///   rho1_rA   row-major (ngrid_cap, nsph_cap, 3) - gradient of rho (anchor-to-surface distance) (optional)
/// The two optional outputs accept NULL. They are computed only for a cavity
/// built with do_fine; a buffer for one that was not computed is an error and
/// nothing is written
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
                          double* r_iI1_rA /* : row-major (ngrid_cap, nsph_cap, 3), or NULL */,
                          double* rho1_rA /* : row-major (ngrid_cap, nsph_cap, 3), or NULL */) moist_API_SUFFIX__V_1_0;

/// Contract DROP surface weights to per-grid LSF adjoint weights (DROP-specific)
/// Includes the projected-coordinate response from w_xyz and the electronic
/// xi path through wleb, cpjac, and the critical-gradient/focus switches
/// Exported w_f is the anchor-only iSwiG overlap and has no electronic
/// LSF response for fixed nuclei
/// Uses cavity-internal ngrid for array extents
/// Arrays (row-major order; see the layout convention at the top of this header):
///   w_xi[ngrid]                        - xi derivative weights
///   w_f[ngrid]                         - anchor switch derivative weights
///   w_xyz     row-major (ngrid, 3)        - coordinate derivative weights
///   w_lsf0[ngrid]                      - LSF value adjoint weights
///   w_lsf1    row-major (ngrid, 3)        - LSF gradient adjoint weights
///   w_lsf2    row-major (ngrid, 3, 3)      - LSF Hessian adjoint weights
/// w_lsf2[i][b][a] is conjugate to the native Hessian(a,b,i)
/// Note: w_lsf2 is NOT symmetric
/// This entry point covers the xi/f/xyz channels only: a thin wrapper
/// forwarding NULL for the extended channels below
moist_API_ENTRY void moist_API_CALL
moist_contract_surface_lsf_weights(moist_error error,
                                   moist_cavity cavity,
                                   const double* w_xi /* [ngrid] */,
                                   const double* w_f /* [ngrid] */,
                                   const double* w_xyz /* : row-major (ngrid, 3) */,
                                   double* w_lsf0 /* [ngrid] */,
                                   double* w_lsf1 /* : row-major (ngrid, 3) */,
                                   double* w_lsf2 /* : row-major (ngrid, 3, 3) */) moist_API_SUFFIX__V_1_0;

/// Extended DROP surface-to-LSF contraction with optional normal/curvature channels
/// Same contract as moist_contract_surface_lsf_weights plus:
///   w_n       row-major (ngrid, 3)        - outward-normal derivative weights (or NULL)
///   w_k1[ngrid]                        - first principal curvature weights (or NULL)
///   w_k2[ngrid]                        - second principal curvature weights (or NULL)
/// Pass NULL for any optional channel to skip it
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

/// Contract direct nuclear and electronic PCM terms
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


/* Ownership: each new_* handle has one owner; constructors copy configuration
 * - Delete couplings before their model
 * - Borrowed model cavities must not outlive that model; deleting their
 *   wrapper does not delete the model's cavity
 * - Responses own their results
 * - Deletion sets the supplied handle to NULL
 */
#ifndef moist_CFFI
#ifdef __cplusplus
inline void moist_delete(moist_context& handle) { moist_delete_context(&handle); }
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
    moist_context: moist_delete_context, \
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
