#include <limits.h>
#include <math.h>
#include <stdatomic.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "moist.h"

/* Fixture builders for numerical tests; public API exercised in test_v1_contract
 * below. Each copies and releases its LSF configuration */
static moist_cavity fixture_drop(moist_error error, const int *nleb, const bool *debug, const int *verbose, const double *blendk, const double *blend1b, const double *blend2b, const double *blend3b, const bool *do_fine, const double *tolerance, const int *proj_maxiter, const int *proj_level, const double *branch_weight_s, const double *rho_grid_h, const int *wleb_prune_level)
{
    moist_drop_options options;
    moist_svdw_options surface;
    moist_init_drop_options(error, &options, sizeof options);
    if (moist_check_error(error)) return NULL;
    moist_init_svdw_options(error, &surface, sizeof surface);
    if (moist_check_error(error)) return NULL;
    if (nleb) options.nleb = *nleb;
    if (debug) options.debug = *debug;
    if (verbose) options.verbosity = *verbose;
    if (do_fine) options.do_fine = *do_fine;
    if (tolerance) options.tolerance = *tolerance;
    if (proj_maxiter) options.proj_maxiter = *proj_maxiter;
    if (proj_level) options.proj_level = *proj_level;
    if (branch_weight_s) options.branch_weight_s = *branch_weight_s;
    if (rho_grid_h) options.rho_grid_h = *rho_grid_h;
    if (wleb_prune_level) options.wleb_prune_level = *wleb_prune_level;
    if (blendk) surface.blend_k = *blendk;
    if (blend1b) surface.blend_1b = *blend1b;
    if (blend2b) surface.blend_2b = *blend2b;
    if (blend3b) surface.blend_3b = *blend3b;
    moist_lsf lsf = moist_new_svdw_lsf(error, &surface);
    if (moist_check_error(error)) return NULL;
    moist_cavity cavity = moist_new_drop_cavity(error, lsf, NULL, &options);
    moist_delete(lsf);
    return cavity;
}

static moist_cavity fixture_drop_with_radii(moist_error error, moist_radii radii, const int *nleb, const bool *debug, const int *verbose, const double *blendk, const double *blend1b, const double *blend2b, const double *blend3b, const bool *do_fine, const double *tolerance, const int *proj_maxiter, const int *proj_level, const double *branch_weight_s, const double *rho_grid_h, const int *wleb_prune_level)
{
    moist_drop_options options;
    moist_svdw_options surface;
    moist_init_drop_options(error, &options, sizeof options);
    if (moist_check_error(error)) return NULL;
    moist_init_svdw_options(error, &surface, sizeof surface);
    if (moist_check_error(error)) return NULL;
    if (nleb) options.nleb = *nleb;
    if (debug) options.debug = *debug;
    if (verbose) options.verbosity = *verbose;
    if (do_fine) options.do_fine = *do_fine;
    if (tolerance) options.tolerance = *tolerance;
    if (proj_maxiter) options.proj_maxiter = *proj_maxiter;
    if (proj_level) options.proj_level = *proj_level;
    if (branch_weight_s) options.branch_weight_s = *branch_weight_s;
    if (rho_grid_h) options.rho_grid_h = *rho_grid_h;
    if (wleb_prune_level) options.wleb_prune_level = *wleb_prune_level;
    if (blendk) surface.blend_k = *blendk;
    if (blend1b) surface.blend_1b = *blend1b;
    if (blend2b) surface.blend_2b = *blend2b;
    if (blend3b) surface.blend_3b = *blend3b;
    moist_lsf lsf = moist_new_svdw_lsf(error, &surface);
    if (moist_check_error(error)) return NULL;
    moist_cavity cavity = moist_new_drop_cavity(error, lsf, radii, &options);
    moist_delete(lsf);
    return cavity;
}

static moist_cavity fixture_iswig(moist_error error, const int *nleb,
    const bool *debug, const int *verbose, const double *cut_a, const double *cut_f)
{
    moist_iswig_options options;
    moist_init_iswig_options(error, &options, sizeof options);
    if (moist_check_error(error)) return NULL;
    if (nleb) options.nleb = *nleb;
    if (debug) options.debug = *debug;
    if (verbose) options.verbosity = *verbose;
    if (cut_a) options.cut_a = *cut_a;
    if (cut_f) options.cut_f = *cut_f;
    return moist_new_iswig_cavity(error, NULL, &options);
}

static moist_cavity fixture_isodensity_internal(moist_error error, int nshell, const int *shell_atom, const int *shell_l, const int *shell_nprim, const double *exps, const double *coeffs, double rho_iso, const double *scale, const int *nleb, const bool *debug, const int *verbose, const bool *do_fine, const int *wleb_prune_level, const double *tolerance)
{
    moist_drop_options options;
    moist_isodensity_options surface;
    moist_init_drop_options(error, &options, sizeof options);
    if (moist_check_error(error)) return NULL;
    moist_init_isodensity_options(error, &surface, sizeof surface);
    if (moist_check_error(error)) return NULL;
    if (nleb) options.nleb = *nleb;
    if (debug) options.debug = *debug;
    if (verbose) options.verbosity = *verbose;
    if (do_fine) options.do_fine = *do_fine;
    if (wleb_prune_level) options.wleb_prune_level = *wleb_prune_level;
    if (tolerance) options.tolerance = *tolerance;
    surface.rho_iso = rho_iso;
    if (scale) surface.scale = *scale;
    moist_lsf lsf = moist_new_isodensity_lsf(error, nshell, shell_atom, shell_l, shell_nprim, exps, coeffs, &surface);
    if (moist_check_error(error)) return NULL;
    moist_cavity cavity = moist_new_drop_cavity(error, lsf, NULL, &options);
    moist_delete(lsf);
    return cavity;
}

static moist_cavity fixture_isodensity_callback(moist_error error, moist_isodensity_lsf_callback callback, void *context, double rho_iso, const double *scale, const int *nleb, const bool *debug, const int *verbose, const bool *do_fine, const int *wleb_prune_level, const double *tolerance)
{
    moist_drop_options options;
    moist_isodensity_options surface;
    moist_init_drop_options(error, &options, sizeof options);
    if (moist_check_error(error)) return NULL;
    moist_init_isodensity_options(error, &surface, sizeof surface);
    if (moist_check_error(error)) return NULL;
    if (nleb) options.nleb = *nleb;
    if (debug) options.debug = *debug;
    if (verbose) options.verbosity = *verbose;
    if (do_fine) options.do_fine = *do_fine;
    if (wleb_prune_level) options.wleb_prune_level = *wleb_prune_level;
    if (tolerance) options.tolerance = *tolerance;
    surface.rho_iso = rho_iso;
    if (scale) surface.scale = *scale;
    moist_lsf lsf = moist_new_isodensity_callback_lsf(error, callback, context, &surface);
    if (moist_check_error(error)) return NULL;
    moist_cavity cavity = moist_new_drop_cavity(error, lsf, NULL, &options);
    moist_delete(lsf);
    return cavity;
}

/* Report a failed handle without terminating the host */
static inline void
show_error(moist_error error)
{
    if (!moist_check_error(error)) return;

    char message[512] = {0};
    const int message_len = (int)sizeof(message);
    moist_get_error(error, message, &message_len);
    printf("[Message] %s\n", message);
}

/* Caller-side abort-on-failure policy; the library itself never terminates the
 * host. Tests below report and return a status instead, so main() runs the
 * rest */
static void
die_on_error(moist_error error, const char* context)
{
    if (!moist_check_error(error)) return;

    char message[512];
    const int message_len = (int)sizeof(message);
    moist_get_error(error, message, &message_len);
    fprintf(stderr, "[moist Error] %s: %s\n", context, message);
    exit(EXIT_FAILURE);
}

/* Normalization coefficient of a single s primitive, (2a/pi)^(3/4); uses
 * acos(-1.0), not M_PI, which is not standard C */
static inline double gto_s_norm(double alpha)
{
    return pow(2.0 * alpha / acos(-1.0), 0.75);
}

/* Read one named cavity field, reporting rather than aborting on failure.
 * moist_get_cavity_field_info returns the element count; the accessor writes
 * exactly that many elements */
static int read_field_real(moist_error error, moist_cavity cav, const char* name,
                           double* values)
{
    moist_get_cavity_field_real(error, cav, name, values);
    if (moist_check_error(error)) {
        printf("  Error: cannot read field '%s'\n", name);
        show_error(error);
        return 1;
    }
    return 0;
}

static int read_field_int(moist_error error, moist_cavity cav, const char* name,
                          int* values)
{
    moist_get_cavity_field_int(error, cav, name, values);
    if (moist_check_error(error)) {
        printf("  Error: cannot read field '%s'\n", name);
        show_error(error);
        return 1;
    }
    return 0;
}

/* DROP fields most callers want, read in one go by name */
static int read_drop_fields(moist_error error, moist_cavity cav, int ngrid,
                            int* nmax, double* normal, double* wleb,
                            double* r_iI0, double* f, double* rho)
{
    return read_field_int(error, cav, "nmax", nmax)
        || read_field_real(error, cav, "normal0", normal)
        || read_field_real(error, cav, "wleb", wleb)
        || read_field_real(error, cav, "r_iI0", r_iI0)
        || read_field_real(error, cav, "f", f)
        || read_field_real(error, cav, "rho", rho);
}

/* Agreement between two independently summed results, relative to the
 * magnitude involved */
static int agrees_to(double a, double b, double rel_tol)
{
    double scale = fabs(a) > fabs(b) ? fabs(a) : fabs(b);
    if (scale < 1.0) {
        scale = 1.0;
    }
    return fabs(a - b) <= rel_tol * scale;
}

/* Native-axis helpers for diagnostic buffers; C axes are the exact reverse */
static inline size_t idx_f2(int i0, int i1, int d0)
{
    return (size_t)i0 + (size_t)d0 * (size_t)i1;
}

static inline size_t idx_f4(int i0, int i1, int i2, int i3, int d0, int d1, int d2)
{
    return (size_t)i0
         + (size_t)d0 * ((size_t)i1
         + (size_t)d1 * ((size_t)i2
         + (size_t)d2 * (size_t)i3));
}

/* H2O geometry most tests below build on, in Bohr. Isodensity callbacks hold a
 * pointer to these centers for the lifetime of a build */
#define H2O_NATOMS 3

static int h2o_numbers[H2O_NATOMS] = {8, 1, 1};
static double h2o_positions[3 * H2O_NATOMS] = {
    0.0000,  0.0000,  0.1173,   /* O */
    0.0000,  1.4309, -0.9370,   /* H */
    0.0000, -1.4309, -0.9370    /* H */
};

static moist_structure make_h2o(moist_error error)
{
    return moist_new_structure(error, H2O_NATOMS, h2o_numbers, h2o_positions,
                               NULL, NULL);
}

int test_version(void)
{
    printf("Start test: version query\n");
    return moist_get_version() > 0 ? 0 : 1;
}

int test_uninitialized_error(void)
{
    printf("Start test: uninitialized error\n");
    moist_error error = moist_new_error();
    int status = moist_check_error(error);
    /* A clean handle must be a no-op for the caller-side abort helper */
    die_on_error(error, "uninitialized error handle");
    moist_delete_error(&error);
    return status == 0 ? 0 : 1;
}

int test_null_handle(void)
{
    printf("Start test: null handle\n");
    return moist_check_error(NULL) == 2 ? 0 : 1;
}

int test_delete_resets_handle(void)
{
    printf("Start test: delete resets handle\n");
    moist_error error = moist_new_error();
    moist_delete_error(&error);
    return error == NULL ? 0 : 1;
}

int test_drop_cavity(void)
{
    printf("Start test: DROP cavity build\n");

    int result = 1;

    moist_error error = moist_new_error();
    moist_radii radii_model = NULL;
    moist_structure mol = NULL;
    moist_cavity cav = NULL;

    /* Heap, not VLAs.
     *
     * Vector arrays are C [ngrid][3]: read component k of point i as
     * xyz[3*i + k] (see the layout convention in moist.h) */
    double* xyz = NULL;
    double* a = NULL;
    int* owner = NULL;
    bool* converged = NULL;
    double* vradii = NULL;
    double* asph = NULL;
    double* normal = NULL;
    double* wleb = NULL;
    double* r_iI0 = NULL;
    double* xi = NULL;
    double* f = NULL;
    double* rho = NULL;
    double* amat0 = NULL;

    const int natoms = 1;
    int numbers[1] = {1};
    double positions[3] = {0.0, 0.0, 0.0};

    mol = moist_new_structure(error, natoms, numbers, positions, NULL, NULL);
    if (moist_check_error(error)) {
        show_error(error);
        goto cleanup;
    }

    // Create explicit radii model and cavity handle (does not build yet)
    radii_model = moist_new_cpcm_radii(error);
    if (moist_check_error(error)) {
        show_error(error);
        goto cleanup;
    }

    cav = fixture_drop_with_radii(error, radii_model,
                                           NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
                                           NULL, NULL, NULL, NULL, NULL, NULL);
    if (moist_check_error(error)) {
        show_error(error);
        goto cleanup;
    }

    // Build the cavity with structure
    moist_update_cavity(error, cav, mol);
    if (moist_check_error(error)) {
        show_error(error);
        goto cleanup;
    }

    // Get generic cavity sizes for array allocation
    double area = 0.0, volume = 0.0;
    int ngrid = 0, nsph = 0, nmax = 0;
    moist_get_cavity_sizes(error, cav, &ngrid, &nsph);
    if (moist_check_error(error)) {
        show_error(error);
        goto cleanup;
    }
    if (ngrid <= 0 || nsph <= 0) {
        printf("  Error: cavity reports empty sizes (ngrid = %d, nsph = %d)\n", ngrid, nsph);
        goto cleanup;
    }

    xyz = (double*)malloc((size_t)3 * ngrid * sizeof(double));
    a = (double*)malloc((size_t)ngrid * sizeof(double));
    owner = (int*)malloc((size_t)ngrid * sizeof(int));
    converged = (bool*)malloc((size_t)ngrid * sizeof(bool));
    vradii = (double*)malloc((size_t)nsph * sizeof(double));
    asph = (double*)malloc((size_t)nsph * sizeof(double));

    // DROP-specific arrays
    normal = (double*)malloc((size_t)3 * ngrid * sizeof(double));
    wleb = (double*)malloc((size_t)ngrid * sizeof(double));
    r_iI0 = (double*)malloc((size_t)ngrid * sizeof(double));
    xi = (double*)malloc((size_t)ngrid * sizeof(double));
    f = (double*)malloc((size_t)ngrid * sizeof(double));
    rho = (double*)malloc((size_t)ngrid * sizeof(double));
    amat0 = (double*)malloc((size_t)ngrid * (size_t)ngrid * sizeof(double));

    if (!xyz || !a || !owner || !converged || !vradii || !asph || !normal ||
        !wleb || !r_iI0 || !xi || !f || !rho || !amat0) {
        printf("Error: Memory allocation failed\n");
        goto cleanup;
    }

    // Get generic cavity results
    moist_get_cavity_results(error, cav, ngrid, nsph,
                             &area, &volume, &ngrid, &nsph,
                             xyz, a, owner, converged, vradii, asph);
    if (moist_check_error(error)) {
        show_error(error);
        goto cleanup;
    }

    // Get DROP results by name
    if (read_drop_fields(error, cav, ngrid, &nmax, normal, wleb, r_iI0, f, rho)) {
        goto cleanup;
    }

    // Assemble A-matrix and get xi values (DROP-specific)
    moist_assemble_amat(error, cav, ngrid, amat0, xi);
    if (moist_check_error(error)) {
        show_error(error);
        goto cleanup;
    }

    // Validate results
    result = area > 0.0 && ngrid > 0 &&
             nmax >= ngrid && wleb[0] > 0.0 && a[0] > 0.0 &&
             r_iI0[0] > 0.0 && xi[0] > 0.0 &&
             f[0] > 0.0 && rho[0] >= 0.0 &&
             converged[0] == true && asph[0] > 0.0 &&
             vradii[0] > 0.0
         ? 0
         : 1;

cleanup:
    free(xyz); free(a); free(owner); free(converged); free(vradii); free(asph);
    free(normal); free(wleb); free(r_iI0); free(xi); free(f); free(rho);
    free(amat0);

    // Free handles and verify they reset to NULL (generic delete)
    moist_delete_cavity(&cav);
    moist_delete_radii(&radii_model);
    moist_delete_structure(&mol);
    moist_delete_error(&error);

    if (cav != NULL || radii_model != NULL || mol != NULL || error != NULL) {
        printf("Error: handles not reset to NULL after deletion\n");
        return 1;
    }

    return result;
}

int test_custom_radii(void)
{
    printf("Start test: custom radii API\n");

    int result = 1;

    moist_error error = moist_new_error();
    moist_radii radii_model = NULL;
    moist_structure mol = NULL;
    moist_cavity cav = NULL;

    double* xyz = NULL;
    double* a = NULL;
    int* owner = NULL;
    bool* converged = NULL;
    double* vradii = NULL;
    double* asph = NULL;

    const int natoms = H2O_NATOMS;

    mol = make_h2o(error);
    if (moist_check_error(error)) {
        show_error(error);
        goto cleanup;
    }

    radii_model = moist_new_custom_radii(error);
    if (moist_check_error(error)) {
        show_error(error);
        goto cleanup;
    }

    // Element-specific custom radii
    int z_list[2] = {1, 8};
    double r_elem[2] = {2.0, 3.2};
    moist_set_custom_radii_elements(error, radii_model, 2, z_list, r_elem);
    if (moist_check_error(error)) {
        show_error(error);
        goto cleanup;
    }

    cav = fixture_drop_with_radii(error, radii_model,
                                           NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
                                           NULL, NULL, NULL, NULL, NULL, NULL);
    if (moist_check_error(error)) {
        show_error(error);
        goto cleanup;
    }

    moist_update_cavity(error, cav, mol);
    if (moist_check_error(error)) {
        show_error(error);
        goto cleanup;
    }

    int ngrid = 0, nsph = 0;
    moist_get_cavity_sizes(error, cav, &ngrid, &nsph);
    if (moist_check_error(error) || nsph != natoms) {
        show_error(error);
        goto cleanup;
    }

    double area = 0.0, volume = 0.0;
    xyz = (double*)malloc((size_t)3 * ngrid * sizeof(double));
    a = (double*)malloc((size_t)ngrid * sizeof(double));
    owner = (int*)malloc((size_t)ngrid * sizeof(int));
    converged = (bool*)malloc((size_t)ngrid * sizeof(bool));
    vradii = (double*)malloc((size_t)nsph * sizeof(double));
    asph = (double*)malloc((size_t)nsph * sizeof(double));
    if (!xyz || !a || !owner || !converged || !vradii || !asph) {
        printf("  Error: memory allocation failed\n");
        goto cleanup;
    }

    moist_get_cavity_results(error, cav, ngrid, nsph, &area, &volume, &ngrid, &nsph,
                             xyz, a, owner, converged, vradii, asph);
    if (moist_check_error(error)) {
        show_error(error);
        goto cleanup;
    }

    int element_ok = fabs(vradii[0] - r_elem[1]) < 1e-12 &&
                     fabs(vradii[1] - r_elem[0]) < 1e-12 &&
                     fabs(vradii[2] - r_elem[0]) < 1e-12;

    /* Grid size can change between rounds; free the first round's buffers before resizing */
    free(xyz); free(a); free(owner); free(converged); free(vradii); free(asph);
    xyz = NULL; a = NULL; owner = NULL; converged = NULL; vradii = NULL; asph = NULL;
    moist_delete_cavity(&cav);

    // Atom-specific custom radii
    double r_atom[3] = {1.7, 2.1, 2.4};
    moist_set_custom_radii_atoms(error, radii_model, natoms, r_atom);
    if (moist_check_error(error)) {
        show_error(error);
        goto cleanup;
    }

    cav = fixture_drop_with_radii(error, radii_model,
                                           NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
                                           NULL, NULL, NULL, NULL, NULL, NULL);
    if (moist_check_error(error)) {
        show_error(error);
        goto cleanup;
    }

    moist_update_cavity(error, cav, mol);
    if (moist_check_error(error)) {
        show_error(error);
        goto cleanup;
    }

    moist_get_cavity_sizes(error, cav, &ngrid, &nsph);
    if (moist_check_error(error) || nsph != natoms) {
        show_error(error);
        goto cleanup;
    }

    xyz = (double*)malloc((size_t)3 * ngrid * sizeof(double));
    a = (double*)malloc((size_t)ngrid * sizeof(double));
    owner = (int*)malloc((size_t)ngrid * sizeof(int));
    converged = (bool*)malloc((size_t)ngrid * sizeof(bool));
    vradii = (double*)malloc((size_t)nsph * sizeof(double));
    asph = (double*)malloc((size_t)nsph * sizeof(double));
    if (!xyz || !a || !owner || !converged || !vradii || !asph) {
        printf("  Error: memory allocation failed\n");
        goto cleanup;
    }

    moist_get_cavity_results(error, cav, ngrid, nsph, &area, &volume, &ngrid, &nsph,
                             xyz, a, owner, converged, vradii, asph);
    if (moist_check_error(error)) {
        show_error(error);
        goto cleanup;
    }

    int atom_ok = fabs(vradii[0] - r_atom[0]) < 1e-12 &&
                  fabs(vradii[1] - r_atom[1]) < 1e-12 &&
                  fabs(vradii[2] - r_atom[2]) < 1e-12;

    result = (element_ok && atom_ok) ? 0 : 1;

cleanup:
    free(xyz); free(a); free(owner); free(converged); free(vradii); free(asph);

    moist_delete_cavity(&cav);
    moist_delete_radii(&radii_model);
    moist_delete_structure(&mol);
    moist_delete_error(&error);
    return result;
}

int test_header_and_version(void)
{
    printf("Start test: header and version printing\n");
    moist_error error = moist_new_error();
    char banner[2048];
    const size_t capacity = sizeof banner;
    size_t length;
    moist_get_banner(error, moist_banner_full, banner, capacity, &length);
    die_on_error(error, "banner");
    fputs(banner, stdout);
    moist_get_version_string(error, banner, capacity, &length);
    die_on_error(error, "version string");
    puts(banner);
    moist_delete(error);
    printf("\n");
    return 0;
}

int test_h2o_cavity(void)
{
    printf("Start test: H2O cavity build\n");

    int result = 1;

    moist_error error = moist_new_error();
    moist_structure mol = NULL;
    moist_cavity cav = NULL;

    // Allocate arrays on the heap, no VLAs
    double* xyz = NULL;
    double* a = NULL;
    int* owner = NULL;
    bool* converged = NULL;
    double* vradii = NULL;
    double* asph = NULL;
    double* normal = NULL;
    double* wleb = NULL;
    double* r_iI0 = NULL;
    double* f = NULL;
    double* gaussian_xi = NULL;
    double* gaussian_f = NULL;
    double* rho = NULL;

    mol = make_h2o(error);
    if (moist_check_error(error)) {
        show_error(error);
        goto cleanup;
    }

    // Create DROP cavity handle (does not build yet)
    cav = fixture_drop(error, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
                                NULL, NULL, NULL, NULL, NULL, NULL);
    if (moist_check_error(error)) {
        show_error(error);
        goto cleanup;
    }

    // Build cavity with structure
    printf("  Building cavity...\n");
    moist_update_cavity(error, cav, mol);
    if (moist_check_error(error)) {
        show_error(error);
        goto cleanup;
    }
    printf("  Cavity built successfully\n");

    // Get cavity sizes for array allocation
    int ngrid = 0, nsph = 0, nmax = 0;
    moist_get_cavity_sizes(error, cav, &ngrid, &nsph);
    if (moist_check_error(error)) {
        show_error(error);
        goto cleanup;
    }

    printf("  Cavity grid points: %d\n", ngrid);
    printf("  Number of spheres: %d\n", nsph);

    xyz = (double*)malloc((size_t)3 * ngrid * sizeof(double));
    a = (double*)malloc((size_t)ngrid * sizeof(double));
    owner = (int*)malloc((size_t)ngrid * sizeof(int));
    converged = (bool*)malloc((size_t)ngrid * sizeof(bool));
    vradii = (double*)malloc((size_t)nsph * sizeof(double));
    asph = (double*)malloc((size_t)nsph * sizeof(double));
    normal = (double*)malloc((size_t)3 * ngrid * sizeof(double));
    wleb = (double*)malloc((size_t)ngrid * sizeof(double));
    r_iI0 = (double*)malloc((size_t)ngrid * sizeof(double));
    f = (double*)malloc((size_t)ngrid * sizeof(double));
    gaussian_xi = (double*)malloc((size_t)ngrid * sizeof(double));
    gaussian_f = (double*)malloc((size_t)ngrid * sizeof(double));
    rho = (double*)malloc((size_t)ngrid * sizeof(double));

    if (!xyz || !a || !owner || !converged || !vradii || !asph || !normal ||
        !wleb || !r_iI0 || !f || !gaussian_xi || !gaussian_f || !rho) {
        printf("Error: Memory allocation failed\n");
        goto cleanup;
    }

    // Get generic cavity results
    double area = 0.0, volume = 0.0;
    moist_get_cavity_results(error, cav, ngrid, nsph, &area, &volume, &ngrid, &nsph,
                             xyz, a, owner, converged, vradii, asph);
    if (moist_check_error(error)) {
        show_error(error);
        goto cleanup;
    }

    printf("  Cavity surface area: %.4f Bohr²\n", area);
    printf("  Cavity volume: %.4f Bohr³\n", volume);

    if (read_drop_fields(error, cav, ngrid, &nmax, normal, wleb, r_iI0, f, rho)) {
        goto cleanup;
    }

    moist_get_cavity_gaussian(error, cav, ngrid, gaussian_xi, gaussian_f);
    if (moist_check_error(error)) {
        show_error(error);
        goto cleanup;
    }

    printf("  Raw grid size (nmax): %d\n", nmax);

    // Validate results
    result = area > 0.0 && volume > 0.0 && ngrid > 0 && nsph == H2O_NATOMS &&
             nmax >= ngrid && wleb[0] > 0.0 && a[0] > 0.0 &&
             gaussian_xi[0] > 0.0 && fabs(gaussian_f[0] - f[0]) < 1e-14 &&
             vradii[0] > 0.0 && vradii[1] > 0.0 && vradii[2] > 0.0
         ? 0
         : 1;

    if (result == 0) {
        printf("  H2O cavity computation completed successfully!\n");
    }

cleanup:
    free(xyz); free(a); free(owner); free(converged); free(vradii); free(asph);
    free(normal); free(wleb); free(r_iI0); free(f);
    free(gaussian_xi); free(gaussian_f); free(rho);

    moist_delete_cavity(&cav);
    moist_delete_structure(&mol);
    moist_delete_error(&error);

    if (cav != NULL || mol != NULL || error != NULL) {
        printf("Error: handles not reset to NULL after deletion\n");
        return 1;
    }

    return result;
}

int test_cavity_gradient(void)
{
    printf("Start test: cavity gradient computation\n");

    /* Single-exit cleanup: every pointer starts NULL here and is freed once
     * under `cleanup:`. `result` starts at failure; only the path that runs
     * to completion clears it */
    int result = 1;

    moist_error error = moist_new_error();
    moist_structure mol = NULL;
    moist_cavity cav = NULL;

    double* A_tot1_rA = NULL;        // (3, nsph)
    double* V_tot1_rA = NULL;        // (3, nsph)
    double* asph1_rA = NULL;         // (3, nsph, nsph)
    double* vsph1_rA = NULL;         // (3, nsph, nsph)
    double* xyz1_rA = NULL;          // (3, 3, nsph, ngrid)
    double* r_iI1_rA = NULL;         // (3, nsph, ngrid)
    double* rho1_rA = NULL;          // (3, nsph, ngrid)
    double* Amat0 = NULL;            // (ngrid, ngrid)
    double* Amat1_rA = NULL;         // (3, nsph, ngrid, ngrid)
    double* xi = NULL;               // (ngrid)
    double* q1 = NULL;               // (ngrid)
    double* q2 = NULL;               // (ngrid)
    double* grad_contract = NULL;    // (3, nsph)
    double* grad_ref = NULL;         // (3, nsph)
    double* surface_q = NULL;        // (ngrid)
    double* surface_q_scaled = NULL; // (ngrid)
    double* qefield = NULL;          // (3, ngrid)
    double* qefield_scaled = NULL;   // (3, ngrid)
    double* za = NULL;               // (nsph)
    double* grad_ne_zero = NULL;     // (3, nsph)
    double* grad_ne = NULL;          // (3, nsph)
    double* grad_ne_scaled = NULL;   // (3, nsph)
    /* Outputs of the surface-to-LSF contraction, own buffers rather than
     * reused from the A-matrix. w_lsf2 is Fortran (3,3,ngrid) */
    double* w_lsf0 = NULL;           // (ngrid)
    double* w_lsf1 = NULL;           // (3, ngrid)
    double* w_lsf2 = NULL;           // (3, 3, ngrid)

    const int natoms = H2O_NATOMS;
    const int* numbers = h2o_numbers;

    mol = make_h2o(error);
    if (moist_check_error(error)) {
        show_error(error);
        goto cleanup;
    }

    // Create DROP cavity
    cav = fixture_drop(error, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
                                NULL, NULL, NULL, NULL, NULL, NULL);
    if (moist_check_error(error)) {
        show_error(error);
        goto cleanup;
    }

    // Build cavity
    moist_update_cavity(error, cav, mol);
    if (moist_check_error(error)) {
        show_error(error);
        goto cleanup;
    }

    // Get cavity sizes
    int ngrid = 0, nsph = 0;
    moist_get_cavity_sizes(error, cav, &ngrid, &nsph);
    if (moist_check_error(error)) {
        show_error(error);
        goto cleanup;
    }

    printf("  Grid: ngrid=%d, nsph=%d\n", ngrid, nsph);

    // Compute cavity gradient w.r.t. nuclear coordinates
    printf("  Computing cavity gradient...\n");
    moist_compute_cavity_gradient(error, cav);
    if (moist_check_error(error)) {
        show_error(error);
        goto cleanup;
    }

    // Allocate gradient arrays
    A_tot1_rA = (double*)malloc(3 * nsph * sizeof(double));
    V_tot1_rA = (double*)malloc(3 * nsph * sizeof(double));
    asph1_rA = (double*)malloc(3 * nsph * nsph * sizeof(double));
    vsph1_rA = (double*)malloc(3 * nsph * nsph * sizeof(double));
    xyz1_rA = (double*)malloc(3 * 3 * nsph * ngrid * sizeof(double));
    r_iI1_rA = (double*)malloc(3 * nsph * ngrid * sizeof(double));
    rho1_rA = (double*)malloc(3 * nsph * ngrid * sizeof(double));

    if (!A_tot1_rA || !V_tot1_rA || !asph1_rA || !vsph1_rA ||
        !xyz1_rA || !r_iI1_rA || !rho1_rA) {
        printf("Error: Memory allocation failed for gradient arrays\n");
        goto cleanup;
    }

    // Get cavity gradient arrays
    moist_get_cavity_gradient(error, cav, nsph, ngrid,
                              A_tot1_rA, V_tot1_rA,
                              asph1_rA, vsph1_rA,
                              xyz1_rA, r_iI1_rA, rho1_rA);
    if (moist_check_error(error)) {
        show_error(error);
        goto cleanup;
    }

    // Print some gradient values
    printf("  Area gradient (dA/dR) for atom 0:\n");
    printf("    x: %12.6f\n", A_tot1_rA[0]);
    printf("    y: %12.6f\n", A_tot1_rA[1]);
    printf("    z: %12.6f\n", A_tot1_rA[2]);

    printf("  Volume gradient (dV/dR) for atom 0:\n");
    printf("    x: %12.6f\n", V_tot1_rA[0]);
    printf("    y: %12.6f\n", V_tot1_rA[1]);
    printf("    z: %12.6f\n", V_tot1_rA[2]);

    // Get A-matrix with gradient
    printf("  Computing A-matrix gradient...\n");
    Amat0 = (double*)malloc(ngrid * ngrid * sizeof(double));
    Amat1_rA = (double*)malloc(3 * nsph * ngrid * ngrid * sizeof(double));
    xi = (double*)malloc(ngrid * sizeof(double));

    if (!Amat0 || !Amat1_rA || !xi) {
        printf("Error: Memory allocation failed for A-matrix arrays\n");
        goto cleanup;
    }

    moist_get_amat_gradient(error, cav, nsph, ngrid, Amat0, Amat1_rA, xi);
    if (moist_check_error(error)) {
        show_error(error);
        goto cleanup;
    }

    printf("  A-matrix diagonal element A[0,0]: %12.6f\n", Amat0[0]);
    printf("  xi[0]: %12.6f\n", xi[0]);

    // Validate Fortran column-major layout and symmetry for Amat0/Amat1_rA
    const double sym_tol = 1e-12;
    int max_check = ngrid < 4 ? ngrid : 4;
    int layout_ok = 1;

    for (int i = 0; i < max_check && layout_ok; i++) {
        for (int j = 0; j < max_check; j++) {
            size_t ij = idx_f2(i, j, ngrid);
            size_t ji = idx_f2(j, i, ngrid);
            if (fabs(Amat0[ij] - Amat0[ji]) > sym_tol) {
                printf("  Error: A-matrix symmetry failed at (%d,%d)\n", i, j);
                layout_ok = 0;
                break;
            }
        }
    }

    if (layout_ok && nsph > 0) {
        for (int i = 0; i < max_check && layout_ok; i++) {
            for (int j = 0; j < max_check; j++) {
                for (int axis = 0; axis < 3; axis++) {
                    size_t ij = idx_f4(axis, 0, i, j, 3, nsph, ngrid);
                    size_t ji = idx_f4(axis, 0, j, i, 3, nsph, ngrid);
                    if (fabs(Amat1_rA[ij] - Amat1_rA[ji]) > sym_tol) {
                        printf("  Error: Amat1_rA symmetry failed at axis=%d i=%d j=%d\n", axis, i, j);
                        layout_ok = 0;
                        break;
                    }
                }
                if (!layout_ok) break;
            }
        }
    }

    // Test API: contract_amat1_q1q2_rA
    q1 = (double*)malloc(ngrid * sizeof(double));
    q2 = (double*)malloc(ngrid * sizeof(double));
    grad_contract = (double*)malloc(3 * nsph * sizeof(double));
    grad_ref = (double*)malloc(3 * nsph * sizeof(double));

    // Test API: contract_pcm_nuclear_gradient
    surface_q = (double*)malloc(ngrid * sizeof(double));
    surface_q_scaled = (double*)malloc(ngrid * sizeof(double));
    qefield = (double*)malloc(3 * ngrid * sizeof(double));
    qefield_scaled = (double*)malloc(3 * ngrid * sizeof(double));
    za = (double*)malloc(nsph * sizeof(double));
    grad_ne_zero = (double*)malloc(3 * nsph * sizeof(double));
    grad_ne = (double*)malloc(3 * nsph * sizeof(double));
    grad_ne_scaled = (double*)malloc(3 * nsph * sizeof(double));

    // Test API: contract_surface_lsf_weights -- one buffer per documented extent
    w_lsf0 = (double*)malloc((size_t)ngrid * sizeof(double));
    w_lsf1 = (double*)malloc((size_t)3 * ngrid * sizeof(double));
    w_lsf2 = (double*)malloc((size_t)9 * ngrid * sizeof(double));

    if (!q1 || !q2 || !grad_contract || !grad_ref ||
        !surface_q || !surface_q_scaled || !qefield || !qefield_scaled ||
        !za || !grad_ne_zero || !grad_ne || !grad_ne_scaled ||
        !w_lsf0 || !w_lsf1 || !w_lsf2) {
        printf("Error: Memory allocation failed for contraction API test arrays\n");
        goto cleanup;
    }

    // Both the original and extended surface-contraction ABI must remain callable
    for (int i = 0; i < ngrid; i++) {
        surface_q[i] = 0.0;
        surface_q_scaled[i] = 0.0;
        for (int iaxis = 0; iaxis < 3; iaxis++) {
            qefield[idx_f2(iaxis, i, 3)] = 0.0;
        }
    }
    moist_contract_surface_lsf_weights(error, cav, surface_q, surface_q_scaled,
                                       qefield, w_lsf0, w_lsf1, w_lsf2);
    moist_contract_surface_lsf_weights_extended(error, cav, surface_q, surface_q_scaled,
                                                qefield, w_lsf0, w_lsf1, w_lsf2,
                                                NULL, NULL, NULL);
    if (moist_check_error(error)) {
        show_error(error);
        goto cleanup;
    }

    for (int i = 0; i < ngrid; i++) {
        q1[i] = xi[i];
        q2[i] = 0.5 * xi[i] + 1.0e-3 * (double)(i + 1);
    }

    moist_contract_amat1_q1q2_rA(error, cav, q1, q2, grad_contract);
    if (moist_check_error(error)) {
        show_error(error);
        goto cleanup;
    }

    // Symmetry check: for symmetric dA, q1^T(dA)q2 == q2^T(dA)q1
    moist_contract_amat1_q1q2_rA(error, cav, q2, q1, grad_ref);
    if (moist_check_error(error)) {
        show_error(error);
        goto cleanup;
    }

    int contract_ok = 1;
    const double contract_tol = 1e-12;
    double contract_norm = 0.0;
    for (int iatom = 0; iatom < nsph && contract_ok; iatom++) {
        for (int iaxis = 0; iaxis < 3; iaxis++) {
            size_t idx = idx_f2(iaxis, iatom, 3);
            if (!agrees_to(grad_contract[idx], grad_ref[idx], contract_tol)) {
                printf("  Error: contract_amat1_q1q2_rA symmetry check failed at axis=%d atom=%d\n", iaxis, iatom);
                contract_ok = 0;
                break;
            }
            if (!isfinite(grad_contract[idx])) {
                printf("  Error: contract_amat1_q1q2_rA produced non-finite values\n");
                contract_ok = 0;
                break;
            }
            contract_norm += grad_contract[idx] * grad_contract[idx];
        }
    }
    /* An all-zero result passes symmetry and finiteness too; check the
     * magnitude as well */
    if (contract_ok && contract_norm < 1e-20) {
        printf("  Error: contract_amat1_q1q2_rA norm is essentially zero\n");
        contract_ok = 0;
    }

    // Zero-input check for contract_pcm_nuclear_gradient: should return zero gradient
    for (int i = 0; i < ngrid; i++) {
        surface_q[i] = 0.0;
        for (int iaxis = 0; iaxis < 3; iaxis++) {
            qefield[idx_f2(iaxis, i, 3)] = 0.0;
        }
    }
    for (int iatom = 0; iatom < nsph; iatom++) {
        za[iatom] = (iatom < natoms) ? (double)numbers[iatom] : 0.0;
    }

    moist_contract_pcm_nuclear_gradient(error, cav, surface_q, qefield, za, grad_ne_zero);
    if (moist_check_error(error)) {
        show_error(error);
        goto cleanup;
    }

    int nuc_elec_ok = 1;
    const double zero_tol = 1e-14;
    for (int iatom = 0; iatom < nsph && nuc_elec_ok; iatom++) {
        for (int iaxis = 0; iaxis < 3; iaxis++) {
            if (fabs(grad_ne_zero[idx_f2(iaxis, iatom, 3)]) > zero_tol) {
                printf("  Error: contract_pcm_nuclear_gradient zero-input check failed\n");
                nuc_elec_ok = 0;
                break;
            }
        }
    }

    // Homogeneity check: f(2*q, 2*qefield) = 2*f(q, qefield)
    for (int i = 0; i < ngrid; i++) {
        surface_q[i] = q1[i];
        surface_q_scaled[i] = 2.0 * surface_q[i];
        for (int iaxis = 0; iaxis < 3; iaxis++) {
            double val = 0.02 * (double)(iaxis + 1) * q2[i];
            qefield[idx_f2(iaxis, i, 3)] = val;
            qefield_scaled[idx_f2(iaxis, i, 3)] = 2.0 * val;
        }
    }

    moist_contract_pcm_nuclear_gradient(error, cav, surface_q, qefield, za, grad_ne);
    if (moist_check_error(error)) {
        show_error(error);
        goto cleanup;
    }

    moist_contract_pcm_nuclear_gradient(error, cav, surface_q_scaled, qefield_scaled, za, grad_ne_scaled);
    if (moist_check_error(error)) {
        show_error(error);
        goto cleanup;
    }

    const double homo_tol = 1e-12;
    for (int iatom = 0; iatom < nsph && nuc_elec_ok; iatom++) {
        for (int iaxis = 0; iaxis < 3; iaxis++) {
            size_t idx = idx_f2(iaxis, iatom, 3);
            if (!isfinite(grad_ne[idx]) || !isfinite(grad_ne_scaled[idx])) {
                printf("  Error: contract_pcm_nuclear_gradient produced non-finite values\n");
                nuc_elec_ok = 0;
                break;
            }
            if (!agrees_to(grad_ne_scaled[idx], 2.0 * grad_ne[idx], homo_tol)) {
                printf("  Error: contract_pcm_nuclear_gradient homogeneity check failed\n");
                nuc_elec_ok = 0;
                break;
            }
        }
    }

    /* A routine that writes nothing but zeros also passes zero-input,
     * finiteness and homogeneity; pin the magnitude too. With nonzero surface
     * charges, a nonzero field and real nuclear charges, the nuclear-electron
     * gradient cannot vanish */
    double nuc_elec_norm = 0.0;
    for (int i = 0; i < 3 * nsph; i++) {
        nuc_elec_norm += grad_ne[i] * grad_ne[i];
    }
    if (nuc_elec_ok && nuc_elec_norm < 1e-20) {
        printf("  Error: contract_pcm_nuclear_gradient gradient is essentially zero\n");
        nuc_elec_ok = 0;
    }

    // Validate results: gradients computed and non-trivial
    result = 0;
    // A_tot1_rA should have non-zero values (area depends on nuclear coords)
    double grad_sum = 0.0;
    for (int i = 0; i < 3 * nsph; i++) {
        grad_sum += A_tot1_rA[i] * A_tot1_rA[i];
    }
    if (grad_sum < 1e-20) {
        printf("  Error: Area gradient is essentially zero\n");
        result = 1;
    }
    if (!layout_ok) {
        result = 1;
    }
    if (!contract_ok || !nuc_elec_ok) {
        result = 1;
    }

    if (result == 0) {
        printf("  Gradient computation completed successfully!\n");
    }

cleanup:
    free(A_tot1_rA); free(V_tot1_rA); free(asph1_rA); free(vsph1_rA);
    free(xyz1_rA); free(r_iI1_rA); free(rho1_rA);
    free(Amat0); free(Amat1_rA); free(xi);
    free(q1); free(q2); free(grad_contract); free(grad_ref);
    free(surface_q); free(surface_q_scaled); free(qefield); free(qefield_scaled);
    free(za); free(grad_ne_zero); free(grad_ne); free(grad_ne_scaled);
    free(w_lsf0); free(w_lsf1); free(w_lsf2);

    moist_delete_cavity(&cav);
    moist_delete_structure(&mol);
    moist_delete_error(&error);
    return result;
}

/* Internal isodensity DROP cavity: basis in, density in, cavity out.
 *
 * One-primitive s shell per atom; level set S = scale * (rho_iso - rho) is an
 * analytic sum of spherical Gaussians. Exercises the full internal-backend
 * chain: constructor -> two-pass layout query -> density install -> build ->
 * numbering / anchor gradients / tolerance getter.
 */
int test_isodensity_internal_cavity(void)
{
    printf("Start test: internal isodensity cavity\n");

    int result = 1;

    moist_error error = moist_new_error();
    moist_structure mol = NULL;
    moist_cavity cav = NULL;

    int* shell_off = NULL;
    int* comp_lx = NULL;
    int* comp_ly = NULL;
    int* comp_lz = NULL;
    double* dcart = NULL;
    int* numbering = NULL;

    const int natoms = H2O_NATOMS;

    /* One normalized s primitive per atom */
    const int nshell_in = 3;
    int shell_atom[3] = {0, 1, 2};   /* 0-based atom indices */
    int shell_l[3] = {0, 0, 0};
    int shell_nprim[3] = {1, 1, 1};
    double exps[3];
    double coeffs[3];
    const double alpha = 0.3;
    for (int i = 0; i < nshell_in; i++) {
        exps[i] = alpha;
        coeffs[i] = gto_s_norm(alpha);
    }
    const double rho_iso = 1.0e-3;

    mol = make_h2o(error);
    if (moist_check_error(error)) {
        show_error(error);
        goto cleanup;
    }

    cav = fixture_isodensity_internal(
        error, nshell_in, shell_atom, shell_l, shell_nprim, exps, coeffs,
        rho_iso, NULL, NULL, NULL, NULL, NULL, NULL, NULL);
    if (moist_check_error(error)) {
        show_error(error);
        goto cleanup;
    }

    /* Master tolerance getter (default value, must be positive) */
    double tolerance = 0.0;
    moist_get_drop_cavity_tolerance(error, cav, &tolerance);
    if (moist_check_error(error) || !(tolerance > 0.0)) {
        show_error(error);
        printf("  Error: non-positive DROP tolerance (%.3e)\n", tolerance);
        goto cleanup;
    }
    printf("  DROP tolerance: %.3e\n", tolerance);

    /* First layout pass: sizes only */
    int ncart = 0, nshell_out = 0;
    moist_get_isodensity_cart_layout(error, cav, &ncart, &nshell_out,
                                     NULL, NULL, NULL, NULL);
    if (moist_check_error(error)) {
        show_error(error);
        goto cleanup;
    }
    printf("  Cartesian layout: ncart = %d, nshell = %d\n", ncart, nshell_out);

    int layout_ok = (ncart == nshell_in) && (nshell_out == nshell_in);

    /* Second layout pass: offsets and monomial powers */
    shell_off = (int*)malloc((size_t)(nshell_out + 1) * sizeof(int));
    comp_lx = (int*)malloc((size_t)ncart * sizeof(int));
    comp_ly = (int*)malloc((size_t)ncart * sizeof(int));
    comp_lz = (int*)malloc((size_t)ncart * sizeof(int));
    if (!shell_off || !comp_lx || !comp_ly || !comp_lz) {
        printf("  Error: memory allocation failed\n");
        goto cleanup;
    }

    moist_get_isodensity_cart_layout(error, cav, &ncart, &nshell_out,
                                     shell_off, comp_lx, comp_ly, comp_lz);
    if (moist_check_error(error)) {
        show_error(error);
        goto cleanup;
    }
    for (int s = 0; s <= nshell_out; s++) {
        if (shell_off[s] != s) layout_ok = 0;   /* one s component per shell */
    }
    for (int c = 0; c < ncart; c++) {
        if (comp_lx[c] != 0 || comp_ly[c] != 0 || comp_lz[c] != 0) layout_ok = 0;
    }
    if (!layout_ok) {
        printf("  Error: unexpected cartesian layout for a pure-s basis\n");
    }

    /* A wrong-size density must be rejected before anything is installed.
     * Uses a scratch handle; the expected failure does not poison `error` */
    double reject = 0.0;
    moist_error probe = moist_new_error();
    moist_set_isodensity_density(probe, cav, ncart + 1, &reject);
    int reject_ok = moist_check_error(probe) != 0;
    moist_delete_error(&probe);
    if (!reject_ok) {
        printf("  Error: mismatched density dimension was not rejected\n");
    }

    /* Diagonal density: rho(r) = sum_A 2 * coeff^2 * exp(-2 alpha |r - R_A|^2) */
    dcart = (double*)calloc((size_t)ncart * (size_t)ncart, sizeof(double));
    if (!dcart) {
        printf("  Error: memory allocation failed\n");
        goto cleanup;
    }
    for (int c = 0; c < ncart; c++) {
        dcart[idx_f2(c, c, ncart)] = 2.0;
    }
    moist_set_isodensity_density(error, cav, ncart, dcart);
    if (moist_check_error(error)) {
        show_error(error);
        goto cleanup;
    }

    printf("  Building isodensity cavity...\n");
    moist_update_cavity(error, cav, mol);
    if (moist_check_error(error)) {
        show_error(error);
        goto cleanup;
    }

    int ngrid = 0, nsph = 0;
    moist_get_cavity_sizes(error, cav, &ngrid, &nsph);
    if (moist_check_error(error)) {
        show_error(error);
        goto cleanup;
    }
    printf("  ngrid = %d, nsph = %d\n", ngrid, nsph);

    result = (layout_ok && reject_ok) ? 0 : 1;
    if (ngrid <= 0 || nsph != natoms) {
        printf("  Error: unexpected cavity sizes\n");
        result = 1;
        goto cleanup;
    }

    /* Stable per-point numbering must be available after a build */
    numbering = (int*)malloc((size_t)ngrid * sizeof(int));
    if (!numbering) {
        printf("  Error: memory allocation failed\n");
        result = 1;
        goto cleanup;
    }
    if (read_field_int(error, cav, "numbering", numbering)) {
        result = 1;
    } else {
        int numbering_ok = 1;
        for (int i = 0; i < ngrid; i++) {
            if (numbering[i] <= 0) numbering_ok = 0;
        }
        if (!numbering_ok) {
            printf("  Error: numbering contains non-positive ids\n");
            result = 1;
        }
    }

    /* Anchor-only nuclear derivatives (the isodensity gradient route) */
    moist_compute_anchor_gradient(error, cav);
    if (moist_check_error(error)) {
        show_error(error);
        result = 1;
    } else {
        double* xyz1_rA = (double*)malloc((size_t)9 * nsph * ngrid * sizeof(double));
        double* xi1_rA = (double*)malloc((size_t)3 * nsph * ngrid * sizeof(double));
        double* a_i1_rA = (double*)malloc((size_t)3 * nsph * ngrid * sizeof(double));
        double* v_i1_rA = (double*)malloc((size_t)3 * nsph * ngrid * sizeof(double));
        double* A_tot1_rA = (double*)malloc((size_t)3 * nsph * sizeof(double));
        double* V_tot1_rA = (double*)malloc((size_t)3 * nsph * sizeof(double));

        if (!xyz1_rA || !xi1_rA || !a_i1_rA || !v_i1_rA || !A_tot1_rA || !V_tot1_rA) {
            printf("  Error: memory allocation failed\n");
            result = 1;
        } else {
            moist_get_anchor_gradient(error, cav, nsph, ngrid, xyz1_rA, xi1_rA,
                                      a_i1_rA, v_i1_rA, A_tot1_rA, V_tot1_rA);
            if (moist_check_error(error)) {
                show_error(error);
                result = 1;
            } else {
                /* The un-summed area elements must add up to the total-area gradient */
                double max_dev = 0.0;
                for (int k = 0; k < 3 * nsph; k++) {
                    double acc = 0.0;
                    for (int i = 0; i < ngrid; i++) {
                        acc += a_i1_rA[(size_t)k + (size_t)3 * nsph * (size_t)i];
                    }
                    double dev = fabs(acc - A_tot1_rA[k]);
                    if (dev > max_dev) max_dev = dev;
                }
                printf("  max |sum_i a_i1_rA - A_tot1_rA| = %.3e\n", max_dev);
                if (!(max_dev < 1.0e-8)) {
                    printf("  Error: per-point area derivatives do not sum to the total\n");
                    result = 1;
                }
            }
        }

        free(xyz1_rA); free(xi1_rA); free(a_i1_rA); free(v_i1_rA);
        free(A_tot1_rA); free(V_tot1_rA);
    }

    if (result == 0) {
        printf("  Internal isodensity cavity completed successfully!\n");
    }

cleanup:
    free(numbering);
    free(dcart); free(shell_off); free(comp_lx); free(comp_ly); free(comp_lz);

    moist_delete_cavity(&cav);
    moist_delete_structure(&mol);
    moist_delete_error(&error);
    return result;
}

/* ---------------------------------------------------------------------------
 * Isodensity level set supplied through the C callback ABI
 *
 * Reference implementation of the V_1_0 callback contract: `hess` and `third`
 * may be NULL, meaning "do not compute this derivative", not merely "do not
 * store it". Check for NULL before dereferencing them.
 *
 * Mirrors the diagonal-density case built by test_isodensity_internal_cavity();
 * the two backends must agree:
 *   rho(r) = sum_A c exp(-a |r - R_A|^2)
 *
 * Returns that density and its spatial derivatives; moist forms
 * S(r) = scale * (rho_iso - rho(r)) from the rho_iso passed to the factory,
 * with c = 2 coeff^2 and a = 2 alpha.
 * ------------------------------------------------------------------------- */

/* moist_get_cavity_results writes every array unconditionally; the caller
 * must size and allocate them even when only the scalars are wanted */
static int cavity_area_volume(moist_error error, moist_cavity cav,
                              double* area, double* volume, int* ngrid, int* nsph)
{
    moist_get_cavity_sizes(error, cav, ngrid, nsph);
    if (moist_check_error(error)) return 1;
    /* This failure path does not set the error handle; report what went
     * wrong here */
    if (*ngrid <= 0 || *nsph <= 0) {
        printf("  Error: cavity reports empty sizes (ngrid = %d, nsph = %d)\n",
               *ngrid, *nsph);
        return 1;
    }

    double* xyz = (double*)malloc((size_t)(3 * *ngrid) * sizeof(double));
    double* a = (double*)malloc((size_t)*ngrid * sizeof(double));
    int* owner = (int*)malloc((size_t)*ngrid * sizeof(int));
    bool* converged = (bool*)malloc((size_t)*ngrid * sizeof(bool));
    double* radii = (double*)malloc((size_t)*nsph * sizeof(double));
    double* asph = (double*)malloc((size_t)*nsph * sizeof(double));

    int stat = 1;
    if (!xyz || !a || !owner || !converged || !radii || !asph) {
        printf("  Error: memory allocation failed\n");
    } else {
        moist_get_cavity_results(error, cav, *ngrid, *nsph, area, volume, ngrid, nsph,
                                 xyz, a, owner, converged, radii, asph);
        stat = moist_check_error(error) ? 1 : 0;
    }

    free(xyz); free(a); free(owner); free(converged); free(radii); free(asph);
    return stat;
}

struct iso_callback_ctx {
    int natoms;
    const double* centers;   /* flat, Fortran (3,natoms): centers[3*A + k], Bohr */
    double c;                /* Gaussian prefactor */
    double a;                /* Gaussian exponent  */
    double rho_iso;
    /* Contract instrumentation: how often each derivative order was requested.
     * The projection evaluates grid points in an OpenMP loop; every counter
     * here is shared across threads */
    atomic_int calls;
    atomic_int calls_with_hess;
    atomic_int calls_with_third;
};

static int iso_gaussian_callback(void* context, const double* point,
                                 double* rho_out, double* drho_out,
                                 double* d2rho_out, double* d3rho_out)
{
    struct iso_callback_ctx* ctx = (struct iso_callback_ctx*)context;
    const int want_hess = (d2rho_out != NULL);
    const int want_third = (d3rho_out != NULL);

    atomic_fetch_add(&ctx->calls, 1);
    if (want_hess) atomic_fetch_add(&ctx->calls_with_hess, 1);
    if (want_third) atomic_fetch_add(&ctx->calls_with_third, 1);

    double rho = 0.0;
    double drho[3] = {0.0, 0.0, 0.0};
    double d2rho[9];
    double d3rho[27];
    if (want_hess) for (int i = 0; i < 9; i++) d2rho[i] = 0.0;
    if (want_third) for (int i = 0; i < 27; i++) d3rho[i] = 0.0;

    for (int A = 0; A < ctx->natoms; A++) {
        double d[3];
        double r2 = 0.0;
        for (int i = 0; i < 3; i++) {
            d[i] = point[i] - ctx->centers[3 * A + i];
            r2 += d[i] * d[i];
        }
        const double g = ctx->c * exp(-ctx->a * r2);
        rho += g;
        for (int i = 0; i < 3; i++) {
            drho[i] += -2.0 * ctx->a * d[i] * g;
        }
        if (want_hess) {
            for (int i = 0; i < 3; i++) {
                for (int j = 0; j < 3; j++) {
                    const double dij = (i == j) ? 1.0 : 0.0;
                    d2rho[3 * j + i] +=
                        (4.0 * ctx->a * ctx->a * d[i] * d[j] - 2.0 * ctx->a * dij) * g;
                }
            }
        }
        if (want_third) {
            for (int i = 0; i < 3; i++) {
                for (int j = 0; j < 3; j++) {
                    for (int k = 0; k < 3; k++) {
                        const double dij = (i == j) ? 1.0 : 0.0;
                        const double dik = (i == k) ? 1.0 : 0.0;
                        const double djk = (j == k) ? 1.0 : 0.0;
                        d3rho[9 * k + 3 * j + i] +=
                            (-8.0 * ctx->a * ctx->a * ctx->a * d[i] * d[j] * d[k]
                             + 4.0 * ctx->a * ctx->a * (d[i] * djk + d[j] * dik + d[k] * dij)) * g;
                    }
                }
            }
        }
    }

    /* Bare density; moist subtracts ctx->rho_iso (handed at construction) and
     * applies the DROP sign convention itself */
    *rho_out = rho;
    for (int i = 0; i < 3; i++) drho_out[i] = drho[i];
    if (want_hess) for (int i = 0; i < 9; i++) d2rho_out[i] = d2rho[i];
    if (want_third) for (int i = 0; i < 27; i++) d3rho_out[i] = d3rho[i];

    return 0;   /* success; any nonzero value would abort the cavity build */
}

int test_isodensity_callback_cavity(void)
{
    printf("Start test: callback isodensity cavity\n");
#ifndef moist_API_SUFFIX__V_1_0
#  error "this example implements the V_1_0 isodensity callback contract"
#endif

    int result = 1;

    moist_error error = moist_new_error();
    moist_structure mol = NULL;
    moist_cavity cav = NULL;
    moist_cavity ref = NULL;
    double* dcart = NULL;

    const int natoms = H2O_NATOMS;

    const double alpha = 0.3;
    const double coeff = gto_s_norm(alpha);
    const double rho_iso = 1.0e-3;

    struct iso_callback_ctx ctx;
    ctx.natoms = natoms;
    ctx.centers = h2o_positions;
    ctx.c = 2.0 * coeff * coeff;   /* matches the diagonal density D_cc = 2 */
    ctx.a = 2.0 * alpha;
    ctx.rho_iso = rho_iso;
    atomic_init(&ctx.calls, 0);
    atomic_init(&ctx.calls_with_hess, 0);
    atomic_init(&ctx.calls_with_third, 0);

    mol = make_h2o(error);
    if (moist_check_error(error)) {
        show_error(error);
        goto cleanup;
    }

    cav = fixture_isodensity_callback(
        error, iso_gaussian_callback, &ctx, ctx.rho_iso,
        NULL, NULL, NULL, NULL, NULL, NULL, NULL);
    if (moist_check_error(error)) {
        show_error(error);
        goto cleanup;
    }

    moist_update_cavity(error, cav, mol);
    if (moist_check_error(error)) {
        show_error(error);
        goto cleanup;
    }

    double area = 0.0, volume = 0.0;
    int ngrid = 0, nsph = 0;
    if (cavity_area_volume(error, cav, &area, &volume, &ngrid, &nsph) != 0) {
        show_error(error);
        goto cleanup;
    }
    printf("  ngrid = %d, nsph = %d, area = %.10f, volume = %.10f\n",
           ngrid, nsph, area, volume);

    /* Snapshot the counters once, after the build; the report and the
     * assertions below use the same numbers */
    const int n_calls = atomic_load(&ctx.calls);
    const int n_hess = atomic_load(&ctx.calls_with_hess);
    const int n_third = atomic_load(&ctx.calls_with_third);
    printf("  callback invocations: %d total, %d with Hessian, %d with third derivative\n",
           n_calls, n_hess, n_third);

    result = 0;
    if (ngrid <= 0 || nsph != natoms) {
        printf("  FAIL: unexpected cavity sizes\n");
        result = 1;
    }

    /* The contract is only exercised if both phases occurred: some calls
     * skipped the higher derivatives, some did not */
    if (n_hess <= 0 || n_hess >= n_calls) {
        printf("  FAIL: the value+gradient-only phase was never exercised\n");
        result = 1;
    }

    /* Third-derivative channel only requested at max_deriv >= 3, which a plain
     * cavity build never reaches (see test_isodensity_callback_third_derivative()).
     * Assert the nesting: third derivatives are never requested without a
     * Hessian alongside them */
    if (n_third > n_hess) {
        printf("  FAIL: third derivatives were requested without a Hessian\n");
        result = 1;
    }

    /* Cross-check against the internal backend built from the equivalent
     * basis: same level set, same surface */
    const int nshell_in = 3;
    int shell_atom[3] = {0, 1, 2};
    int shell_l[3] = {0, 0, 0};
    int shell_nprim[3] = {1, 1, 1};
    double exps[3], coeffs[3];
    for (int i = 0; i < nshell_in; i++) {
        exps[i] = alpha;
        coeffs[i] = coeff;
    }

    ref = fixture_isodensity_internal(
        error, nshell_in, shell_atom, shell_l, shell_nprim, exps, coeffs,
        rho_iso, NULL, NULL, NULL, NULL, NULL, NULL, NULL);
    if (moist_check_error(error)) {
        show_error(error);
        result = 1;
        goto cleanup;
    }

    int ncart = 0, nshell_out = 0;
    moist_get_isodensity_cart_layout(error, ref, &ncart, &nshell_out, NULL, NULL, NULL, NULL);
    if (moist_check_error(error)) {
        show_error(error);
        result = 1;
        goto cleanup;
    }

    dcart = (double*)calloc((size_t)ncart * (size_t)ncart, sizeof(double));
    if (!dcart) {
        printf("  Error: memory allocation failed\n");
        result = 1;
        goto cleanup;
    }
    for (int c = 0; c < ncart; c++) {
        dcart[idx_f2(c, c, ncart)] = 2.0;
    }
    moist_set_isodensity_density(error, ref, ncart, dcart);
    moist_update_cavity(error, ref, mol);
    if (moist_check_error(error)) {
        show_error(error);
        result = 1;
        goto cleanup;
    }

    double ref_area = 0.0, ref_volume = 0.0;
    int ref_ngrid = 0, ref_nsph = 0;
    if (cavity_area_volume(error, ref, &ref_area, &ref_volume, &ref_ngrid, &ref_nsph) != 0) {
        show_error(error);
        result = 1;
        goto cleanup;
    }

    const double area_diff = fabs(area - ref_area);
    const double volume_diff = fabs(volume - ref_volume);
    printf("  vs internal backend: |dA| = %.3e, |dV| = %.3e\n", area_diff, volume_diff);
    if (ngrid != ref_ngrid || area_diff > 1.0e-8 || volume_diff > 1.0e-8) {
        printf("  FAIL: callback and internal backends disagree\n");
        result = 1;
    }

    if (result == 0) {
        printf("  Callback isodensity cavity completed successfully!\n");
    }

cleanup:
    free(dcart);
    moist_delete_cavity(&ref);
    moist_delete_cavity(&cav);
    moist_delete_structure(&mol);
    moist_delete_error(&error);
    return result;
}

/* Third-derivative channel of the V_1_0 callback contract.
 *
 * A plain cavity build never reaches it: the projection runs the level-set
 * model at max_deriv 1 and 2, so `third` is always NULL there and the d3rho
 * branch of iso_gaussian_callback never executes.
 *
 * The surface-to-LSF contraction raises the model to max_deriv = 3 (see
 * src/moist/cavity/drop/derivatives/potential.f90); running it against a
 * callback-backed cavity exercises that branch */
int test_isodensity_callback_third_derivative(void)
{
    printf("Start test: callback isodensity third-derivative channel\n");

    int result = 1;

    moist_error error = moist_new_error();
    moist_structure mol = NULL;
    moist_cavity cav = NULL;

    double* w_xi = NULL;     // (ngrid)
    double* w_f = NULL;      // (ngrid)
    double* w_xyz = NULL;    // (3, ngrid)
    double* w_lsf0 = NULL;   // (ngrid)
    double* w_lsf1 = NULL;   // (3, ngrid)
    double* w_lsf2 = NULL;   // (3, 3, ngrid)

    const int natoms = H2O_NATOMS;

    const double alpha = 0.3;
    const double coeff = gto_s_norm(alpha);

    struct iso_callback_ctx ctx;
    ctx.natoms = natoms;
    ctx.centers = h2o_positions;
    ctx.c = 2.0 * coeff * coeff;
    ctx.a = 2.0 * alpha;
    ctx.rho_iso = 1.0e-3;
    atomic_init(&ctx.calls, 0);
    atomic_init(&ctx.calls_with_hess, 0);
    atomic_init(&ctx.calls_with_third, 0);

    mol = make_h2o(error);
    if (moist_check_error(error)) {
        show_error(error);
        goto cleanup;
    }

    cav = fixture_isodensity_callback(
        error, iso_gaussian_callback, &ctx, ctx.rho_iso,
        NULL, NULL, NULL, NULL, NULL, NULL, NULL);
    if (moist_check_error(error)) {
        show_error(error);
        goto cleanup;
    }

    moist_update_cavity(error, cav, mol);
    if (moist_check_error(error)) {
        show_error(error);
        goto cleanup;
    }

    moist_compute_cavity_gradient(error, cav);
    if (moist_check_error(error)) {
        show_error(error);
        goto cleanup;
    }

    int ngrid = 0, nsph = 0;
    moist_get_cavity_sizes(error, cav, &ngrid, &nsph);
    if (moist_check_error(error) || ngrid <= 0 || nsph != natoms) {
        show_error(error);
        printf("  FAIL: unexpected cavity sizes (ngrid = %d, nsph = %d)\n", ngrid, nsph);
        goto cleanup;
    }

    w_xi = (double*)malloc((size_t)ngrid * sizeof(double));
    w_f = (double*)malloc((size_t)ngrid * sizeof(double));
    w_xyz = (double*)malloc((size_t)3 * ngrid * sizeof(double));
    w_lsf0 = (double*)malloc((size_t)ngrid * sizeof(double));
    w_lsf1 = (double*)malloc((size_t)3 * ngrid * sizeof(double));
    w_lsf2 = (double*)malloc((size_t)9 * ngrid * sizeof(double));

    if (!w_xi || !w_f || !w_xyz || !w_lsf0 || !w_lsf1 || !w_lsf2) {
        printf("  Error: memory allocation failed\n");
        goto cleanup;
    }

    /* Nonzero surface adjoints; zero values would let the contraction
     * short-circuit before touching the level-set model */
    for (int i = 0; i < ngrid; i++) {
        w_xi[i] = 1.0e-3 * (double)(i + 1);
        w_f[i] = 2.0e-3;
        for (int iaxis = 0; iaxis < 3; iaxis++) {
            w_xyz[idx_f2(iaxis, i, 3)] = 1.0e-4 * (double)(iaxis + 1);
        }
    }

    /* Count only what the contraction asks for; the build above already ran
     * the value/gradient and Hessian phases */
    atomic_store(&ctx.calls, 0);
    atomic_store(&ctx.calls_with_hess, 0);
    atomic_store(&ctx.calls_with_third, 0);

    moist_contract_surface_lsf_weights(error, cav, w_xi, w_f, w_xyz,
                                       w_lsf0, w_lsf1, w_lsf2);
    if (moist_check_error(error)) {
        show_error(error);
        goto cleanup;
    }

    const int n_calls = atomic_load(&ctx.calls);
    const int n_hess = atomic_load(&ctx.calls_with_hess);
    const int n_third = atomic_load(&ctx.calls_with_third);
    printf("  contraction callbacks: %d total, %d with Hessian, %d with third derivative\n",
           n_calls, n_hess, n_third);

    result = 0;

    if (n_third <= 0) {
        printf("  FAIL: the third-derivative channel was never requested\n");
        result = 1;
    }

    /* max_deriv = 3 implies max_deriv >= 2: every call that got a third
     * derivative also got a Hessian */
    if (n_third > n_hess) {
        printf("  FAIL: third derivatives were requested without a Hessian\n");
        result = 1;
    }

    /* The weights the third-derivative channel feeds must be usable numbers */
    double w_norm = 0.0;
    int finite_ok = 1;
    for (int i = 0; i < ngrid && finite_ok; i++) {
        if (!isfinite(w_lsf0[i])) finite_ok = 0;
        w_norm += w_lsf0[i] * w_lsf0[i];
    }
    for (int i = 0; i < 3 * ngrid && finite_ok; i++) {
        if (!isfinite(w_lsf1[i])) finite_ok = 0;
        w_norm += w_lsf1[i] * w_lsf1[i];
    }
    for (int i = 0; i < 9 * ngrid && finite_ok; i++) {
        if (!isfinite(w_lsf2[i])) finite_ok = 0;
        w_norm += w_lsf2[i] * w_lsf2[i];
    }
    if (!finite_ok) {
        printf("  FAIL: the contraction produced non-finite LSF weights\n");
        result = 1;
    } else if (w_norm < 1e-30) {
        printf("  FAIL: the contraction produced an all-zero adjoint\n");
        result = 1;
    }

    if (result == 0) {
        printf("  Third-derivative callback channel exercised!\n");
    }

cleanup:
    free(w_xi); free(w_f); free(w_xyz);
    free(w_lsf0); free(w_lsf1); free(w_lsf2);

    moist_delete_cavity(&cav);
    moist_delete_structure(&mol);
    moist_delete_error(&error);
    return result;
}
/* A callback that cannot evaluate must be able to say so.
 *
 * Grid points run in an OpenMP loop; failing on the 51st evaluation aborts
 * from the middle of the loop, exercising the shared failure flag */
struct iso_fail_ctx {
    struct iso_callback_ctx base;
    /* Number of successful evaluations before the callback starts failing */
    int fail_after;
    /* Status the callback reports once it starts failing */
    int status;
    /* Evaluation counter, shared across the projection's OpenMP threads */
    atomic_int calls;
};

static int iso_failing_callback(void* context, const double* point,
                                double* rho_out, double* drho_out,
                                double* d2rho_out, double* d3rho_out)
{
    struct iso_fail_ctx* ctx = (struct iso_fail_ctx*)context;

    const int n = atomic_fetch_add(&ctx->calls, 1) + 1;
    if (n > ctx->fail_after) return ctx->status;

    return iso_gaussian_callback(&ctx->base, point, rho_out, drho_out,
                                 d2rho_out, d3rho_out);
}

int test_isodensity_callback_failure(void)
{
    printf("Start test: callback isodensity failure aborts the build\n");

    int result = 1;

    moist_error error = moist_new_error();
    moist_structure mol = NULL;
    moist_cavity cav = NULL;

    const int natoms = H2O_NATOMS;

    const double alpha = 0.3;
    const double coeff = gto_s_norm(alpha);

    struct iso_fail_ctx ctx;
    ctx.base.natoms = natoms;
    ctx.base.centers = h2o_positions;
    ctx.base.c = 2.0 * coeff * coeff;
    ctx.base.a = 2.0 * alpha;
    ctx.base.rho_iso = 1.0e-3;
    atomic_init(&ctx.base.calls, 0);
    atomic_init(&ctx.base.calls_with_hess, 0);
    atomic_init(&ctx.base.calls_with_third, 0);
    ctx.fail_after = 50;
    ctx.status = 7;
    atomic_init(&ctx.calls, 0);

    mol = make_h2o(error);
    if (moist_check_error(error)) {
        show_error(error);
        goto cleanup;
    }

    cav = fixture_isodensity_callback(
        error, iso_failing_callback, &ctx, ctx.base.rho_iso,
        NULL, NULL, NULL, NULL, NULL, NULL, NULL);
    if (moist_check_error(error)) {
        show_error(error);
        goto cleanup;
    }

    result = 0;

    moist_update_cavity(error, cav, mol);
    if (!moist_check_error(error)) {
        printf("  FAIL: a failing callback still produced a cavity\n");
        result = 1;
    } else {
        char message[512] = {0};
        const int message_len = (int)sizeof(message);
        moist_get_error(error, message, &message_len);
        printf("  reported: %s\n", message);
        if (strstr(message, "External LSF evaluation failed") == NULL ||
            strstr(message, "status 7") == NULL) {
            printf("  FAIL: the error does not identify the callback status\n");
            result = 1;
        }
        /* Reading a message does not clear the handle; swap in a fresh one
         * before the recovery build below */
        moist_delete_error(&error);
        error = moist_new_error();
    }

    /* The callback must not be re-entered indefinitely: once it reports
     * failure, moist stops calling it, and the count stays close to the point
     * of failure, not the whole grid */
    const int total_calls = atomic_load(&ctx.calls);
    printf("  callback invocations: %d (failure requested after %d)\n",
           total_calls, ctx.fail_after);
    if (total_calls <= ctx.fail_after) {
        printf("  FAIL: the failing branch was never reached\n");
        result = 1;
    }

    /* The build must also be repeatable: a callback that stops failing has to
     * produce a normal cavity, i.e. the failure latch is per-build state */
    ctx.fail_after = INT_MAX;
    atomic_store(&ctx.calls, 0);
    moist_update_cavity(error, cav, mol);
    if (moist_check_error(error)) {
        printf("  FAIL: the cavity stayed broken after the callback recovered\n");
        show_error(error);
        result = 1;
    } else {
        int ngrid = 0, nsph = 0;
        moist_get_cavity_sizes(error, cav, &ngrid, &nsph);
        if (moist_check_error(error) || ngrid <= 0 || nsph != natoms) {
            printf("  FAIL: unexpected cavity sizes after recovery\n");
            result = 1;
        }
    }

    /* Upper bound on the aborted build: stopping at evaluation 51 must cost a
     * small fraction of the full traversal (full_calls, from the recovery
     * build above), with only the calls already in flight on other threads
     * added on top */
    const int full_calls = atomic_load(&ctx.calls);
    printf("  full build for comparison: %d invocations\n", full_calls);
    if (result == 0 && total_calls >= full_calls / 2) {
        printf("  FAIL: the aborted build kept evaluating (%d of %d invocations)\n",
               total_calls, full_calls);
        result = 1;
    }

    if (result == 0) {
        printf("  Callback isodensity failure reported cleanly!\n");
    }

cleanup:
    moist_delete_cavity(&cav);
    moist_delete_structure(&mol);
    moist_delete_error(&error);
    return result;
}

/* Deprecated moist_update_drop_cavity rebuilds the cavity from scratch to
 * honour a new Lebedev order. Must not silently reset the cavity's configured
 * numerical settings to their compiled defaults */
int test_update_drop_cavity_keeps_params(void)
{
    printf("Start test: update_drop_cavity preserves configured parameters\n");

    int result = 1;

    moist_error error = moist_new_error();
    moist_structure mol = NULL;
    moist_cavity cav = NULL;

    const int natoms = H2O_NATOMS;

    int nleb = 26;
    const double tolerance_in = 1.0e-8;   /* not the default 1e-10 */

    mol = make_h2o(error);
    if (moist_check_error(error)) {
        show_error(error);
        goto cleanup;
    }

    cav = fixture_drop(error, &nleb, NULL, NULL, NULL, NULL,
                                NULL, NULL, NULL, &tolerance_in,
                                NULL, NULL, NULL, NULL, NULL);
    if (moist_check_error(error)) {
        show_error(error);
        goto cleanup;
    }

    moist_update_cavity(error, cav, mol);
    if (moist_check_error(error)) {
        show_error(error);
        goto cleanup;
    }

    /* Rebuild without changing the configured options */
    moist_update_cavity(error, cav, mol);
    if (moist_check_error(error)) {
        show_error(error);
        goto cleanup;
    }

    double tolerance_out = 0.0;
    moist_get_drop_cavity_tolerance(error, cav, &tolerance_out);
    if (moist_check_error(error)) {
        show_error(error);
        goto cleanup;
    }
    printf("  tolerance before/after rebuild: %.3e / %.3e\n", tolerance_in, tolerance_out);

    result = 0;
    if (tolerance_out != tolerance_in) {
        printf("  FAIL: rebuild reset the master tolerance to its default\n");
        result = 1;
    }

    /* The rebuilt cavity must still be usable; models are re-sourced from
     * detached copies, so a use-after-free would surface as garbage sizes or
     * a crash */
    int ngrid = 0, nsph = 0;
    moist_get_cavity_sizes(error, cav, &ngrid, &nsph);
    if (moist_check_error(error) || ngrid <= 0 || nsph != natoms) {
        show_error(error);
        printf("  FAIL: cavity unusable after rebuild (ngrid = %d, nsph = %d)\n", ngrid, nsph);
        result = 1;
    }

    if (result == 0) {
        printf("  update_drop_cavity preserved the configured parameters!\n");
    }

cleanup:
    moist_delete_cavity(&cav);
    moist_delete_structure(&mol);
    moist_delete_error(&error);
    return result;
}

/* Array-writing getters must refuse a capacity smaller than the cavity needs,
 * and accept one that is larger. Without the check, an undersized call writes
 * the cavity's full ngrid elements into a buffer the caller declared smaller.
 *
 * Buffers here are allocated at full size and handed a smaller logical
 * capacity: stamp every byte with a sentinel beforehand, then verify the
 * whole region afterward. A rejected call must write nothing at all,
 * anywhere */

/* Not a value any getter would plausibly write: neither a small integer, nor a
 * zero double, nor a valid bool */
#define CAPACITY_SENTINEL 0xA5

/* Named result field API: enumerate what the cavity holds, then pull fields
 * by name. Replaces the per-field entry points; anchor_id, branch and
 * branch_count come straight from the cavity */
int test_cavity_fields(void)
{
    printf("\n=== Test: named cavity result fields ===\n");

    moist_error error = moist_new_error();
    moist_structure mol = NULL;
    moist_cavity cav = NULL;
    int result = 1;
    char *about = NULL;
    int* numbering = NULL;
    int* anchor_id = NULL;
    int* branch = NULL;
    int* branch_count = NULL;
    bool* converged = NULL;
    bool* converged2 = NULL;

    mol = make_h2o(error);
    if (moist_check_error(error)) { show_error(error); goto cleanup; }

    cav = fixture_drop(error, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
                                NULL, NULL, NULL, NULL, NULL, NULL);
    if (moist_check_error(error)) { show_error(error); goto cleanup; }

    moist_update_cavity(error, cav, mol);
    if (moist_check_error(error)) { show_error(error); goto cleanup; }

    int ngrid = 0, nsph = 0;
    moist_get_cavity_sizes(error, cav, &ngrid, &nsph);
    if (moist_check_error(error)) { show_error(error); goto cleanup; }

    int nfield = 0;
    moist_get_cavity_field_count(error, cav, &nfield);
    if (moist_check_error(error)) { show_error(error); goto cleanup; }
    printf("  Cavity declares %d fields\n", nfield);
    if (nfield <= 0) {
        printf("  FAIL: a built cavity declares no fields\n");
        goto cleanup;
    }

    /* Every declared field must describe itself consistently: the reported
     * count is what an accessor writes, and must equal the product of dims */
    int saw_numbering = 0, saw_xyz = 0;
    for (int i = 0; i < nfield; i++) {
        char name[MOIST_FIELD_NAME_MAX + 1] = {0};
        int dtype = 0, rank = 0, count = 0;
        int dims[MOIST_FIELD_MAX_RANK] = {0};

        moist_get_cavity_field_info(error, cav, i, name,
                                    &dtype, &rank, dims, &count);
        if (moist_check_error(error)) { show_error(error); goto cleanup; }

        int expect = 1;
        for (int d = 0; d < rank; d++) expect *= dims[d];
        if (count != expect) {
            printf("  FAIL: field '%s' reports count %d for shape product %d\n",
                   name, count, expect);
            goto cleanup;
        }
        if (dtype != MOIST_FIELD_REAL && dtype != MOIST_FIELD_INT &&
            dtype != MOIST_FIELD_BOOL) {
            printf("  FAIL: field '%s' reports unknown type tag %d\n", name, dtype);
            goto cleanup;
        }
        if (strcmp(name, "numbering") == 0) saw_numbering = 1;
        if (strcmp(name, "xyz") == 0) {
            saw_xyz = 1;
            if (rank != 2 || dims[0] != ngrid || dims[1] != 3) {
                printf("  FAIL: xyz declared as rank %d (%d, %d)\n",
                       rank, dims[0], dims[1]);
                goto cleanup;
            }
        }
    }
    if (!saw_xyz || !saw_numbering) {
        printf("  FAIL: cavity does not declare xyz and numbering\n");
        goto cleanup;
    }

    size_t about_length = 0;
    moist_get_cavity_field_about(error, cav, "branch_count", NULL, 0, &about_length);
    if (moist_check_error(error) || about_length == 0) goto cleanup;
    size_t about_capacity = about_length + 1;
    about = malloc(about_capacity);
    if (!about) goto cleanup;
    size_t copied_length = 0;
    moist_get_cavity_field_about(error, cav, "branch_count", about, about_capacity, &copied_length);
    if (copied_length != about_length || strlen(about) != about_length) goto cleanup;
    if (moist_check_error(error)) { show_error(error); goto cleanup; }
    printf("  branch_count: %s\n", about);
    about[0] = 'X';
    copied_length = 123;
    moist_get_cavity_field_about(error, cav, "branch_count", about, SIZE_MAX, &copied_length);
    if (!moist_check_error(error) || about[0] != '\0' || copied_length != 123) goto cleanup;
    moist_get_cavity_field_about(error, cav, "branch_count", about, 1, &copied_length);
    if (moist_check_error(error) || about[0] != '\0' || copied_length != about_length) goto cleanup;
    copied_length = 123;
    about[0] = 'X';
    moist_get_cavity_field_about(error, cav, "absent", about, about_capacity, &copied_length);
    if (!moist_check_error(error) || about[0] != '\0' || copied_length != 123) goto cleanup;

    /* The branching data the cavity maintains, read directly */
    numbering = (int*)malloc((size_t)ngrid * sizeof(int));
    anchor_id = (int*)malloc((size_t)ngrid * sizeof(int));
    branch = (int*)malloc((size_t)ngrid * sizeof(int));
    branch_count = (int*)malloc((size_t)ngrid * sizeof(int));
    if (!numbering || !anchor_id || !branch || !branch_count) {
        printf("  Error: memory allocation failed\n");
        goto cleanup;
    }

    int num_leb = 0;
    if (read_field_int(error, cav, "num_leb", &num_leb) ||
        read_field_int(error, cav, "numbering", numbering) ||
        read_field_int(error, cav, "anchor_id", anchor_id) ||
        read_field_int(error, cav, "branch", branch) ||
        read_field_int(error, cav, "branch_count", branch_count)) {
        goto cleanup;
    }

    /* numbering is the packing of the other two: the cavity's own arrays and
     * the packed id must agree point for point */
    int nbranched = 0;
    for (int i = 0; i < ngrid; i++) {
        const int base = nsph * num_leb;
        if (numbering[i] != anchor_id[i] + base * (branch[i] - 1)) {
            printf("  FAIL: numbering[%d] = %d disagrees with anchor %d branch %d\n",
                   i, numbering[i], anchor_id[i], branch[i]);
            goto cleanup;
        }
        if (branch[i] < 1 || branch_count[i] < 1) {
            printf("  FAIL: point %d has branch %d of %d\n",
                   i, branch[i], branch_count[i]);
            goto cleanup;
        }
        if (branch_count[i] > 1) nbranched++;
    }
    printf("  %d of %d points are branched\n", nbranched, ngrid);

    /* Logical accessor, the third of the three: every value it can return is
     * also a legal seed. Read into two oppositely seeded buffers to confirm
     * every entry was written */
    converged = (bool*)malloc((size_t)ngrid * sizeof(bool));
    converged2 = (bool*)malloc((size_t)ngrid * sizeof(bool));
    if (!converged || !converged2) {
        printf("  Error: memory allocation failed\n");
        goto cleanup;
    }
    memset(converged, 0, (size_t)ngrid * sizeof(bool));
    memset(converged2, 1, (size_t)ngrid * sizeof(bool));

    moist_get_cavity_field_bool(error, cav, "converged", converged);
    if (moist_check_error(error)) { show_error(error); goto cleanup; }
    moist_get_cavity_field_bool(error, cav, "converged", converged2);
    if (moist_check_error(error)) { show_error(error); goto cleanup; }

    int nconverged = 0;
    for (int i = 0; i < ngrid; i++) {
        if (converged[i] != converged2[i]) {
            printf("  FAIL: converged[%d] depends on the seed of the buffer\n", i);
            goto cleanup;
        }
        if (converged[i]) nconverged++;
    }
    printf("  %d of %d points converged\n", nconverged, ngrid);

    /* A name the cavity does not hold is an error, never a buffer of zeros */
    double probe = 12345.0;
    moist_get_cavity_field_real(error, cav, "not_a_field", &probe);
    if (!moist_check_error(error)) {
        printf("  FAIL: an unknown field name was accepted\n");
        goto cleanup;
    }
    if (probe != 12345.0) {
        printf("  FAIL: a rejected field read wrote into the buffer\n");
        goto cleanup;
    }
    moist_delete_error(&error);
    error = moist_new_error();

    /* Same for a field whose optional property was never requested: the
     * default cavity computes no curvature, k1 is absent rather than zero */
    moist_get_cavity_field_real(error, cav, "k1", &probe);
    if (!moist_check_error(error)) {
        printf("  FAIL: a field that was never computed was handed out\n");
        goto cleanup;
    }
    if (probe != 12345.0) {
        printf("  FAIL: an uncomputed field read wrote into the buffer\n");
        goto cleanup;
    }
    moist_delete_error(&error);
    error = moist_new_error();

    result = 0;
    printf("  Named field API behaved as declared\n");

cleanup:
    free(about);
    free(numbering); free(anchor_id); free(branch); free(branch_count);
    free(converged); free(converged2);
    moist_delete_cavity(&cav);
    moist_delete_structure(&mol);
    moist_delete_error(&error);
    return result;
}

static void fill_sentinel(void* p, size_t nbytes)
{
    memset(p, CAPACITY_SENTINEL, nbytes);
}

/* Non-zero if every byte of the region still carries the sentinel */
static int sentinel_intact(const void* p, size_t nbytes)
{
    const unsigned char* b = (const unsigned char*)p;
    for (size_t i = 0; i < nbytes; i++) {
        if (b[i] != (unsigned char)CAPACITY_SENTINEL) return 0;
    }
    return 1;
}

static int expect_failure(moist_error* error, const char* fragment, const char* label);

int test_capacity_validation(void)
{
    printf("Start test: array capacity validation\n");

    int result = 1;

    moist_error error = moist_new_error();
    moist_structure mol = NULL;
    moist_cavity cav = NULL;

    double* xyz = NULL;
    double* a = NULL;
    int* owner = NULL;
    bool* converged = NULL;
    double* vradii = NULL;
    double* asph = NULL;
    double* wleb = NULL;
    double* r_iI0 = NULL;
    double* f = NULL;
    double* rho = NULL;
    int* numbering = NULL;
    double* amat0 = NULL;
    double* xi = NULL;
    double* gaussian_f = NULL;
    double* big_xyz = NULL;
    double* big_a = NULL;
    int* big_owner = NULL;
    bool* big_conv = NULL;
    double* big_radii = NULL;
    double* big_asph = NULL;

    const int natoms = H2O_NATOMS;

    mol = make_h2o(error);
    cav = fixture_drop(error, NULL, NULL, NULL, NULL, NULL,
                                NULL, NULL, NULL, NULL,
                                NULL, NULL, NULL, NULL, NULL);
    moist_update_cavity(error, cav, mol);

    int ngrid = 0, nsph = 0;
    moist_get_cavity_sizes(error, cav, &ngrid, &nsph);
    if (moist_check_error(error) || ngrid <= 1 || nsph != natoms) {
        show_error(error);
        goto cleanup;
    }

    /* Physical extents: what an unchecked getter would write */
    const size_t n_xyz = (size_t)3 * (size_t)ngrid * sizeof(double);
    const size_t n_a = (size_t)ngrid * sizeof(double);
    const size_t n_owner = (size_t)ngrid * sizeof(int);
    const size_t n_conv = (size_t)ngrid * sizeof(bool);
    const size_t n_sph = (size_t)nsph * sizeof(double);
    const size_t n_amat = (size_t)ngrid * (size_t)ngrid * sizeof(double);

    xyz = (double*)malloc(n_xyz);
    a = (double*)malloc(n_a);
    owner = (int*)malloc(n_owner);
    converged = (bool*)malloc(n_conv);
    vradii = (double*)malloc(n_sph);
    asph = (double*)malloc(n_sph);
    wleb = (double*)malloc(n_a);
    r_iI0 = (double*)malloc(n_a);
    f = (double*)malloc(n_a);
    rho = (double*)malloc(n_a);
    numbering = (int*)malloc(n_owner);
    amat0 = (double*)malloc(n_amat);
    xi = (double*)malloc(n_a);
    gaussian_f = (double*)malloc(n_a);

    if (!xyz || !a || !owner || !converged || !vradii || !asph || !wleb ||
        !r_iI0 || !f || !rho || !numbering || !amat0 || !xi ||
        !gaussian_f) {
        printf("  Error: memory allocation failed\n");
        goto cleanup;
    }

    /* Logical capacity offered to the library: one short of the truth */
    const int short_ngrid = ngrid - 1;

    double area = 0.0, volume = 0.0;
    int out_ngrid = 0, out_nsph = 0;

    result = 0;

    if (result == 0) {
        fill_sentinel(xyz, n_xyz); fill_sentinel(a, n_a);
        fill_sentinel(owner, n_owner); fill_sentinel(converged, n_conv);
        fill_sentinel(vradii, n_sph); fill_sentinel(asph, n_sph);

        moist_get_cavity_results(error, cav, short_ngrid, nsph,
                                 &area, &volume, &out_ngrid, &out_nsph,
                                 xyz, a, owner, converged, vradii, asph);
        if (!moist_check_error(error)) {
            printf("  FAIL: get_cavity_results accepted ngrid_cap = %d < %d\n",
                   short_ngrid, ngrid);
            result = 1;
        } else {
            show_error(error);
            if (!sentinel_intact(xyz, n_xyz) || !sentinel_intact(a, n_a) ||
                !sentinel_intact(owner, n_owner) ||
                !sentinel_intact(converged, n_conv) ||
                !sentinel_intact(vradii, n_sph) || !sentinel_intact(asph, n_sph)) {
                printf("  FAIL: get_cavity_results wrote into a rejected buffer\n");
                result = 1;
            }
            moist_delete_error(&error);
            error = moist_new_error();
        }
    }

    if (result == 0) {
        /* A short sphere capacity has to be caught just the same; the grid
         * capacity is honest here, isolating the sphere check */
        fill_sentinel(xyz, n_xyz); fill_sentinel(a, n_a);
        fill_sentinel(owner, n_owner); fill_sentinel(converged, n_conv);
        fill_sentinel(vradii, n_sph); fill_sentinel(asph, n_sph);

        moist_get_cavity_results(error, cav, ngrid, nsph - 1,
                                 &area, &volume, &out_ngrid, &out_nsph,
                                 xyz, a, owner, converged, vradii, asph);
        if (!moist_check_error(error)) {
            printf("  FAIL: get_cavity_results accepted nsph_cap = %d < %d\n",
                   nsph - 1, nsph);
            result = 1;
        } else {
            if (!sentinel_intact(xyz, n_xyz) || !sentinel_intact(a, n_a) ||
                !sentinel_intact(owner, n_owner) ||
                !sentinel_intact(converged, n_conv) ||
                !sentinel_intact(vradii, n_sph) || !sentinel_intact(asph, n_sph)) {
                printf("  FAIL: get_cavity_results wrote into a rejected buffer\n");
                result = 1;
            }
            moist_delete_error(&error);
            error = moist_new_error();
        }
    }

    /* Reading a real field through the logical accessor must be refused too:
     * the payload the accessor would copy out is the unallocated one, so a
     * mistyped read must not become a memory error */
    if (result == 0) {
        fill_sentinel(converged, n_conv);

        moist_get_cavity_field_bool(error, cav, "wleb", converged);
        if (!moist_check_error(error)) {
            printf("  FAIL: get_cavity_field_bool accepted a real-valued field\n");
            result = 1;
        } else {
            if (!sentinel_intact(converged, n_conv)) {
                printf("  FAIL: get_cavity_field_bool wrote into a rejected buffer\n");
                result = 1;
            }
            moist_delete_error(&error);
            error = moist_new_error();
        }
    }

    /* Reading a real field through the integer accessor must be refused too:
     * the type tag is part of the contract, not a hint */
    if (result == 0) {
        fill_sentinel(numbering, n_owner);

        moist_get_cavity_field_int(error, cav, "wleb", numbering);
        if (!moist_check_error(error)) {
            printf("  FAIL: get_cavity_field_int accepted a real-valued field\n");
            result = 1;
        } else {
            if (!sentinel_intact(numbering, n_owner)) {
                printf("  FAIL: get_cavity_field_int wrote into a rejected buffer\n");
                result = 1;
            }
            moist_delete_error(&error);
            error = moist_new_error();
        }
    }

    if (result == 0) {
        fill_sentinel(amat0, n_amat); fill_sentinel(xi, n_a);

        moist_assemble_amat(error, cav, short_ngrid, amat0, xi);
        if (!moist_check_error(error)) {
            printf("  FAIL: assemble_amat accepted a short capacity\n");
            result = 1;
        } else {
            if (!sentinel_intact(amat0, n_amat) || !sentinel_intact(xi, n_a)) {
                printf("  FAIL: assemble_amat wrote into a rejected buffer\n");
                result = 1;
            }
            moist_delete_error(&error);
            error = moist_new_error();
        }
    }

    if (result == 0) {
        fill_sentinel(xi, n_a); fill_sentinel(gaussian_f, n_a);

        moist_get_cavity_gaussian(error, cav, short_ngrid, xi, gaussian_f);
        if (!moist_check_error(error)) {
            printf("  FAIL: get_cavity_gaussian accepted a short capacity\n");
            result = 1;
        } else {
            if (!sentinel_intact(xi, n_a) || !sentinel_intact(gaussian_f, n_a)) {
                printf("  FAIL: get_cavity_gaussian wrote into a rejected buffer\n");
                result = 1;
            }
            moist_delete_error(&error);
            error = moist_new_error();
        }
    }

    /* Oversized buffers are accepted; only the returned logical extent is valid */
    if (result == 0) {
        const int big_ngrid = ngrid + 17;
        const int big_nsph = nsph + 5;
        big_xyz = (double*)calloc((size_t)3 * big_ngrid, sizeof(double));
        big_a = (double*)calloc((size_t)big_ngrid, sizeof(double));
        big_owner = (int*)calloc((size_t)big_ngrid, sizeof(int));
        big_conv = (bool*)calloc((size_t)big_ngrid, sizeof(bool));
        big_radii = (double*)calloc((size_t)big_nsph, sizeof(double));
        big_asph = (double*)calloc((size_t)big_nsph, sizeof(double));

        if (!big_xyz || !big_a || !big_owner || !big_conv || !big_radii || !big_asph) {
            printf("  Error: memory allocation failed\n");
            result = 1;
        } else {
            moist_get_cavity_results(error, cav, big_ngrid, big_nsph,
                                     &area, &volume, &out_ngrid, &out_nsph,
                                     big_xyz, big_a, big_owner, big_conv,
                                     big_radii, big_asph);
            if (moist_check_error(error) || out_ngrid != ngrid || out_nsph != nsph)
                result = 1;
            for (int i = 3 * ngrid; i < 3 * big_ngrid; ++i)
                if (big_xyz[i] != 0.0) result = 1;
            for (int i = nsph; i < big_nsph; ++i)
                if (big_radii[i] != 0.0) result = 1;
        }
    }

    if (result == 0) {
        printf("  Array capacity validation behaved as specified!\n");
    }

cleanup:
    free(xyz); free(a); free(owner); free(converged); free(vradii); free(asph);
    free(wleb); free(r_iI0); free(f); free(rho); free(numbering);
    free(amat0); free(xi); free(gaussian_f);
    free(big_xyz); free(big_a); free(big_owner); free(big_conv);
    free(big_radii); free(big_asph);

    moist_delete_cavity(&cav);
    moist_delete_structure(&mol);
    moist_delete_error(&error);
    return result;
}

/* test_registry (below) records each test by name; "Start test:" lines
 * interleave with everything else the tests print */
/* ------------------------------------------------------------------------- */
/* Host coupling protocol                                                    */
/* ------------------------------------------------------------------------- */

/* Raw Gaussian-potential primitives from the water nuclei. A QM host adds
 * the electronic integral contribution to each of these quantities */
static void gaussian_nuclear_primitives(int ngrid, const double* xyz, const double* xi,
                                        double* phi, double* dphi, double* dxi)
{
    const double two_sqrt_pi = 2.0 / sqrt(acos(-1.0));
    for (int i = 0; i < ngrid; ++i) {
        phi[i] = 0.0;
        dxi[i] = 0.0;
        for (int k = 0; k < 3; ++k) dphi[3*i+k] = 0.0;
        for (int a = 0; a < H2O_NATOMS; ++a) {
            double d[3];
            for (int k = 0; k < 3; ++k) d[k] = xyz[3*i+k] - h2o_positions[3*a+k];
            double r = sqrt(d[0]*d[0] + d[1]*d[1] + d[2]*d[2]);
            double x = xi[i]*r;
            double z = (double)h2o_numbers[a];
            double radial = fabs(x) < 1e-3
                ? (2.0/3.0)*two_sqrt_pi*xi[i]*xi[i]*xi[i]*(1.0-0.6*x*x+3.0*x*x*x*x/14.0)
                : (erf(x)-two_sqrt_pi*x*exp(-x*x))/(r*r*r);
            phi[i] += z*(r == 0.0 ? two_sqrt_pi*xi[i] : erf(x)/r);
            dxi[i] += z*two_sqrt_pi*exp(-x*x);
            for (int k = 0; k < 3; ++k) dphi[3*i+k] -= z*radial*d[k];
        }
    }
}

/* Consume an EXPECTED failure: the handle must carry an error whose message
 * names `fragment`; the handle is then replaced by a clean one so the test can
 * go on. Returns 1 when the failure did not happen or has the wrong wording */
static int expect_failure(moist_error* error, const char* fragment, const char* what)
{
    if (!moist_check_error(*error)) {
        printf("  FAIL: %s succeeded, an error was expected\n", what);
        return 1;
    }
    char message[512] = {0};
    const int message_len = (int)sizeof(message);
    moist_get_error(*error, message, &message_len);
    int ok = strstr(message, fragment) != NULL;
    printf("  %s: %s\n    -> %s\n", ok ? "expected failure" : "FAIL: wrong message for",
           what, message);
    moist_delete_error(error);
    *error = moist_new_error();
    return ok ? 0 : 1;
}

/* Every diagnostic identifies the function the C caller actually invoked */
static int error_has_origin(moist_error error, const char *function)
{
    char message[512] = {0}, prefix[128];
    const int capacity = sizeof message;
    if (!moist_check_error(error)) return 0;
    moist_get_error(error, message, &capacity);
    snprintf(prefix, sizeof prefix, "[%s] ", function);
    return strncmp(message, prefix, strlen(prefix)) == 0
        && strstr(message + strlen(prefix), "[moist_") == NULL;
}

int test_error_origins(void)
{
    moist_error error = moist_new_error();
    moist_structure mol = NULL;
    moist_component component = NULL;
    moist_lsf lsf = NULL;
    moist_cavity cavity = NULL;
    const int numbers[2] = {1, 1};
    const double coincident[6] = {0}, positions[6] = {0, 0, 0, 2, 0, 0};
    moist_pcm_options pcm = {0};
    moist_drop_options drop = {0};
    moist_iswig_options iswig = {0};
    moist_svdw_options svdw = {0};
    moist_cfc_options cfc = {0};
    moist_isodensity_options iso = {0};
    moist_model_options model_options = {0};
    moist_model model = NULL;
    moist_coupling coupling = NULL;
    char request_name[MOIST_NAME_MAX + 1];
    bool flag = false;
    double value = 0.0;
    int failed = 1;
#define REQUIRE_ORIGIN(name) do { if (!error_has_origin(error, name)) goto cleanup; } while (0)
#define CHECK_INIT(name) do { moist_init_##name##_options(error, NULL, 0); \
    REQUIRE_ORIGIN("moist_init_" #name "_options"); } while (0)
    CHECK_INIT(drop);
    CHECK_INIT(iswig);
    CHECK_INIT(svdw);
    CHECK_INIT(cfc);
    CHECK_INIT(isodensity);
    CHECK_INIT(model);
    CHECK_INIT(pcm);
    component = moist_new_cpcm_component(error, 32, &pcm);
    REQUIRE_ORIGIN("moist_new_cpcm_component");
    component = moist_new_cosmo_component(error, 32, &pcm);
    REQUIRE_ORIGIN("moist_new_cosmo_component");
    lsf = moist_new_svdw_lsf(error, &svdw);
    REQUIRE_ORIGIN("moist_new_svdw_lsf");
    lsf = moist_new_cfc_lsf(error, &cfc);
    REQUIRE_ORIGIN("moist_new_cfc_lsf");
    lsf = moist_new_isodensity_lsf(error, 0, NULL, NULL, NULL, NULL, NULL, &iso);
    REQUIRE_ORIGIN("moist_new_isodensity_lsf");
    lsf = moist_new_isodensity_callback_lsf(error, NULL, NULL, &iso);
    REQUIRE_ORIGIN("moist_new_isodensity_callback_lsf");
    cavity = moist_new_drop_cavity(error, NULL, NULL, &drop);
    REQUIRE_ORIGIN("moist_new_drop_cavity");
    cavity = moist_new_iswig_cavity(error, NULL, &iswig);
    REQUIRE_ORIGIN("moist_new_iswig_cavity");
    model = moist_new_model(error, NULL, &model_options);
    REQUIRE_ORIGIN("moist_new_model");
    moist_init_isodensity_options(error, &iso, sizeof iso);
    iso.rho_iso = -1;
    lsf = moist_new_isodensity_lsf(error, 0, NULL, NULL, NULL, NULL, NULL, &iso);
    REQUIRE_ORIGIN("moist_new_isodensity_lsf");
    mol = moist_new_structure(error, 2, numbers, coincident, NULL, NULL);
    REQUIRE_ORIGIN("moist_new_structure");
    moist_delete_structure(&mol);
    mol = moist_new_structure(error, 2, numbers, positions, NULL, NULL);
    if (!mol || moist_check_error(error)) goto cleanup;
    moist_update_structure(error, mol, coincident, NULL);
    REQUIRE_ORIGIN("moist_update_structure");
    moist_update_model(error, NULL, mol);
    REQUIRE_ORIGIN("moist_update_model");
    moist_get_model_energy(error, NULL, NULL, NULL);
    REQUIRE_ORIGIN("moist_get_model_energy");
    moist_prepare_model_energy(error, NULL, NULL);
    REQUIRE_ORIGIN("moist_prepare_model_energy");
    moist_prepare_model_response(error, NULL, NULL);
    REQUIRE_ORIGIN("moist_prepare_model_response");
    moist_prepare_model_gradient(error, NULL, NULL);
    REQUIRE_ORIGIN("moist_prepare_model_gradient");
    moist_get_model_response(error, NULL, NULL, NULL);
    REQUIRE_ORIGIN("moist_get_model_response");
    moist_get_model_gradient(error, NULL, NULL, NULL, 0, NULL);
    REQUIRE_ORIGIN("moist_get_model_gradient");
    coupling = moist_new_coupling(error, NULL);
    REQUIRE_ORIGIN("moist_new_coupling");
    /* A failing walk reports and still ends the host loop */
    if (moist_next_coupling_request(error, NULL)) goto cleanup;
    REQUIRE_ORIGIN("moist_next_coupling_request");
    moist_get_coupling_request_name(error, NULL, request_name);
    REQUIRE_ORIGIN("moist_get_coupling_request_name");
    moist_get_coupling_request_missing(error, NULL, "phi", &flag);
    REQUIRE_ORIGIN("moist_get_coupling_request_missing");
    moist_answer_coupling_request(error, NULL, "phi", &value);
    REQUIRE_ORIGIN("moist_answer_coupling_request");
    moist_get_coupling_request_width(error, NULL, &value);
    REQUIRE_ORIGIN("moist_get_coupling_request_width");
    /* So does a failing response walk */
    if (moist_next_response_item(error, NULL)) goto cleanup;
    REQUIRE_ORIGIN("moist_next_response_item");
    moist_get_response_item_name(error, NULL, request_name);
    REQUIRE_ORIGIN("moist_get_response_item_name");
    moist_get_response_array(error, NULL, "w_phi", &value);
    REQUIRE_ORIGIN("moist_get_response_array");
    moist_contract_surface_lsf_weights(error, NULL, NULL, NULL, NULL, NULL, NULL, NULL);
    REQUIRE_ORIGIN("moist_contract_surface_lsf_weights");
    moist_contract_surface_lsf_weights_extended(error, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL);
    REQUIRE_ORIGIN("moist_contract_surface_lsf_weights_extended");
    failed = 0;
cleanup:
    if (failed) show_error(error);
    moist_delete_structure(&mol);
    moist_delete_component(&component);
    moist_delete_cavity(&cavity);
    moist_delete_lsf(&lsf);
    moist_delete_coupling(&coupling);
    moist_delete_model(&model);
    moist_delete_error(&error);
    return failed;
#undef CHECK_INIT
#undef REQUIRE_ORIGIN
}

/* The refusal of every request query outside a next() window */
static const char no_current_request[] = "No current coupling request - call next() first";

/* Native rejection invalidates model state even when Python would reject first */
int test_model_density_rejection(void)
{
    moist_error error = moist_new_error();
    moist_cavity cavity = NULL;
    moist_structure mol = NULL;
    moist_model model = NULL;
    moist_component pcm = NULL;
    moist_coupling coupling = NULL;
    const int atom = 0, angular = 0, primitives = 1, nleb = 26, hydrogen = 1;
    const double exponent = .5, coefficient = 1, density = 1, position[3] = {0};
    const double bad_density[4] = {1, 0, 0, 1};
    double energy = 2;
    bool missing = true;
    int failed = 1;
    cavity = fixture_isodensity_internal(error, 1, &atom, &angular, &primitives,
        &exponent, &coefficient, 1e-3, NULL, &nleb, NULL, NULL, NULL, NULL, NULL);
    if (!cavity || moist_check_error(error)) goto cleanup;
    moist_set_isodensity_density(error, cavity, 1, &density);
    if (moist_check_error(error)) goto cleanup;
    mol = moist_new_structure(error, 1, &hydrogen, position, NULL, NULL);
    if (!mol || moist_check_error(error)) goto cleanup;
    model = moist_new_model(error, cavity, NULL);
    pcm = moist_new_cpcm_component(error, 32, NULL);
    if (!model || !pcm || moist_check_error(error)) goto cleanup;
    moist_add_model_component(error, model, pcm);
    if (moist_check_error(error)) goto cleanup;
    moist_update_model(error, model, mol);
    if (moist_check_error(error)) goto cleanup;
    coupling = moist_new_coupling(error, model);
    moist_prepare_model_energy(error, model, coupling);
    if (moist_check_error(error)) goto cleanup;
    /* Walk to the request; the rejection needs a current one to end */
    if (!moist_next_coupling_request(error, coupling) || moist_check_error(error)) goto cleanup;
    moist_set_model_isodensity_density(error, model, 2, bad_density);
    if (!error_has_origin(error, "moist_set_model_isodensity_density")) goto cleanup;
    if (expect_failure(&error, "shape", "invalid model density")) goto cleanup;
    moist_get_coupling_request_missing(error, coupling, "phi", &missing);
    if (!error_has_origin(error, "moist_get_coupling_request_missing")) goto cleanup;
    if (expect_failure(&error, no_current_request,
                       "density rejection ends the walk")) goto cleanup;
    if (!missing) goto cleanup;
    /* The staging went with it: nothing to walk until the next prepare */
    if (moist_next_coupling_request(error, coupling) || moist_check_error(error)) goto cleanup;
    moist_get_model_energy(error, model, coupling, &energy);
    if (!error_has_origin(error, "moist_get_model_energy") || energy != 2) goto cleanup;
    failed = 0;
cleanup:
    if (failed) show_error(error);
    moist_delete_coupling(&coupling);
    moist_delete_model(&model);
    moist_delete_component(&pcm);
    moist_delete_cavity(&cavity);
    moist_delete_structure(&mol);
    moist_delete_error(&error);
    return failed;
}

/* Leave a protocol test at its `cleanup` label when an expectation fails,
 * naming the condition. Expected failures go through expect_failure */
#define REQUIRE(cond) do { if (!(cond)) { \
    printf("  FAIL: %s (line %d)\n", #cond, __LINE__); goto cleanup; } } while (0)

/* The current request carries the name `expected` */
static int request_is(moist_error error, moist_coupling cpl, const char* expected)
{
    char name[MOIST_NAME_MAX + 1] = {0};
    moist_get_coupling_request_name(error, cpl, name);
    if (moist_check_error(error)) return 0;
    if (strcmp(name, expected) != 0) {
        printf("  FAIL: the current request is '%s', not '%s'\n", name, expected);
        return 0;
    }
    return 1;
}

/* The current response item carries the name `expected` */
static int item_is(moist_error error, moist_response response, const char* expected)
{
    char name[MOIST_NAME_MAX + 1] = {0};
    moist_get_response_item_name(error, response, name);
    if (moist_check_error(error)) return 0;
    if (strcmp(name, expected) != 0) {
        printf("  FAIL: the current response item is '%s', not '%s'\n", name, expected);
        return 0;
    }
    return 1;
}

/* Output `output` of the current request is missing exactly when `expected`.
 * The flag starts at the opposite value; an unwritten result then fails */
static int missing_is(moist_error error, moist_coupling cpl, const char* output, bool expected)
{
    bool missing = !expected;
    moist_get_coupling_request_missing(error, cpl, output, &missing);
    if (moist_check_error(error)) return 0;
    if (missing != expected) {
        printf("  FAIL: '%s' is %smissing\n", output, missing ? "" : "not ");
        return 0;
    }
    return 1;
}

/* Everything the host computes on the cavity grid, by output name: the nuclear
 * potential primitives and stand-ins for the Gaussian moments of its density
 * (NULL when the model asks for no moments) */
struct host_outputs {
    int ngrid;
    const double *phi, *dphi_dr, *dphi_dxi;
    const double *gt, *pt, *mt, *rt;
};

/* One pass of the host loop as moist.h prints it: walk the requests, dispatch
 * on the name and answer every missing output. `visited` receives the request
 * names in visiting order, comma separated. Returns the number of requests
 * visited, or -1 on a failure */
static int answer_pass(moist_error error, moist_coupling cpl, const struct host_outputs* host,
                       char* visited, size_t visited_cap)
{
    static const char* const potential[] = {"phi", "dphi_dr", "dphi_dxi"};
    static const char* const moments[] = {"gt", "pt", "mt", "rt"};
    const double* potential_values[] = {host->phi, host->dphi_dr, host->dphi_dxi};
    const double* moment_values[] = {host->gt, host->pt, host->mt, host->rt};
    int visits = 0;
    visited[0] = '\0';
    while (moist_next_coupling_request(error, cpl)) {
        char name[MOIST_NAME_MAX + 1];
        const char* const* outputs = potential;
        const double* const* values = potential_values;
        size_t noutputs = 3;
        moist_get_coupling_request_name(error, cpl, name);
        if (moist_check_error(error)) return -1;
        if (strcmp(name, "gaussian_moments") == 0) {
            outputs = moments;
            values = moment_values;
            noutputs = 4;
        } else if (strcmp(name, "gaussian_potential") != 0) {
            printf("  FAIL: unexpected request '%s'\n", name);
            return -1;
        }
        const size_t used = strlen(visited);
        snprintf(visited + used, visited_cap - used, "%s%s", visits++ ? "," : "", name);
        for (size_t k = 0; k < noutputs; ++k) {
            bool missing = false;
            moist_get_coupling_request_missing(error, cpl, outputs[k], &missing);
            if (moist_check_error(error)) return -1;
            if (!missing) continue;
            if (!values[k]) {
                printf("  FAIL: '%s' is missing, but the host has no value for it\n", outputs[k]);
                return -1;
            }
            moist_answer_coupling_request(error, cpl, outputs[k], values[k]);
            if (moist_check_error(error)) return -1;
        }
    }
    return moist_check_error(error) ? -1 : visits;
}

/* The refusal of every item query outside a next() window of the response */
static const char no_current_item[] = "No current response item - call next() first";

/* One pass of the response walk as moist.h prints it: `names` receives the
 * item names in visiting order, comma separated. A pass left open is finished
 * first; the listing always starts at the first item. Returns the number of
 * items visited, or -1 on a failure */
static int response_names(moist_error error, moist_response response, char* names,
                          size_t names_cap)
{
    int visits = 0;
    names[0] = '\0';
    while (moist_next_response_item(error, response)) {}
    if (moist_check_error(error)) return -1;
    while (moist_next_response_item(error, response)) {
        char name[MOIST_NAME_MAX + 1];
        moist_get_response_item_name(error, response, name);
        if (moist_check_error(error)) return -1;
        const size_t used = strlen(names);
        snprintf(names + used, names_cap - used, "%s%s", visits++ ? "," : "", name);
    }
    return moist_check_error(error) ? -1 : visits;
}

/* Walk the response to the item named `item`, as a host contracting it does.
 * A pass left open is finished first. Returns 1 with that item current, 0 when
 * the response has no such item (the pass then ended), -1 on a failure */
static int walk_to_item(moist_error error, moist_response response, const char* item)
{
    while (moist_next_response_item(error, response)) {}
    if (moist_check_error(error)) return -1;
    while (moist_next_response_item(error, response)) {
        char name[MOIST_NAME_MAX + 1];
        moist_get_response_item_name(error, response, name);
        if (moist_check_error(error)) return -1;
        if (strcmp(name, item) == 0) return 1;
    }
    return moist_check_error(error) ? -1 : 0;
}

/* Read one array of the current response item as a host should, checking the
 * size contract on the way: a NULL buffer is refused by name, and moist writes
 * exactly the documented ngrid * per_point values -- never past them into a
 * larger buffer. `per_point` counts the values per grid point (1, 3 or 9);
 * `values` receives ngrid * per_point of them */
static int read_response_array(moist_error* error, moist_response response, const char* array,
                               int ngrid, int per_point, double* values)
{
    const int padding = 5;
    const size_t logical = (size_t)ngrid * (size_t)per_point * sizeof(double);
    const size_t allocated = (size_t)(ngrid + padding) * (size_t)per_point * sizeof(double);
    char fragment[128], what[128];
    int failed = 1;
    double* buffer = (double*)malloc(allocated);
    if (!buffer) {
        printf("  Error: memory allocation failed\n");
        return 1;
    }

    fill_sentinel(buffer, allocated);
    moist_get_response_array(*error, response, array, NULL);
    snprintf(fragment, sizeof fragment, "Null array pointer provided for '%s'", array);
    snprintf(what, sizeof what, "NULL buffer for %s", array);
    if (expect_failure(error, fragment, what)) goto done;

    moist_get_response_array(*error, response, array, buffer);
    if (moist_check_error(*error)) goto done;
    if (!sentinel_intact((const unsigned char*)buffer + logical, allocated - logical)) {
        printf("  FAIL: reading %s wrote past the logical grid\n", array);
        goto done;
    }
    for (size_t k = 0; k < (size_t)ngrid * (size_t)per_point; ++k) {
        if (!isfinite(buffer[k])) {
            printf("  FAIL: %s holds a non-finite value\n", array);
            goto done;
        }
    }
    memcpy(values, buffer, logical);
    failed = 0;
done:
    free(buffer);
    return failed;
}

/* One fixed cavity, one CPCM component, one request: the host walk through the
 * three phases, with every refusal the protocol promises on the way */
static int run_coupling_protocol(const char* label, moist_cavity cav)
{
    printf("Start test: raw coupling protocol (%s)\n", label);
    int result = 1, ngrid = 0, nsph = 0, visits = 0;
    moist_error error = moist_new_error();
    moist_structure mol = make_h2o(error);
    moist_pcm_options pcm_options;
    moist_init_pcm_options(error, &pcm_options, sizeof pcm_options);
    pcm_options.solver = moist_pcm_solver_lu;
    moist_component pcm = moist_new_cpcm_component(error, 78.4, &pcm_options);
    moist_model model = moist_new_model(error, cav, NULL);
    moist_coupling cpl = NULL, other = NULL;
    moist_response response = NULL;
    moist_cavity borrowed = NULL;
    double *xyz = NULL, *phi = NULL, *dphi = NULL, *xi = NULL, *dxi = NULL, *w_phi = NULL;
    double energy = 0.0, gradient[3 * (H2O_NATOMS + 1)];
    char name[MOIST_NAME_MAX + 1], visited[128];
    bool more = false, missing = true;
    REQUIRE(mol && pcm && model && !moist_check_error(error));
    moist_add_model_component(error, model, pcm);
    REQUIRE(!moist_check_error(error));
    moist_update_model(error, model, mol);
    REQUIRE(!moist_check_error(error));
    cpl = moist_new_coupling(error, model);
    other = moist_new_coupling(error, model);
    response = moist_new_response(error);
    REQUIRE(cpl && other && response && !moist_check_error(error));
    /* A response nothing has filled has nothing to walk, which is no error */
    more = moist_next_response_item(error, response);
    REQUIRE(!moist_check_error(error) && !more);
    moist_get_model_energy(error, model, cpl, &energy);
    if (expect_failure(&error, "get_energy requires a coupling staged by prepare_energy; "
                       "this coupling is not staged", "evaluation before preparing")) goto cleanup;
    /* Grid inputs are read from the model's cavity */
    borrowed = moist_get_model_cavity(error, model);
    moist_get_cavity_sizes(error, borrowed, &ngrid, &nsph);
    REQUIRE(!moist_check_error(error) && ngrid > 1);
    xyz = malloc(3 * ngrid * sizeof(double));
    phi = malloc(ngrid * sizeof(double));
    dphi = malloc(3 * ngrid * sizeof(double));
    dxi = calloc(ngrid, sizeof(double));
    xi = malloc(ngrid * sizeof(double));
    w_phi = malloc(ngrid * sizeof(double));
    REQUIRE(xyz && phi && dphi && xi && dxi && w_phi);
    moist_get_cavity_field_real(error, borrowed, "xyz", xyz);
    REQUIRE(!moist_check_error(error));
    moist_get_cavity_field_real(error, borrowed, "xi0", xi);
    REQUIRE(!moist_check_error(error));
    moist_delete_cavity(&borrowed);
    gaussian_nuclear_primitives(ngrid, xyz, xi, phi, dphi, dxi);
    const struct host_outputs host = {ngrid, phi, dphi, dxi, NULL, NULL, NULL, NULL};

    /* Before a prepare there is nothing to walk. That is no error: the phase
     * getters report the staging by name */
    more = moist_next_coupling_request(error, cpl);
    REQUIRE(!moist_check_error(error) && !more);
    /* Without a current request nothing can be inspected or answered */
    moist_get_coupling_request_name(error, cpl, name);
    if (expect_failure(&error, no_current_request, "name without a walk")) goto cleanup;
    moist_get_coupling_request_missing(error, cpl, "phi", &missing);
    if (expect_failure(&error, no_current_request, "missing query without a walk")) goto cleanup;
    REQUIRE(missing);
    moist_answer_coupling_request(error, cpl, "phi", phi);
    if (expect_failure(&error, no_current_request, "answer without a walk")) goto cleanup;
    moist_get_coupling_request_width(error, cpl, w_phi);
    if (expect_failure(&error, no_current_request, "width without a walk")) goto cleanup;

    /* Every prepare restarts the walk, even in the middle of a pass */
    moist_prepare_model_energy(error, model, cpl);
    REQUIRE(!moist_check_error(error));
    more = moist_next_coupling_request(error, cpl);
    REQUIRE(!moist_check_error(error) && more);
    moist_prepare_model_energy(error, model, cpl);
    REQUIRE(!moist_check_error(error));
    moist_answer_coupling_request(error, cpl, "phi", phi);
    if (expect_failure(&error, no_current_request, "answer after a new prepare")) goto cleanup;

    /* The energy pass: one request, phi missing and nothing else */
    more = moist_next_coupling_request(error, cpl);
    REQUIRE(!moist_check_error(error) && more && request_is(error, cpl, "gaussian_potential"));
    REQUIRE(missing_is(error, cpl, "phi", true) && missing_is(error, cpl, "dphi_dr", false));
    /* A name the request does not declare is not missing; that is no error */
    REQUIRE(missing_is(error, cpl, "gt", false));
    moist_get_coupling_request_width(error, cpl, w_phi);
    if (expect_failure(&error, "gaussian_potential has no input 'width'",
                       "width of a potential request")) goto cleanup;
    /* Answers are checked by name and value; a refused answer drops the
     * previous one */
    moist_answer_coupling_request(error, cpl, "gt", phi);
    if (expect_failure(&error, "gaussian_potential has no output 'gt'",
                       "unknown output on answer")) goto cleanup;
    moist_answer_coupling_request(error, cpl, "phi", phi);
    REQUIRE(!moist_check_error(error));
    const double saved = phi[0];
    phi[0] = NAN;
    moist_answer_coupling_request(error, cpl, "phi", phi);
    phi[0] = saved;
    if (expect_failure(&error, "gaussian_potential: phi must be finite", "non-finite answer")) goto cleanup;
    REQUIRE(missing_is(error, cpl, "phi", true));
    moist_answer_coupling_request(error, cpl, "phi", NULL);
    if (expect_failure(&error, "Null array pointer provided for 'phi'", "NULL answer buffer")) goto cleanup;
    moist_answer_coupling_request(error, cpl, NULL, phi);
    if (expect_failure(&error, "Output name is missing", "NULL output name")) goto cleanup;
    moist_answer_coupling_request(error, cpl, "phi", phi);
    REQUIRE(!moist_check_error(error) && missing_is(error, cpl, "phi", false));
    /* One request, one visit: the pass ends, and the current request with it */
    more = moist_next_coupling_request(error, cpl);
    REQUIRE(!moist_check_error(error) && !more);
    moist_answer_coupling_request(error, cpl, "phi", phi);
    if (expect_failure(&error, no_current_request, "answer after the pass")) goto cleanup;
    /* A second loop starts a new pass, which has nothing left to retry */
    more = moist_next_coupling_request(error, cpl);
    REQUIRE(!moist_check_error(error) && !more);
    moist_get_model_energy(error, model, cpl, &energy);
    REQUIRE(!moist_check_error(error) && isfinite(energy));
    double contribution = energy;
    printf("  energy = %.12f\n", contribution);
    energy = 7.0;
    moist_get_model_energy(error, model, cpl, &energy);
    moist_get_model_energy(error, model, cpl, &energy);
    REQUIRE(!moist_check_error(error) && fabs(energy - (7.0 + 2.0 * contribution)) <= 1e-12);
    moist_get_model_response(error, model, cpl, response);
    if (expect_failure(&error, "get_response requires a coupling staged by prepare_response; "
                       "this coupling is staged for the energy phase",
                       "evaluation in the wrong phase")) goto cleanup;
    for (int k = 0; k < 3 * (H2O_NATOMS + 1); ++k) gradient[k] = -12345.0;
    moist_get_model_gradient(error, model, cpl, response, H2O_NATOMS + 1, gradient);
    if (expect_failure(&error, "get_gradient requires a coupling staged by prepare_gradient; "
                       "this coupling is staged for the energy phase",
                       "gradient in the wrong phase")) goto cleanup;
    for (int k = 0; k < 3 * (H2O_NATOMS + 1); ++k) REQUIRE(gradient[k] == -12345.0);

    /* The response phase reuses the energy answer: on a fixed cavity there is
     * nothing to walk, and the prepare ended the old walk */
    moist_prepare_model_response(error, model, cpl);
    REQUIRE(!moist_check_error(error));
    moist_get_coupling_request_name(error, cpl, name);
    if (expect_failure(&error, no_current_request, "name after a new prepare")) goto cleanup;
    more = moist_next_coupling_request(error, cpl);
    REQUIRE(!moist_check_error(error) && !more);
    moist_get_model_response(error, model, cpl, response);
    REQUIRE(!moist_check_error(error));
    /* The response is walked like the coupling. A fixed cavity with CPCM
     * hands back the potential adjoint alone */
    REQUIRE(response_names(error, response, visited, sizeof visited) == 1 &&
            strcmp(visited, "potential_adjoint") == 0);
    /* Outside a pass nothing is current, and nothing is written */
    fill_sentinel(w_phi, (size_t)ngrid * sizeof(double));
    moist_get_response_item_name(error, response, name);
    if (expect_failure(&error, no_current_item, "item name outside a pass")) goto cleanup;
    moist_get_response_array(error, response, "w_phi", w_phi);
    if (expect_failure(&error, no_current_item, "array outside a pass")) goto cleanup;
    /* Arrays are read by name from the current item: an array of another item
     * and an unknown one are refused, and nothing is written */
    REQUIRE(walk_to_item(error, response, "potential_adjoint") == 1);
    moist_get_response_array(error, response, "w_rho", w_phi);
    if (expect_failure(&error, "potential_adjoint has no array 'w_rho'",
                       "array of another item")) goto cleanup;
    moist_get_response_array(error, response, "w_bogus", w_phi);
    if (expect_failure(&error, "potential_adjoint has no array 'w_bogus'",
                       "unknown array")) goto cleanup;
    REQUIRE(sentinel_intact(w_phi, (size_t)ngrid * sizeof(double)));
    REQUIRE(!read_response_array(&error, response, "w_phi", ngrid, 1, w_phi));
    /* Past the only item the pass ends, and nothing is current any more */
    more = moist_next_response_item(error, response);
    REQUIRE(!moist_check_error(error) && !more);
    moist_get_response_array(error, response, "w_phi", w_phi);
    if (expect_failure(&error, no_current_item, "array after the pass")) goto cleanup;
    /* For CPCM E = q.phi / 2 with q = dE/dphi: the adjoint read back by name
     * has to reproduce the energy of the answer */
    double half_q_phi = 0.0;
    for (int i = 0; i < ngrid; ++i) half_q_phi += 0.5 * w_phi[i] * phi[i];
    printf("  q.phi/2 = %.12f\n", half_q_phi);
    REQUIRE(agrees_to(half_q_phi, contribution, 1e-10));

    /* The gradient phase keeps phi and adds the derivatives. A pass left early
     * is resumed by the next call; past the only request, that ends it */
    moist_prepare_model_gradient(error, model, cpl);
    REQUIRE(!moist_check_error(error));
    visits = 0;
    while (moist_next_coupling_request(error, cpl)) {
        visits++;
        REQUIRE(request_is(error, cpl, "gaussian_potential"));
        REQUIRE(missing_is(error, cpl, "phi", false) && missing_is(error, cpl, "dphi_dr", true) &&
                missing_is(error, cpl, "dphi_dxi", true));
        break;
    }
    REQUIRE(!moist_check_error(error) && visits == 1);
    more = moist_next_coupling_request(error, cpl);
    REQUIRE(!moist_check_error(error) && !more);
    moist_get_model_gradient(error, model, cpl, response, H2O_NATOMS + 1, gradient);
    if (expect_failure(&error, "gaussian_potential: missing required outputs: dphi_dr, dphi_dxi",
                       "gradient with unanswered derivatives")) goto cleanup;
    for (int k = 0; k < 3 * (H2O_NATOMS + 1); ++k) REQUIRE(gradient[k] == -12345.0);
    /* The next pass retries the request that still misses outputs */
    visits = answer_pass(error, cpl, &host, visited, sizeof visited);
    REQUIRE(visits == 1 && strcmp(visited, "gaussian_potential") == 0);
    moist_get_model_gradient(error, model, cpl, response, H2O_NATOMS - 1, gradient);
    if (expect_failure(&error, "capacity", "short gradient capacity")) goto cleanup;
    for (int k = 0; k < 3 * (H2O_NATOMS + 1); ++k) REQUIRE(gradient[k] == -12345.0);
    moist_get_model_gradient(error, model, cpl, response, H2O_NATOMS + 1, gradient);
    REQUIRE(!moist_check_error(error));
    for (int k = 3 * H2O_NATOMS; k < 3 * (H2O_NATOMS + 1); ++k) REQUIRE(gradient[k] == -12345.0);
    double first_gradient[3 * H2O_NATOMS];
    for (int k = 0; k < 3 * H2O_NATOMS; ++k) first_gradient[k] = gradient[k] + 12345.0;
    moist_get_model_gradient(error, model, cpl, response, H2O_NATOMS + 1, gradient);
    REQUIRE(!moist_check_error(error));
    for (int k = 0; k < 3 * H2O_NATOMS; ++k)
        REQUIRE(fabs(gradient[k] - (-12345.0 + 2.0 * first_gradient[k])) <= 1e-10);
    for (int a = 0; a < H2O_NATOMS; ++a)
        printf("  gradient[%d] = % .10f % .10f % .10f\n", a, first_gradient[3 * a],
               first_gradient[3 * a + 1], first_gradient[3 * a + 2]);
    /* The potential is the nuclear one, so this is the complete gradient and a
     * rigid translation leaves the energy alone. A slip in the [ngrid][3]
     * layout of the dphi_dr answer breaks this by orders of magnitude */
    for (int k = 0; k < 3; ++k) {
        double net = 0.0;
        for (int a = 0; a < H2O_NATOMS; ++a) net += first_gradient[3 * a + k];
        REQUIRE(fabs(net) < 1e-8);
    }
    REQUIRE(response_names(error, response, visited, sizeof visited) == 1 &&
            strcmp(visited, "potential_adjoint") == 0);

    /* A model update ends the walk and the staging of every coupling */
    moist_prepare_model_energy(error, model, other);
    REQUIRE(!moist_check_error(error));
    more = moist_next_coupling_request(error, other);
    REQUIRE(!moist_check_error(error) && more);
    moist_update_model(error, model, mol);
    REQUIRE(!moist_check_error(error));
    moist_answer_coupling_request(error, other, "phi", phi);
    if (expect_failure(&error, no_current_request, "model update invalidates all couplings")) goto cleanup;
    more = moist_next_coupling_request(error, other);
    REQUIRE(!moist_check_error(error) && !more);
    energy = 5.0;
    moist_get_model_energy(error, model, other, &energy);
    if (expect_failure(&error, "get_energy requires a coupling staged by prepare_energy; "
                       "this coupling is not staged", "energy after a model update")) goto cleanup;
    REQUIRE(energy == 5.0);
    borrowed = moist_get_model_cavity(error, model);
    moist_set_isodensity_density(error, borrowed, 0, NULL);
    if (expect_failure(&error, "model-owned", "borrowed cavity density mutation")) goto cleanup;
    result = 0;
cleanup:
    if (moist_check_error(error)) show_error(error);
    free(xyz); free(phi); free(dphi); free(xi); free(dxi); free(w_phi);
    moist_delete_coupling(&cpl); moist_delete_coupling(&other);
    moist_delete_cavity(&borrowed);
    moist_delete_response(&response); moist_delete_model(&model);
    moist_delete_component(&pcm); moist_delete_structure(&mol); moist_delete_error(&error);
    return result;
}

int test_coupling_protocol_fixed_cavity(void)
{
    moist_error error = moist_new_error();
    moist_cavity cav = fixture_iswig(error, NULL, NULL, NULL, NULL, NULL);
    int result = 1;
    if (moist_check_error(error)) {
        show_error(error);
    } else {
        result = run_coupling_protocol("iSwiG cavity, CPCM", cav);
    }
    moist_delete(cav);
    moist_delete(error);
    return result;
}

int test_coupling_protocol_drop_cavity(void)
{
    moist_error error = moist_new_error();
    moist_cavity cav = fixture_drop(error, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
                                             NULL, NULL, NULL, NULL, NULL, NULL, NULL);
    int result = 1;
    if (moist_check_error(error)) {
        show_error(error);
    } else {
        result = run_coupling_protocol("DROP cavity, CPCM", cav);
    }
    moist_delete(cav);
    moist_delete(error);
    return result;
}

/* A density-backed cavity follows the host density: its geometry has a
 * response of its own. The response phase asks for the potential derivatives
 * on top of the energy answer and returns the density weights next to the
 * potential adjoint; the gradient phase then has nothing left to ask */
int test_coupling_protocol_density_cavity(void)
{
    printf("Start test: raw coupling protocol (isodensity DROP cavity, CPCM)\n");
    int result = 1, ngrid = 0, nsph = 0, visits = 0, ncart = 0, nshell = 0;
    moist_error error = moist_new_error();
    moist_structure mol = make_h2o(error);
    moist_pcm_options pcm_options;
    moist_cavity cav = NULL, borrowed = NULL;
    moist_component pcm = NULL;
    moist_model model = NULL;
    moist_coupling cpl = NULL;
    moist_response response = NULL;
    double *dcart = NULL, *xyz = NULL, *xi = NULL, *phi = NULL, *dphi = NULL, *dxi = NULL;
    double *w_phi = NULL, *w_rho = NULL, *w_grad_rho = NULL, *w_hess_rho = NULL;
    double energy = 0.0, gradient[3 * H2O_NATOMS] = {0};
    char visited[128];
    bool more = false;
    /* One normalized s primitive per atom and a diagonal density, as built by
     * test_isodensity_internal_cavity */
    const int shell_atom[3] = {0, 1, 2}, shell_l[3] = {0, 0, 0}, shell_nprim[3] = {1, 1, 1};
    const double alpha = 0.3;
    const double coeff = gto_s_norm(alpha);
    const double exps[3] = {alpha, alpha, alpha}, coeffs[3] = {coeff, coeff, coeff};
    REQUIRE(mol && !moist_check_error(error));
    cav = fixture_isodensity_internal(error, 3, shell_atom, shell_l, shell_nprim, exps, coeffs,
                                      1.0e-3, NULL, NULL, NULL, NULL, NULL, NULL, NULL);
    REQUIRE(cav && !moist_check_error(error));
    moist_get_isodensity_cart_layout(error, cav, &ncart, &nshell, NULL, NULL, NULL, NULL);
    REQUIRE(!moist_check_error(error) && ncart > 0);
    dcart = (double*)calloc((size_t)ncart * (size_t)ncart, sizeof(double));
    REQUIRE(dcart);
    for (int c = 0; c < ncart; ++c) dcart[idx_f2(c, c, ncart)] = 2.0;
    moist_set_isodensity_density(error, cav, ncart, dcart);
    REQUIRE(!moist_check_error(error));
    moist_init_pcm_options(error, &pcm_options, sizeof pcm_options);
    pcm_options.solver = moist_pcm_solver_lu;
    pcm = moist_new_cpcm_component(error, 78.4, &pcm_options);
    model = moist_new_model(error, cav, NULL);
    REQUIRE(pcm && model && !moist_check_error(error));
    moist_add_model_component(error, model, pcm);
    REQUIRE(!moist_check_error(error));
    moist_update_model(error, model, mol);
    REQUIRE(!moist_check_error(error));
    cpl = moist_new_coupling(error, model);
    response = moist_new_response(error);
    REQUIRE(cpl && response && !moist_check_error(error));
    borrowed = moist_get_model_cavity(error, model);
    moist_get_cavity_sizes(error, borrowed, &ngrid, &nsph);
    REQUIRE(!moist_check_error(error) && ngrid > 1);
    xyz = malloc((size_t)3 * ngrid * sizeof(double));
    xi = malloc((size_t)ngrid * sizeof(double));
    phi = malloc((size_t)ngrid * sizeof(double));
    dphi = malloc((size_t)3 * ngrid * sizeof(double));
    dxi = malloc((size_t)ngrid * sizeof(double));
    w_phi = malloc((size_t)ngrid * sizeof(double));
    w_rho = malloc((size_t)ngrid * sizeof(double));
    w_grad_rho = malloc((size_t)3 * ngrid * sizeof(double));
    w_hess_rho = malloc((size_t)9 * ngrid * sizeof(double));
    REQUIRE(xyz && xi && phi && dphi && dxi && w_phi && w_rho && w_grad_rho && w_hess_rho);
    moist_get_cavity_field_real(error, borrowed, "xyz", xyz);
    REQUIRE(!moist_check_error(error));
    moist_get_cavity_field_real(error, borrowed, "xi0", xi);
    REQUIRE(!moist_check_error(error));
    moist_delete_cavity(&borrowed);
    gaussian_nuclear_primitives(ngrid, xyz, xi, phi, dphi, dxi);
    const struct host_outputs host = {ngrid, phi, dphi, dxi, NULL, NULL, NULL, NULL};

    moist_prepare_model_energy(error, model, cpl);
    REQUIRE(!moist_check_error(error));
    visits = answer_pass(error, cpl, &host, visited, sizeof visited);
    REQUIRE(visits == 1 && strcmp(visited, "gaussian_potential") == 0);
    moist_get_model_energy(error, model, cpl, &energy);
    REQUIRE(!moist_check_error(error) && isfinite(energy));
    printf("  energy = %.12f\n", energy);

    /* Unlike a fixed cavity, the response phase walks the potential again */
    moist_prepare_model_response(error, model, cpl);
    REQUIRE(!moist_check_error(error));
    more = moist_next_coupling_request(error, cpl);
    REQUIRE(!moist_check_error(error) && more && request_is(error, cpl, "gaussian_potential"));
    REQUIRE(missing_is(error, cpl, "phi", false) && missing_is(error, cpl, "dphi_dr", true) &&
            missing_is(error, cpl, "dphi_dxi", true));
    more = moist_next_coupling_request(error, cpl);
    REQUIRE(!moist_check_error(error) && !more);
    moist_get_model_response(error, model, cpl, response);
    if (expect_failure(&error, "gaussian_potential: missing required outputs: dphi_dr, dphi_dxi",
                       "response with unanswered derivatives")) goto cleanup;
    visits = answer_pass(error, cpl, &host, visited, sizeof visited);
    REQUIRE(visits == 1 && strcmp(visited, "gaussian_potential") == 0);
    moist_get_model_response(error, model, cpl, response);
    REQUIRE(!moist_check_error(error));
    /* Each item is visited once per pass, in the order the model produced
     * them; the pass ends after the last and rewinds, and the next one starts
     * over */
    for (int pass = 0; pass < 2; ++pass) {
        more = moist_next_response_item(error, response);
        REQUIRE(!moist_check_error(error) && more && item_is(error, response, "potential_adjoint"));
        more = moist_next_response_item(error, response);
        REQUIRE(!moist_check_error(error) && more && item_is(error, response, "density"));
        more = moist_next_response_item(error, response);
        REQUIRE(!moist_check_error(error) && !more);
    }
    /* Refilling the response starts a new walk, even in the middle of one */
    more = moist_next_response_item(error, response);
    REQUIRE(!moist_check_error(error) && more);
    moist_get_model_response(error, model, cpl, response);
    REQUIRE(!moist_check_error(error));
    more = moist_next_response_item(error, response);
    REQUIRE(!moist_check_error(error) && more && item_is(error, response, "potential_adjoint"));
    REQUIRE(!read_response_array(&error, response, "w_phi", ngrid, 1, w_phi));
    more = moist_next_response_item(error, response);
    REQUIRE(!moist_check_error(error) && more && item_is(error, response, "density"));
    REQUIRE(!read_response_array(&error, response, "w_rho", ngrid, 1, w_rho));
    REQUIRE(!read_response_array(&error, response, "w_grad_rho", ngrid, 3, w_grad_rho));
    REQUIRE(!read_response_array(&error, response, "w_hess_rho", ngrid, 9, w_hess_rho));
    moist_get_response_array(error, response, "w_phi", w_phi);
    if (expect_failure(&error, "density has no array 'w_phi'", "array of another item")) goto cleanup;
    more = moist_next_response_item(error, response);
    REQUIRE(!moist_check_error(error) && !more);
    /* The adjoint still reproduces the energy, and the density weights carry
     * the surface response and cannot all vanish */
    double half_q_phi = 0.0, density_norm = 0.0;
    for (int i = 0; i < ngrid; ++i) half_q_phi += 0.5 * w_phi[i] * phi[i];
    for (int i = 0; i < ngrid; ++i) density_norm += w_rho[i] * w_rho[i];
    for (int i = 0; i < 3 * ngrid; ++i) density_norm += w_grad_rho[i] * w_grad_rho[i];
    for (int i = 0; i < 9 * ngrid; ++i) density_norm += w_hess_rho[i] * w_hess_rho[i];
    printf("  q.phi/2 = %.12f, |w_rho, w_grad_rho, w_hess_rho|^2 = %.6e\n", half_q_phi, density_norm);
    REQUIRE(agrees_to(half_q_phi, energy, 1e-10) && density_norm > 0.0);

    /* The response answers cover the gradient phase too. No translation check
     * here: the host's share, the density weights contracted with its density
     * derivatives, is not part of the model gradient */
    moist_prepare_model_gradient(error, model, cpl);
    REQUIRE(!moist_check_error(error));
    more = moist_next_coupling_request(error, cpl);
    REQUIRE(!moist_check_error(error) && !more);
    moist_get_model_gradient(error, model, cpl, response, H2O_NATOMS, gradient);
    REQUIRE(!moist_check_error(error));
    for (int k = 0; k < 3 * H2O_NATOMS; ++k) REQUIRE(isfinite(gradient[k]));
    result = 0;
cleanup:
    if (moist_check_error(error)) show_error(error);
    free(dcart); free(xyz); free(xi); free(phi); free(dphi); free(dxi);
    free(w_phi); free(w_rho); free(w_grad_rho); free(w_hess_rho);
    moist_delete_coupling(&cpl);
    moist_delete_cavity(&borrowed);
    moist_delete_response(&response);
    moist_delete_model(&model);
    moist_delete_component(&pcm);
    moist_delete_cavity(&cav);
    moist_delete_structure(&mol);
    moist_delete_error(&error);
    return result;
}

/* Two requests in declaration order: GOSTSHYP adds its Gaussian moments after
 * the CPCM potential. A pass visits each request at most once, a pass left
 * early resumes where it was left, every prepare restarts the walk, and a later
 * pass retries only the requests that still miss an output */
int test_coupling_protocol_gostshyp(void)
{
    printf("Start test: raw coupling protocol (iSwiG cavity, CPCM + GOSTSHYP)\n");
    int result = 1, ngrid = 0, nsph = 0, visits = 0;
    const int padding = 5;
    moist_error error = moist_new_error();
    moist_structure mol = make_h2o(error);
    moist_cavity cav = fixture_iswig(error, NULL, NULL, NULL, NULL, NULL);
    moist_component pcm = moist_new_cpcm_component(error, 78.4, NULL);
    moist_component gostshyp = moist_new_gostshyp_component(error, 1.0e-4);
    moist_model model = NULL;
    moist_coupling cpl = NULL;
    moist_response response = NULL;
    moist_cavity borrowed = NULL;
    double *xyz = NULL, *xi = NULL, *phi = NULL, *dphi = NULL, *dxi = NULL;
    double *gt = NULL, *pt = NULL, *mt = NULL, *rt = NULL, *width = NULL;
    double *w_overlap = NULL, *w_normal_deriv = NULL;
    double energy = 0.0, gradient[3 * H2O_NATOMS] = {0};
    char visited[128];
    bool more = false;
    REQUIRE(mol && cav && pcm && gostshyp && !moist_check_error(error));
    model = moist_new_model(error, cav, NULL);
    REQUIRE(model && !moist_check_error(error));
    moist_add_model_component(error, model, pcm);
    REQUIRE(!moist_check_error(error));
    moist_add_model_component(error, model, gostshyp);
    REQUIRE(!moist_check_error(error));
    moist_update_model(error, model, mol);
    REQUIRE(!moist_check_error(error));
    cpl = moist_new_coupling(error, model);
    response = moist_new_response(error);
    REQUIRE(cpl && response && !moist_check_error(error));
    borrowed = moist_get_model_cavity(error, model);
    moist_get_cavity_sizes(error, borrowed, &ngrid, &nsph);
    REQUIRE(!moist_check_error(error) && ngrid > 1);
    xyz = malloc((size_t)3 * ngrid * sizeof(double));
    xi = malloc((size_t)ngrid * sizeof(double));
    phi = malloc((size_t)ngrid * sizeof(double));
    dphi = malloc((size_t)3 * ngrid * sizeof(double));
    dxi = malloc((size_t)ngrid * sizeof(double));
    gt = malloc((size_t)ngrid * sizeof(double));
    pt = malloc((size_t)3 * ngrid * sizeof(double));
    mt = malloc((size_t)9 * ngrid * sizeof(double));
    rt = malloc((size_t)3 * ngrid * sizeof(double));
    width = malloc((size_t)(ngrid + padding) * sizeof(double));
    w_overlap = malloc((size_t)ngrid * sizeof(double));
    w_normal_deriv = malloc((size_t)ngrid * sizeof(double));
    REQUIRE(xyz && xi && phi && dphi && dxi && gt && pt && mt && rt && width &&
            w_overlap && w_normal_deriv);
    moist_get_cavity_field_real(error, borrowed, "xyz", xyz);
    REQUIRE(!moist_check_error(error));
    moist_get_cavity_field_real(error, borrowed, "xi0", xi);
    REQUIRE(!moist_check_error(error));
    moist_delete_cavity(&borrowed);
    gaussian_nuclear_primitives(ngrid, xyz, xi, phi, dphi, dxi);
    /* Stand-ins for the density traces a QM host integrates */
    for (int i = 0; i < ngrid; ++i) {
        gt[i] = 1.0e-3;
        for (int k = 0; k < 3; ++k) {
            pt[3 * i + k] = 1.0e-4 * (k + 1);
            rt[3 * i + k] = -2.0e-4 * (k + 1);
            for (int l = 0; l < 3; ++l) mt[9 * i + 3 * k + l] = k == l ? 1.0e-3 : 0.0;
        }
    }
    const struct host_outputs host = {ngrid, phi, dphi, dxi, gt, pt, mt, rt};

    /* Declaration order: the potential, then the moments. A prepare in the
     * middle of a pass restarts the walk at the first request */
    moist_prepare_model_energy(error, model, cpl);
    REQUIRE(!moist_check_error(error));
    more = moist_next_coupling_request(error, cpl);
    REQUIRE(more && request_is(error, cpl, "gaussian_potential"));
    more = moist_next_coupling_request(error, cpl);
    REQUIRE(more && request_is(error, cpl, "gaussian_moments"));
    moist_prepare_model_energy(error, model, cpl);
    REQUIRE(!moist_check_error(error));
    /* Leave the pass at the potential with phi unanswered... */
    while (moist_next_coupling_request(error, cpl)) {
        visits++;
        REQUIRE(request_is(error, cpl, "gaussian_potential") && missing_is(error, cpl, "phi", true));
        break;
    }
    REQUIRE(!moist_check_error(error) && visits == 1);
    /* ...and the next call resumes at the moments instead of rewinding */
    more = moist_next_coupling_request(error, cpl);
    REQUIRE(more && request_is(error, cpl, "gaussian_moments"));
    REQUIRE(missing_is(error, cpl, "gt", true) && missing_is(error, cpl, "pt", true));
    REQUIRE(missing_is(error, cpl, "mt", false) && missing_is(error, cpl, "rt", false));
    REQUIRE(missing_is(error, cpl, "phi", false));
    /* The exponents are an input of this request: read them, never recompute */
    const size_t width_bytes = (size_t)(ngrid + padding) * sizeof(double);
    fill_sentinel(width, width_bytes);
    moist_get_coupling_request_width(error, cpl, NULL);
    if (expect_failure(&error, "Null array pointer provided for 'width'", "NULL width buffer")) goto cleanup;
    moist_get_coupling_request_width(error, cpl, width);
    REQUIRE(!moist_check_error(error));
    REQUIRE(sentinel_intact(width + ngrid, (size_t)padding * sizeof(double)));
    for (int i = 0; i < ngrid; ++i) REQUIRE(isfinite(width[i]) && width[i] > 0.0);
    /* gt is taken; pt is refused for a non-finite value and stays missing */
    moist_answer_coupling_request(error, cpl, "gt", gt);
    REQUIRE(!moist_check_error(error));
    const double saved_pt = pt[0];
    pt[0] = NAN;
    moist_answer_coupling_request(error, cpl, "pt", pt);
    pt[0] = saved_pt;
    if (expect_failure(&error, "gaussian_moments: pt must be finite",
                       "non-finite moment answer")) goto cleanup;
    REQUIRE(missing_is(error, cpl, "gt", false) && missing_is(error, cpl, "pt", true));
    /* Each request once per pass: the skipped potential is not revisited */
    more = moist_next_coupling_request(error, cpl);
    REQUIRE(!moist_check_error(error) && !more);
    energy = 3.0;
    moist_get_model_energy(error, model, cpl, &energy);
    if (expect_failure(&error, "gaussian_potential: missing required outputs: phi",
                       "energy with a skipped request")) goto cleanup;
    REQUIRE(energy == 3.0);
    /* The next pass retries both requests; the one after has nothing left */
    visits = answer_pass(error, cpl, &host, visited, sizeof visited);
    REQUIRE(visits == 2 && strcmp(visited, "gaussian_potential,gaussian_moments") == 0);
    visits = answer_pass(error, cpl, &host, visited, sizeof visited);
    REQUIRE(visits == 0);
    energy = 0.0;
    moist_get_model_energy(error, model, cpl, &energy);
    REQUIRE(!moist_check_error(error) && isfinite(energy));
    printf("  energy = %.12f\n", energy);

    /* On a fixed cavity the response phase needs nothing new; the response
     * carries the amplitudes next to the potential adjoint */
    moist_prepare_model_response(error, model, cpl);
    REQUIRE(!moist_check_error(error));
    visits = answer_pass(error, cpl, &host, visited, sizeof visited);
    REQUIRE(visits == 0);
    moist_get_model_response(error, model, cpl, response);
    REQUIRE(!moist_check_error(error));
    REQUIRE(response_names(error, response, visited, sizeof visited) == 2 &&
            strcmp(visited, "potential_adjoint,gostshyp_amplitude") == 0);
    REQUIRE(walk_to_item(error, response, "gostshyp_amplitude") == 1);
    REQUIRE(!read_response_array(&error, response, "w_overlap", ngrid, 1, w_overlap));
    REQUIRE(!read_response_array(&error, response, "w_normal_deriv", ngrid, 1, w_normal_deriv));
    double amplitude_norm = 0.0;
    for (int i = 0; i < ngrid; ++i)
        amplitude_norm += w_overlap[i] * w_overlap[i] + w_normal_deriv[i] * w_normal_deriv[i];
    printf("  |w_overlap, w_normal_deriv|^2 = %.6e\n", amplitude_norm);
    REQUIRE(amplitude_norm > 0.0);
    moist_get_response_array(error, response, "w_phi", w_overlap);
    if (expect_failure(&error, "gostshyp_amplitude has no array 'w_phi'",
                       "array of another item")) goto cleanup;

    /* The gradient phase asks both requests for more. A refused output is
     * retried alone: the next pass visits its request only */
    moist_prepare_model_gradient(error, model, cpl);
    REQUIRE(!moist_check_error(error));
    more = moist_next_coupling_request(error, cpl);
    REQUIRE(more && request_is(error, cpl, "gaussian_potential"));
    REQUIRE(missing_is(error, cpl, "phi", false) && missing_is(error, cpl, "dphi_dr", true) &&
            missing_is(error, cpl, "dphi_dxi", true));
    moist_answer_coupling_request(error, cpl, "dphi_dr", dphi);
    REQUIRE(!moist_check_error(error));
    moist_answer_coupling_request(error, cpl, "dphi_dxi", dxi);
    REQUIRE(!moist_check_error(error));
    more = moist_next_coupling_request(error, cpl);
    REQUIRE(more && request_is(error, cpl, "gaussian_moments"));
    REQUIRE(missing_is(error, cpl, "gt", false) && missing_is(error, cpl, "pt", false) &&
            missing_is(error, cpl, "mt", true) && missing_is(error, cpl, "rt", true));
    moist_answer_coupling_request(error, cpl, "mt", mt);
    REQUIRE(!moist_check_error(error));
    const double saved = rt[0];
    rt[0] = NAN;
    moist_answer_coupling_request(error, cpl, "rt", rt);
    rt[0] = saved;
    if (expect_failure(&error, "gaussian_moments: rt must be finite", "non-finite moment")) goto cleanup;
    more = moist_next_coupling_request(error, cpl);
    REQUIRE(!moist_check_error(error) && !more);
    moist_get_model_gradient(error, model, cpl, response, H2O_NATOMS, gradient);
    if (expect_failure(&error, "gaussian_moments: missing required outputs: rt",
                       "gradient with a refused output")) goto cleanup;
    visits = answer_pass(error, cpl, &host, visited, sizeof visited);
    REQUIRE(visits == 1 && strcmp(visited, "gaussian_moments") == 0);
    moist_get_model_gradient(error, model, cpl, response, H2O_NATOMS, gradient);
    REQUIRE(!moist_check_error(error));
    for (int a = 0; a < H2O_NATOMS; ++a)
        printf("  gradient[%d] = % .10f % .10f % .10f\n", a, gradient[3 * a],
               gradient[3 * a + 1], gradient[3 * a + 2]);
    for (int k = 0; k < 3 * H2O_NATOMS; ++k) REQUIRE(isfinite(gradient[k]));
    REQUIRE(response_names(error, response, visited, sizeof visited) == 2 &&
            strcmp(visited, "potential_adjoint,gostshyp_amplitude") == 0);
    result = 0;
cleanup:
    if (moist_check_error(error)) show_error(error);
    free(xyz); free(xi); free(phi); free(dphi); free(dxi);
    free(gt); free(pt); free(mt); free(rt); free(width);
    free(w_overlap); free(w_normal_deriv);
    moist_delete_coupling(&cpl);
    moist_delete_cavity(&borrowed);
    moist_delete_response(&response);
    moist_delete_model(&model);
    moist_delete_component(&gostshyp);
    moist_delete_component(&pcm);
    moist_delete_cavity(&cav);
    moist_delete_structure(&mol);
    moist_delete_error(&error);
    return result;
}

#undef REQUIRE

/* -1 selects version text; nonnegative values select one of the four banners */
static void query_text(moist_error error, int style, char *buffer,
                       size_t capacity, size_t *length)
{
    if (style == -1) moist_get_version_string(error, buffer, capacity, length);
    else moist_get_banner(error, (moist_banner)style, buffer, capacity, length);
}

int test_string_queries(void)
{
    moist_error error = moist_new_error();
    char *full = NULL;
    char small[8];
    int failed = 1;
    for (int style = -1; style <= moist_banner_build; ++style) {
        size_t needed = 0, length = 123;
        query_text(error, style, NULL, 0, &needed);
        if (moist_check_error(error) || needed < sizeof small) goto cleanup;
        full = malloc(needed + 1);
        if (!full) goto cleanup;
        query_text(error, style, full, needed + 1, &length);
        if (moist_check_error(error) || length != needed || strlen(full) != needed) goto cleanup;
        memset(small, 'X', sizeof small);
        query_text(error, style, small, sizeof small, &length);
        if (moist_check_error(error) || length != needed || small[7] != '\0' ||
            memcmp(small, full, 7)) goto cleanup;
        memset(small, 'X', sizeof small);
        query_text(error, style, small, 1, &length);
        if (moist_check_error(error) || length != needed || small[0] != '\0' || small[1] != 'X') goto cleanup;
        small[0] = 'X';
        query_text(error, style, small, 0, &length);
        if (moist_check_error(error) || length != needed || small[0] != 'X') goto cleanup;
        length = 123;
        small[0] = 'X';
        query_text(error, style, small, SIZE_MAX, &length);
        if (!moist_check_error(error) || small[0] != '\0' || length != 123) goto cleanup;
        query_text(error, style, NULL, 1, &length);
        if (!moist_check_error(error) || length != 123) goto cleanup;
        query_text(error, style, small, sizeof small, NULL);
        if (!moist_check_error(error) || small[0] != '\0') goto cleanup;
        small[0] = 'X';
        query_text(NULL, style, small, sizeof small, &length);
        if (small[0] != '\0' || length != 123) goto cleanup;
        free(full);
        full = NULL;
    }
    size_t length = 123;
    small[0] = 'X';
    moist_get_banner(error, (moist_banner)99, small, sizeof small, &length);
    if (!moist_check_error(error) || length != 123 || small[0] != '\0') goto cleanup;
    failed = 0;
cleanup:
    if (failed) show_error(error);
    free(full);
    moist_delete(error);
    return failed;
}

/* Historical 1.0 wire layout. Never append fields to this fixture */
struct drop_options_1_0 {
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
    int reserved0;
};

struct cfc_options_1_0 {
    size_t struct_size;
    double a1, a2, c;
    int m, reserved0;
};

struct pcm_options_1_0 {
    size_t struct_size;
    moist_pcm_solver solver;
    int solver_maxiter;
    double solver_tol;
};

/* Simulate an old binary, even after the installed header gains fields */
int test_options_1_0_prefix(void)
{
    moist_error error = moist_new_error();
    moist_lsf lsf = NULL;
    moist_cavity cavity = NULL;
    moist_component pcm = NULL;
    struct { struct cfc_options_1_0 prefix; double appended; } cfc;
    struct { struct pcm_options_1_0 prefix; double appended; } pcm_options;
    struct { struct drop_options_1_0 prefix; double appended; } future;
    int failed = 1;
    memset(&future, 0x5a, sizeof future);
    moist_init_drop_options(error, (moist_drop_options *)&future.prefix, sizeof future.prefix);
    if (moist_check_error(error) || future.prefix.struct_size != sizeof future.prefix) goto cleanup;
    for (size_t i = sizeof future.prefix; i < sizeof future; ++i)
        if (((unsigned char *)&future)[i] != 0x5a) goto cleanup;
    lsf = moist_new_svdw_lsf(error, NULL);
    if (!lsf || moist_check_error(error)) goto cleanup;
    cavity = moist_new_drop_cavity(error, lsf, NULL, (const moist_drop_options *)&future.prefix);
    if (!cavity || moist_check_error(error)) goto cleanup;
    /* Reserved inputs stay ignored, including residue from pre-reservation callers */
    moist_delete(cavity);
    future.prefix.reserved0 = 123;
    cavity = moist_new_drop_cavity(error, lsf, NULL, (const moist_drop_options *)&future.prefix);
    if (!cavity || moist_check_error(error)) goto cleanup;

    memset(&cfc, 0x5a, sizeof cfc);
    moist_init_cfc_options(error, (moist_cfc_options *)&cfc.prefix, sizeof cfc.prefix);
    if (moist_check_error(error) || cfc.prefix.struct_size != sizeof cfc.prefix) goto cleanup;
    for (size_t i = sizeof cfc.prefix; i < sizeof cfc; ++i)
        if (((unsigned char *)&cfc)[i] != 0x5a) goto cleanup;
    moist_delete(lsf);
    cfc.prefix.reserved0 = 123;
    lsf = moist_new_cfc_lsf(error, (const moist_cfc_options *)&cfc.prefix);
    if (!lsf || moist_check_error(error)) goto cleanup;

    memset(&pcm_options, 0x5a, sizeof pcm_options);
    moist_init_pcm_options(error, (moist_pcm_options *)&pcm_options.prefix, sizeof pcm_options.prefix);
    if (moist_check_error(error) || pcm_options.prefix.struct_size != sizeof pcm_options.prefix) goto cleanup;
    for (size_t i = sizeof pcm_options.prefix; i < sizeof pcm_options; ++i)
        if (((unsigned char *)&pcm_options)[i] != 0x5a) goto cleanup;
    pcm = moist_new_cpcm_component(error, 80.0, (const moist_pcm_options *)&pcm_options.prefix);
    if (!pcm || moist_check_error(error)) goto cleanup;
    failed = 0;
cleanup:
    moist_delete(pcm);
    if (failed) show_error(error);
    moist_delete(cavity);
    moist_delete(lsf);
    moist_delete(error);
    return failed;
}

/* Check every byte outside meaningful fields, including reserved members */
static int padding_is_zero(const void *value, const unsigned char *fields, size_t size)
{
    const unsigned char *bytes = value;
    for (size_t i = 0; i < size; ++i)
        if (!fields[i] && bytes[i] != 0) return 0;
    return 1;
}

int test_options_padding(void)
{
    moist_error error = moist_new_error();
    int failed = 1;
#define MARK_FIELD(kind, member) \
    memset(fields + offsetof(moist_##kind##_options, member), 1, sizeof value.member)
    {
        moist_drop_options value;
        unsigned char fields[sizeof value] = {0};
        MARK_FIELD(drop, struct_size);
        MARK_FIELD(drop, nleb);
        MARK_FIELD(drop, debug);
        MARK_FIELD(drop, verbosity);
        MARK_FIELD(drop, do_fine);
        MARK_FIELD(drop, tolerance);
        MARK_FIELD(drop, proj_maxiter);
        MARK_FIELD(drop, proj_level);
        MARK_FIELD(drop, branch_weight_s);
        MARK_FIELD(drop, rho_grid_h);
        MARK_FIELD(drop, wleb_prune_level);
        for (int pass = 0; pass < 2; ++pass) {
            memset(&value, pass ? 0x5a : 0xa5, sizeof value);
            moist_init_drop_options(error, &value, sizeof value);
            if (moist_check_error(error) || !padding_is_zero(&value, fields, sizeof value)) goto cleanup;
        }
    }
    {
        moist_iswig_options value;
        unsigned char fields[sizeof value] = {0};
        MARK_FIELD(iswig, struct_size);
        MARK_FIELD(iswig, nleb);
        MARK_FIELD(iswig, debug);
        MARK_FIELD(iswig, verbosity);
        MARK_FIELD(iswig, cut_a);
        MARK_FIELD(iswig, cut_f);
        for (int pass = 0; pass < 2; ++pass) {
            memset(&value, pass ? 0x5a : 0xa5, sizeof value);
            moist_init_iswig_options(error, &value, sizeof value);
            if (moist_check_error(error) || !padding_is_zero(&value, fields, sizeof value)) goto cleanup;
        }
    }
    {
        moist_svdw_options value;
        unsigned char fields[sizeof value] = {0};
        MARK_FIELD(svdw, struct_size);
        MARK_FIELD(svdw, blend_k);
        MARK_FIELD(svdw, blend_1b);
        MARK_FIELD(svdw, blend_2b);
        MARK_FIELD(svdw, blend_3b);
        for (int pass = 0; pass < 2; ++pass) {
            memset(&value, pass ? 0x5a : 0xa5, sizeof value);
            moist_init_svdw_options(error, &value, sizeof value);
            if (moist_check_error(error) || !padding_is_zero(&value, fields, sizeof value)) goto cleanup;
        }
    }
    {
        moist_cfc_options value;
        unsigned char fields[sizeof value] = {0};
        MARK_FIELD(cfc, struct_size);
        MARK_FIELD(cfc, a1);
        MARK_FIELD(cfc, a2);
        MARK_FIELD(cfc, c);
        MARK_FIELD(cfc, m);
        for (int pass = 0; pass < 2; ++pass) {
            memset(&value, pass ? 0x5a : 0xa5, sizeof value);
            moist_init_cfc_options(error, &value, sizeof value);
            if (moist_check_error(error) || !padding_is_zero(&value, fields, sizeof value)) goto cleanup;
        }
    }
    {
        moist_isodensity_options value;
        unsigned char fields[sizeof value] = {0};
        MARK_FIELD(isodensity, struct_size);
        MARK_FIELD(isodensity, rho_iso);
        MARK_FIELD(isodensity, scale);
        for (int pass = 0; pass < 2; ++pass) {
            memset(&value, pass ? 0x5a : 0xa5, sizeof value);
            moist_init_isodensity_options(error, &value, sizeof value);
            if (moist_check_error(error) || !padding_is_zero(&value, fields, sizeof value)) goto cleanup;
        }
    }
    {
        moist_model_options value;
        unsigned char fields[sizeof value] = {0};
        MARK_FIELD(model, struct_size);
        MARK_FIELD(model, debug);
        MARK_FIELD(model, verbosity);
        for (int pass = 0; pass < 2; ++pass) {
            memset(&value, pass ? 0x5a : 0xa5, sizeof value);
            moist_init_model_options(error, &value, sizeof value);
            if (moist_check_error(error) || !padding_is_zero(&value, fields, sizeof value)) goto cleanup;
        }
    }
    {
        moist_pcm_options value;
        unsigned char fields[sizeof value] = {0};
        MARK_FIELD(pcm, struct_size);
        MARK_FIELD(pcm, solver);
        MARK_FIELD(pcm, solver_maxiter);
        MARK_FIELD(pcm, solver_tol);
        for (int pass = 0; pass < 2; ++pass) {
            memset(&value, pass ? 0x5a : 0xa5, sizeof value);
            moist_init_pcm_options(error, &value, sizeof value);
            if (moist_check_error(error) || !padding_is_zero(&value, fields, sizeof value)) goto cleanup;
        }
    }
#undef MARK_FIELD
    failed = 0;
cleanup:
    if (failed) show_error(error);
    moist_delete(error);
    return failed;
}

/* Iterative-solver controls reach the native validator through the options struct */
int test_pcm_solver_options(void)
{
    moist_error error = moist_new_error();
    moist_component pcm = NULL;
    moist_pcm_options options;
    int failed = 1;
    moist_init_pcm_options(error, &options, sizeof options);
    if (moist_check_error(error)) goto cleanup;
    /* Defaults are the native ones: a positive threshold and a positive cap */
    if (!(options.solver_tol > 0.0) || options.solver_maxiter < 1) goto cleanup;
    options.solver = moist_pcm_solver_iterative;
    options.solver_tol = 1.0e-8;
    options.solver_maxiter = 200;
    pcm = moist_new_cpcm_component(error, 78.4, &options);
    if (!pcm || moist_check_error(error)) goto cleanup;
    moist_delete(pcm);
    pcm = moist_new_cosmo_component(error, 78.4, &options);
    if (!pcm || moist_check_error(error)) goto cleanup;
    moist_delete(pcm);
    /* A zero cap and a NaN threshold are rejected, for either variant */
    options.solver_maxiter = 0;
    pcm = moist_new_cpcm_component(error, 78.4, &options);
    if (pcm || !moist_check_error(error)) goto cleanup;
    options.solver_maxiter = 200;
    options.solver_tol = NAN;
    pcm = moist_new_cosmo_component(error, 78.4, &options);
    if (pcm || !moist_check_error(error)) goto cleanup;
    failed = 0;
cleanup:
    if (failed) show_error(error);
    moist_delete(pcm);
    moist_delete(error);
    return failed;
}

/* ABI contracts independent of a particular numerical model */
int test_v1_contract(void)
{
    moist_error error = moist_new_error();
    moist_lsf lsf = NULL;
    moist_cavity cavity = NULL;
    moist_structure mol = NULL;
    int failed = 1;
    char buffer[64];
    const int capacity = sizeof buffer;
    const int one = 1;
    double energy = 123.0;
    bool flag = true;
    int ngrid = 123, nsph = 456;
    moist_drop_options options;
    struct { moist_drop_options options; unsigned char tail[16]; } future;

    memset(buffer, 'X', sizeof buffer);
    moist_get_error(error, buffer, &capacity);
    if (buffer[0] != '\0' || buffer[1] != 'X') goto cleanup;
    moist_get_error(NULL, buffer, &capacity);
    if (!memchr(buffer, '\0', sizeof buffer) || !strstr(buffer, "Invalid")) goto cleanup;
    memset(buffer, 'X', sizeof buffer);
    moist_get_error(NULL, buffer, &one);
    if (buffer[0] != '\0' || buffer[1] != 'X') goto cleanup;
    memset(buffer, 'X', sizeof buffer);
    size_t length = 123;
    moist_get_banner(NULL, moist_banner_full, buffer, sizeof buffer, &length);
    if (buffer[0] != '\0' || buffer[1] != 'X') goto cleanup;
    moist_get_version_string(error, buffer, sizeof buffer, &length);
    if (strcmp(buffer, "1.0.0-alpha.1")) goto cleanup;
    moist_get_version_string(error, buffer, 1, &length);
    if (buffer[0] != '\0') goto cleanup;

    /* Every C layout must match its frozen Fortran minimum, including padding */
#define CHECK_OPTIONS_LAYOUT(kind) do { \
    moist_##kind##_options value; \
    moist_init_##kind##_options(error, &value, sizeof value); \
    if (moist_check_error(error) || value.struct_size != sizeof value) goto cleanup; \
    memset(&value, 0x5a, sizeof value); \
    moist_init_##kind##_options(error, &value, sizeof value - 1); \
    if (!moist_check_error(error)) goto cleanup; \
    for (size_t i = 0; i < sizeof value; ++i) \
        if (((unsigned char *)&value)[i] != 0x5a) goto cleanup; \
} while (0)
    CHECK_OPTIONS_LAYOUT(drop);
    CHECK_OPTIONS_LAYOUT(iswig);
    CHECK_OPTIONS_LAYOUT(svdw);
    CHECK_OPTIONS_LAYOUT(cfc);
    CHECK_OPTIONS_LAYOUT(isodensity);
    CHECK_OPTIONS_LAYOUT(model);
    CHECK_OPTIONS_LAYOUT(pcm);
#undef CHECK_OPTIONS_LAYOUT

    memset(&options, 0x5a, sizeof options);
    moist_init_drop_options(error, &options, sizeof options - 1);
    if (!moist_check_error(error)) goto cleanup;
    for (size_t i = 0; i < sizeof options; ++i)
        if (((unsigned char *)&options)[i] != 0x5a) goto cleanup;
    memset(&future, 0x5a, sizeof future);
    moist_init_drop_options(error, &future.options, sizeof future);
    if (moist_check_error(error) || future.options.struct_size != sizeof future) goto cleanup;
    for (size_t i = 0; i < sizeof future.tail; ++i)
        if (future.tail[i] != 0x5a) goto cleanup;
    options = future.options;
    options.struct_size = sizeof options;
    if (options.nleb != 194 || options.tolerance != 1e-10) goto cleanup;

    cavity = moist_new_drop_cavity(error, NULL, NULL, NULL);
    if (cavity || !moist_check_error(error)) goto cleanup;
    lsf = moist_new_svdw_lsf(error, NULL);
    if (!lsf || moist_check_error(error)) goto cleanup;
    cavity = moist_new_drop_cavity(error, lsf, NULL, &options);
    if (!cavity || moist_check_error(error)) goto cleanup;
    moist_delete(lsf);  /* The cavity owns its LSF copy */
    mol = make_h2o(error);
    moist_update_cavity(error, cavity, mol);
    if (moist_check_error(error)) goto cleanup;
    moist_get_cavity_sizes(error, cavity, NULL, &nsph);
    if (!moist_check_error(error) || nsph != 456) goto cleanup;
    moist_get_cavity_sizes(error, NULL, &ngrid, &nsph);
    if (!moist_check_error(error) || ngrid != 123 || nsph != 456) goto cleanup;
    moist_get_model_energy(error, NULL, NULL, &energy);
    if (!moist_check_error(error) || energy != 123.0) goto cleanup;
    moist_get_model_energy(error, NULL, NULL, NULL);
    if (!moist_check_error(error)) goto cleanup;
    /* The coupling walk ends on any failure, even without an error handle */
    if (moist_next_coupling_request(NULL, NULL)) goto cleanup;
    if (moist_next_coupling_request(error, NULL) || !moist_check_error(error)) goto cleanup;
    moist_get_coupling_request_missing(error, NULL, "phi", &flag);
    if (!moist_check_error(error) || !flag) goto cleanup;
    moist_get_coupling_request_missing(error, NULL, "phi", NULL);
    if (!moist_check_error(error)) goto cleanup;
    moist_get_coupling_request_name(error, NULL, NULL);
    if (!moist_check_error(error)) goto cleanup;
    /* So does the response walk */
    if (moist_next_response_item(NULL, NULL)) goto cleanup;
    if (moist_next_response_item(error, NULL) || !moist_check_error(error)) goto cleanup;
    moist_get_response_item_name(error, NULL, NULL);
    if (!moist_check_error(error)) goto cleanup;
    moist_get_cavity_field_real(error, cavity, "xyz", NULL);
    if (!moist_check_error(error)) goto cleanup;
    moist_get_cavity_field_info(error, cavity, 0, NULL, &ngrid, &nsph, NULL, NULL);
    if (!moist_check_error(error) || ngrid != 123 || nsph != 456) goto cleanup;
    moist_get_cavity_sizes(error, cavity, &ngrid, &nsph);
    if (moist_check_error(error) || ngrid <= 0 || nsph != H2O_NATOMS) goto cleanup;
    failed = 0;
cleanup:
    if (failed) show_error(error);
    moist_delete(lsf);
    moist_delete(cavity);
    moist_delete(mol);
    moist_delete(error);
    return failed;
}

static const struct {
    const char* name;
    int (*fn)(void);
} test_registry[] = {
    {"v1_contract", test_v1_contract},
    {"options_1_0_prefix", test_options_1_0_prefix},
    {"options_padding", test_options_padding},
    {"pcm_solver_options", test_pcm_solver_options},
    {"string_queries", test_string_queries},
    {"header_and_version",                   test_header_and_version},
    {"version",                              test_version},
    {"uninitialized_error",                  test_uninitialized_error},
    {"error_origins",                        test_error_origins},
    {"model_density_rejection",              test_model_density_rejection},
    {"null_handle",                          test_null_handle},
    {"delete_resets_handle",                 test_delete_resets_handle},
    {"drop_cavity",                          test_drop_cavity},
    {"custom_radii",                         test_custom_radii},
    {"h2o_cavity",                           test_h2o_cavity},
    {"cavity_gradient",                      test_cavity_gradient},
    {"isodensity_internal_cavity",           test_isodensity_internal_cavity},
    {"isodensity_callback_cavity",           test_isodensity_callback_cavity},
    {"isodensity_callback_third_derivative", test_isodensity_callback_third_derivative},
    {"isodensity_callback_failure",          test_isodensity_callback_failure},
    {"update_drop_cavity_keeps_params",      test_update_drop_cavity_keeps_params},
    {"cavity_fields",                         test_cavity_fields},
    {"capacity_validation",                  test_capacity_validation},
    {"coupling_protocol_fixed_cavity",        test_coupling_protocol_fixed_cavity},
    {"coupling_protocol_drop_cavity",        test_coupling_protocol_drop_cavity},
    {"coupling_protocol_density_cavity",     test_coupling_protocol_density_cavity},
    {"coupling_protocol_gostshyp",           test_coupling_protocol_gostshyp},
};

int main(void)
{
    const size_t ntests = sizeof(test_registry) / sizeof(test_registry[0]);
    int failed[sizeof(test_registry) / sizeof(test_registry[0])];
    size_t nfailed = 0;

    for (size_t i = 0; i < ntests; i++) {
        failed[i] = test_registry[i].fn() != 0;
        if (failed[i]) nfailed++;
    }

    printf("\n=== %zu of %zu tests passed ===\n", ntests - nfailed, ntests);
    if (nfailed > 0) {
        printf("failed:\n");
        for (size_t i = 0; i < ntests; i++) {
            if (failed[i]) printf("  - %s\n", test_registry[i].name);
        }
    }

    return nfailed == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}
