#include <limits.h>
#include <math.h>
#include <stdatomic.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "moist.h"

/* Report the message behind a failed handle.
 *
 * moist_get_error leaves the buffer untouched -- not even a terminating NUL --
 * when the handle carries no error, so a caller that prints unconditionally
 * would print stack garbage. Both halves of the guard matter: the early return
 * keeps a clean handle silent, and the initializer keeps the buffer printable
 * if the getter ever declines to write for some other reason. */
static inline void
show_error(moist_error error)
{
    if (!moist_check_error(error)) return;

    char message[512] = {0};
    const int message_len = (int)sizeof(message);
    moist_get_error(error, message, &message_len);
    printf("[Message] %s\n", message);
}

/* Caller-side "abort on failure" policy.
 *
 * The library deliberately offers no entry point that terminates the process
 * -- a library linked into a host program must never do that. Deciding to die
 * is the driver's business, and it is these few lines: check the handle, pull
 * the message out, report, exit. Each test below instead reports and returns a
 * status so that main() can run the remaining tests, which is what a real host
 * should do too. */
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

/* Normalization coefficient of a single s primitive, (2a/pi)^(3/4).
 *
 * M_PI is not in the C standard; glibc hides it once __STRICT_ANSI__ is set,
 * which is exactly what this project's c_std=c11 does. acos(-1.0) is portable. */
static inline double gto_s_norm(double alpha)
{
    return pow(2.0 * alpha / acos(-1.0), 0.75);
}

/* Agreement between two independently summed results, measured relative to the
 * magnitude involved.*/
static int agrees_to(double a, double b, double rel_tol)
{
    double scale = fabs(a) > fabs(b) ? fabs(a) : fabs(b);
    if (scale < 1.0) {
        scale = 1.0;
    }
    return fabs(a - b) <= rel_tol * scale;
}

/* Fortran column-major flattening helpers */
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

/* The H2O geometry nearly every test below builds on, in Bohr.
 *
 * Static storage rather than a literal repeated in each test: the isodensity
 * callbacks hold a pointer to the centers for the lifetime of a build, and one
 * shared definition keeps the callback and internal backends from drifting
 * apart. A real host would of course pass its own arrays. */
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
    /* A clean handle must be a no-op for the caller-side abort helper. */
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
     * moist.h is C++-clean and C++ has no variable-length arrays, so an example
     * meant to be copied must not use them. The A-matrix is the reason it also
     * matters here and not only for portability: it is ngrid*ngrid doubles,
     * which is nothing to put on a stack once the solute grows.
     *
     * Vector arrays are Fortran (3,ngrid): read component k of point i as
     * xyz[3*i + k] (see the layout convention in moist.h). */
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

    cav = moist_new_drop_cavity_with_radii(error, radii_model,
                                           NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
                                           NULL, NULL, NULL, NULL, NULL, NULL);
    if (moist_check_error(error)) {
        show_error(error);
        goto cleanup;
    }

    // Now build the cavity with structure
    moist_update_cavity(error, cav, mol);
    if (moist_check_error(error)) {
        show_error(error);
        goto cleanup;
    }

    // First get generic cavity sizes to allocate arrays correctly
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

    // Get generic cavity results (Tier 1 - works for all cavity types)
    moist_get_cavity_results(error, cav, ngrid, nsph,
                             &area, &volume, &ngrid, &nsph,
                             xyz, a, owner, converged, vradii, asph);
    if (moist_check_error(error)) {
        show_error(error);
        goto cleanup;
    }

    // Get DROP-specific data (Tier 2 - only works for DROP cavities)
    moist_get_drop_specific(error, cav, ngrid, &nmax,
                            normal, wleb, r_iI0,
                            f, rho);
    if (moist_check_error(error)) {
        show_error(error);
        goto cleanup;
    }

    // Assemble A-matrix and get xi values (DROP-specific)
    moist_assemble_amat(error, cav, ngrid, amat0, xi);
    if (moist_check_error(error)) {
        show_error(error);
        goto cleanup;
    }

    // Validate results before cleanup
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

    // Cleanup and verify handles are reset to NULL (using generic delete)
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

    // Element-specific custom radii.
    int z_list[2] = {1, 8};
    double r_elem[2] = {2.0, 3.2};
    moist_set_custom_radii_elements(error, radii_model, 2, z_list, r_elem);
    if (moist_check_error(error)) {
        show_error(error);
        goto cleanup;
    }

    cav = moist_new_drop_cavity_with_radii(error, radii_model,
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

    /* The second round rebuilds at a different radius set, so the grid size can
     * change: release the first round's buffers before resizing. */
    free(xyz); free(a); free(owner); free(converged); free(vradii); free(asph);
    xyz = NULL; a = NULL; owner = NULL; converged = NULL; vradii = NULL; asph = NULL;
    moist_delete_cavity(&cav);

    // Atom-specific custom radii.
    double r_atom[3] = {1.7, 2.1, 2.4};
    moist_set_custom_radii_atoms(error, radii_model, natoms, r_atom);
    if (moist_check_error(error)) {
        show_error(error);
        goto cleanup;
    }

    cav = moist_new_drop_cavity_with_radii(error, radii_model,
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
    moist_print_header(6);  // Print to stdout (Fortran unit 6)
    moist_print_version(6);
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

    // Allocate arrays dynamically (C++ compatible, no VLAs)
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
    cav = moist_new_drop_cavity(error, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
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

    moist_get_drop_specific(error, cav, ngrid, &nmax, normal, wleb, r_iI0,
                            f, rho);
    if (moist_check_error(error)) {
        show_error(error);
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

    /* Single-exit cleanup.
     *
     * Every owned pointer is declared here and starts NULL, so the one ladder
     * under `cleanup:` is correct no matter how early the test leaves -- adding
     * an allocation means touching two places, not seventeen. `result` starts
     * at failure; only the path that runs to the end clears it. */
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
    /* Outputs of the surface-to-LSF contraction. They get their own buffers on
     * purpose: w_lsf2 is Fortran (3,3,ngrid), so reusing the ngrid*ngrid
     * A-matrix for it would only fit when ngrid >= 9 and would destroy the
     * matrix validated just above. Capacity discipline is what this file is
     * meant to demonstrate. */
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
    cav = moist_new_drop_cavity(error, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
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

    // Now get A-matrix with gradient
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

    // Test API: contract_nuc_elec_qefield_rA
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

    // Both the original and extended surface-contraction ABI must remain callable.
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
    /* Symmetry and finiteness are both satisfied by an all-zero result, so the
     * magnitude has to be part of the assertion, not a printed remark. */
    if (contract_ok && contract_norm < 1e-20) {
        printf("  Error: contract_amat1_q1q2_rA norm is essentially zero\n");
        contract_ok = 0;
    }

    // Zero-input check for contract_nuc_elec_qefield_rA: should return zero gradient
    for (int i = 0; i < ngrid; i++) {
        surface_q[i] = 0.0;
        for (int iaxis = 0; iaxis < 3; iaxis++) {
            qefield[idx_f2(iaxis, i, 3)] = 0.0;
        }
    }
    for (int iatom = 0; iatom < nsph; iatom++) {
        za[iatom] = (iatom < natoms) ? (double)numbers[iatom] : 0.0;
    }

    moist_contract_nuc_elec_qefield_rA(error, cav, surface_q, qefield, za, grad_ne_zero);
    if (moist_check_error(error)) {
        show_error(error);
        goto cleanup;
    }

    int nuc_elec_ok = 1;
    const double zero_tol = 1e-14;
    for (int iatom = 0; iatom < nsph && nuc_elec_ok; iatom++) {
        for (int iaxis = 0; iaxis < 3; iaxis++) {
            if (fabs(grad_ne_zero[idx_f2(iaxis, iatom, 3)]) > zero_tol) {
                printf("  Error: contract_nuc_elec_qefield_rA zero-input check failed\n");
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

    moist_contract_nuc_elec_qefield_rA(error, cav, surface_q, qefield, za, grad_ne);
    if (moist_check_error(error)) {
        show_error(error);
        goto cleanup;
    }

    moist_contract_nuc_elec_qefield_rA(error, cav, surface_q_scaled, qefield_scaled, za, grad_ne_scaled);
    if (moist_check_error(error)) {
        show_error(error);
        goto cleanup;
    }

    const double homo_tol = 1e-12;
    for (int iatom = 0; iatom < nsph && nuc_elec_ok; iatom++) {
        for (int iaxis = 0; iaxis < 3; iaxis++) {
            size_t idx = idx_f2(iaxis, iatom, 3);
            if (!isfinite(grad_ne[idx]) || !isfinite(grad_ne_scaled[idx])) {
                printf("  Error: contract_nuc_elec_qefield_rA produced non-finite values\n");
                nuc_elec_ok = 0;
                break;
            }
            if (!agrees_to(grad_ne_scaled[idx], 2.0 * grad_ne[idx], homo_tol)) {
                printf("  Error: contract_nuc_elec_qefield_rA homogeneity check failed\n");
                nuc_elec_ok = 0;
                break;
            }
        }
    }

    /* Zero-input, finiteness and homogeneity are all satisfied by a routine
     * that writes nothing but zeros. Pin the magnitude too: with nonzero
     * surface charges, a nonzero field and real nuclear charges, the
     * nuclear-electron gradient cannot vanish. */
    double nuc_elec_norm = 0.0;
    for (int i = 0; i < 3 * nsph; i++) {
        nuc_elec_norm += grad_ne[i] * grad_ne[i];
    }
    if (nuc_elec_ok && nuc_elec_norm < 1e-20) {
        printf("  Error: contract_nuc_elec_qefield_rA gradient is essentially zero\n");
        nuc_elec_ok = 0;
    }

    // Validate results - check that gradients are computed and non-trivial
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
 * Uses a minimal one-primitive s shell per atom so the level set
 * S = scale * (rho_iso - rho) is an analytic sum of spherical Gaussians. The
 * test exercises the full internal-backend entry-point chain: constructor ->
 * two-pass layout query -> density install -> build -> numbering / anchor
 * gradients / tolerance getter.
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

    cav = moist_new_drop_cavity_isodensity_internal(
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
     * Uses a scratch handle so the expected failure does not poison `error`. */
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
    moist_get_drop_numbering(error, cav, ngrid, numbering);
    if (moist_check_error(error)) {
        show_error(error);
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
 * Isodensity level set supplied through the C callback ABI.
 *
 * This is the reference implementation of the V_0_6 callback contract: `hess`
 * and `third` may be NULL, and NULL means "do not compute this derivative",
 * not merely "do not store it". A callback that dereferences them
 * unconditionally -- which the V_0_5 contract allowed -- crashes here.
 *
 * The level set mirrors the diagonal-density case built by
 * test_isodensity_internal_cavity() exactly, so the two backends must agree:
 *   rho(r) = sum_A c exp(-a |r - R_A|^2),  S(r) = rho_iso - rho(r)
 * with c = 2 coeff^2 and a = 2 alpha.
 * ------------------------------------------------------------------------- */

/* moist_get_cavity_results writes every array unconditionally, so the caller
 * must size and allocate them even when only the scalars are wanted. */
static int cavity_area_volume(moist_error error, moist_cavity cav,
                              double* area, double* volume, int* ngrid, int* nsph)
{
    moist_get_cavity_sizes(error, cav, ngrid, nsph);
    if (moist_check_error(error)) return 1;
    /* Failing without setting the error handle would leave the caller with
     * nothing to report, so say what went wrong here. */
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
     * The projection evaluates its grid points in an OpenMP loop, so every
     * counter the callback touches is shared across threads. */
    atomic_int calls;
    atomic_int calls_with_hess;
    atomic_int calls_with_third;
};

static int iso_gaussian_callback(void* context, const double* point,
                                 double* value, double* grad,
                                 double* hess, double* third)
{
    struct iso_callback_ctx* ctx = (struct iso_callback_ctx*)context;
    const int want_hess = (hess != NULL);
    const int want_third = (third != NULL);

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

    /* DROP sign convention: interior negative, exterior positive */
    *value = ctx->rho_iso - rho;
    for (int i = 0; i < 3; i++) grad[i] = -drho[i];
    if (want_hess) for (int i = 0; i < 9; i++) hess[i] = -d2rho[i];
    if (want_third) for (int i = 0; i < 27; i++) third[i] = -d3rho[i];

    return 0;   /* success; any nonzero value would abort the cavity build */
}

int test_isodensity_callback_cavity(void)
{
    printf("Start test: callback isodensity cavity\n");
#ifndef moist_API_SUFFIX__V_0_6
#  error "this example implements the V_0_6 isodensity callback contract"
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

    cav = moist_new_drop_cavity_isodensity_callback(
        error, iso_gaussian_callback, &ctx, NULL, NULL, NULL, NULL, NULL, NULL, NULL);
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

    /* Snapshot the counters once: the build is over, so they are stable, and
     * the report and the assertions below must agree on the same numbers. */
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

    /* The contract is only exercised if both phases actually occurred: some
     * calls must have skipped the higher derivatives and some must not have. */
    if (n_hess <= 0 || n_hess >= n_calls) {
        printf("  FAIL: the value+gradient-only phase was never exercised\n");
        result = 1;
    }

    /* The third-derivative channel is only requested at max_deriv >= 3, which a
     * plain cavity build never reaches -- see the dedicated coverage in
     * test_isodensity_callback_third_derivative(). What must hold here is the
     * nesting of the phases: third derivatives are never requested without a
     * Hessian alongside them. */
    if (n_third > n_hess) {
        printf("  FAIL: third derivatives were requested without a Hessian\n");
        result = 1;
    }

    /* Cross-check against the internal backend built from the equivalent basis:
     * same level set, so the same surface. */
    const int nshell_in = 3;
    int shell_atom[3] = {0, 1, 2};
    int shell_l[3] = {0, 0, 0};
    int shell_nprim[3] = {1, 1, 1};
    double exps[3], coeffs[3];
    for (int i = 0; i < nshell_in; i++) {
        exps[i] = alpha;
        coeffs[i] = coeff;
    }

    ref = moist_new_drop_cavity_isodensity_internal(
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

/* The third-derivative channel of the V_0_6 callback contract.
 *
 * A plain cavity build never reaches it: the projection runs the level-set
 * model at max_deriv 1 and 2, so `third` is always NULL and the d3rho branch of
 * iso_gaussian_callback -- which moist.h advertises as the reference
 * implementation of this contract -- would otherwise never execute.
 *
 * The surface-to-LSF contraction is the path that raises the model to
 * max_deriv = 3 (see src/moist/cavity/drop/derivatives/potential.f90), so
 * running it against a callback-backed cavity is what puts the remaining half
 * of the reference implementation under test. */
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

    cav = moist_new_drop_cavity_isodensity_callback(
        error, iso_gaussian_callback, &ctx, NULL, NULL, NULL, NULL, NULL, NULL, NULL);
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

    /* Nonzero surface adjoints, so the contraction has something to propagate
     * and cannot short-circuit before touching the level-set model. */
    for (int i = 0; i < ngrid; i++) {
        w_xi[i] = 1.0e-3 * (double)(i + 1);
        w_f[i] = 2.0e-3;
        for (int iaxis = 0; iaxis < 3; iaxis++) {
            w_xyz[idx_f2(iaxis, i, 3)] = 1.0e-4 * (double)(iaxis + 1);
        }
    }

    /* Count only what the contraction itself asks for: the build above already
     * ran the value/gradient and Hessian phases. */
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

    /* max_deriv = 3 implies max_deriv >= 2, so every call that got a third
     * derivative must have got a Hessian with it. */
    if (n_third > n_hess) {
        printf("  FAIL: third derivatives were requested without a Hessian\n");
        result = 1;
    }

    /* The weights the third-derivative channel feeds must be usable numbers. */
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
 * The failure is reported on the 51st evaluation rather than the first: the
 * projection runs its grid points in an OpenMP loop, so an abort raised from
 * the middle of the loop is the case that actually exercises the shared
 * failure flag. A first-evaluation failure would pass even if the parallel
 * handling were wrong. */
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
                                double* value, double* grad,
                                double* hess, double* third)
{
    struct iso_fail_ctx* ctx = (struct iso_fail_ctx*)context;

    const int n = atomic_fetch_add(&ctx->calls, 1) + 1;
    if (n > ctx->fail_after) return ctx->status;

    return iso_gaussian_callback(&ctx->base, point, value, grad, hess, third);
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

    cav = moist_new_drop_cavity_isodensity_callback(
        error, iso_failing_callback, &ctx, NULL, NULL, NULL, NULL, NULL, NULL, NULL);
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
        /* Reading a message does not clear the handle, so swap in a fresh one
         * before the recovery build below. */
        moist_delete_error(&error);
        error = moist_new_error();
    }

    /* The callback must not have been re-entered indefinitely: once it reports
     * failure moist stops calling it, so the count stays close to the point of
     * failure rather than covering the whole grid. */
    const int total_calls = atomic_load(&ctx.calls);
    printf("  callback invocations: %d (failure requested after %d)\n",
           total_calls, ctx.fail_after);
    if (total_calls <= ctx.fail_after) {
        printf("  FAIL: the failing branch was never reached\n");
        result = 1;
    }

    /* The build must also be repeatable: a callback that stops failing has to
     * produce a normal cavity, i.e. the failure latch is per-build state. */
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

    /* Upper bound on the aborted build.
     *
     * The lower bound above only proves the failing branch was entered; an
     * implementation that reported the first error and then kept invoking the
     * callback across the whole grid would satisfy it too, while breaking the
     * latch contract. The recovery build just above is that same grid evaluated
     * to completion, which makes it the natural yardstick: stopping at
     * evaluation 51 must cost a small fraction of a full traversal, with only
     * the calls already in flight on the other threads added on top. */
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

/* The deprecated moist_update_drop_cavity rebuilds the cavity from scratch to
 * honour a new Lebedev order. Rebuilding must not silently roll the cavity's
 * configured numerical settings back to their compiled defaults, which is what
 * re-running the constructor does unless they are carried across explicitly. */
int test_update_drop_cavity_keeps_params(void)
{
    printf("Start test: update_drop_cavity preserves configured parameters\n");

    int result = 1;

    moist_error error = moist_new_error();
    moist_structure mol = NULL;
    moist_cavity cav = NULL;

    const int natoms = H2O_NATOMS;

    int nleb = 26;
    const double tolerance_in = 1.0e-8;   /* deliberately not the default 1e-10 */

    mol = make_h2o(error);
    if (moist_check_error(error)) {
        show_error(error);
        goto cleanup;
    }

    cav = moist_new_drop_cavity(error, &nleb, NULL, NULL, NULL, NULL,
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

    /* Rebuild at a different Lebedev order through the legacy entry point */
    int nleb_new = 50;
    moist_update_drop_cavity(error, cav, mol, &nleb_new);
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

    /* The rebuilt cavity must still be usable. The models are re-sourced from
     * detached copies, so a use-after-free here would surface as garbage sizes
     * or a crash. */
    int ngrid = 0, nmax = 0, nsph = 0;
    moist_get_drop_sizes(error, cav, &ngrid, &nmax, &nsph);
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

/* The array-writing getters must refuse a capacity smaller than the cavity
 * needs, and must accept one that is larger. The undersized calls below are the
 * heap corruption this contract exists to prevent: without the capacity check
 * they write the cavity's full ngrid elements into a buffer the caller declared
 * smaller.
 *
 * The buffers here are therefore allocated at *full* size and handed a smaller
 * *logical* capacity. That is the only arrangement in which the test can both
 * survive a missing check and detect it: a genuinely short allocation would
 * corrupt the heap before any assertion ran, and a single trailing canary would
 * miss a write that lands anywhere else. Every byte is stamped with a sentinel
 * beforehand and the whole region is verified afterwards, so a rejected call
 * has to have written nothing at all, anywhere. */

/* Not a value any getter would plausibly write: neither a small integer, nor a
 * zero double, nor a valid bool. */
#define CAPACITY_SENTINEL 0xA5

static void fill_sentinel(void* p, size_t nbytes)
{
    memset(p, CAPACITY_SENTINEL, nbytes);
}

/* Non-zero if every byte of the region still carries the sentinel. */
static int sentinel_intact(const void* p, size_t nbytes)
{
    const unsigned char* b = (const unsigned char*)p;
    for (size_t i = 0; i < nbytes; i++) {
        if (b[i] != (unsigned char)CAPACITY_SENTINEL) return 0;
    }
    return 1;
}

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
    double* normal = NULL;
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
    cav = moist_new_drop_cavity(error, NULL, NULL, NULL, NULL, NULL,
                                NULL, NULL, NULL, NULL,
                                NULL, NULL, NULL, NULL, NULL);
    moist_update_cavity(error, cav, mol);

    int ngrid = 0, nsph = 0;
    moist_get_cavity_sizes(error, cav, &ngrid, &nsph);
    if (moist_check_error(error) || ngrid <= 1 || nsph != natoms) {
        show_error(error);
        goto cleanup;
    }

    /* Physical extents: what an unchecked getter would write. */
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
    normal = (double*)malloc(n_xyz);
    numbering = (int*)malloc(n_owner);
    amat0 = (double*)malloc(n_amat);
    xi = (double*)malloc(n_a);
    gaussian_f = (double*)malloc(n_a);

    if (!xyz || !a || !owner || !converged || !vradii || !asph || !wleb ||
        !r_iI0 || !f || !rho || !normal || !numbering || !amat0 || !xi ||
        !gaussian_f) {
        printf("  Error: memory allocation failed\n");
        goto cleanup;
    }

    /* Logical capacity offered to the library: one short of the truth. */
    const int short_ngrid = ngrid - 1;

    double area = 0.0, volume = 0.0;
    int out_ngrid = 0, out_nsph = 0, nmax = 0;

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
        /* A short sphere capacity has to be caught just the same. The grid
         * capacity is honest here so that only the sphere check can fire. */
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

    if (result == 0) {
        fill_sentinel(normal, n_xyz); fill_sentinel(wleb, n_a);
        fill_sentinel(r_iI0, n_a); fill_sentinel(f, n_a); fill_sentinel(rho, n_a);

        moist_get_drop_specific(error, cav, short_ngrid, &nmax,
                                normal, wleb, r_iI0, f, rho);
        if (!moist_check_error(error)) {
            printf("  FAIL: get_drop_specific accepted a short capacity\n");
            result = 1;
        } else {
            if (!sentinel_intact(normal, n_xyz) || !sentinel_intact(wleb, n_a) ||
                !sentinel_intact(r_iI0, n_a) || !sentinel_intact(f, n_a) ||
                !sentinel_intact(rho, n_a)) {
                printf("  FAIL: get_drop_specific wrote into a rejected buffer\n");
                result = 1;
            }
            moist_delete_error(&error);
            error = moist_new_error();
        }
    }

    if (result == 0) {
        fill_sentinel(numbering, n_owner);

        moist_get_drop_numbering(error, cav, short_ngrid, numbering);
        if (!moist_check_error(error)) {
            printf("  FAIL: get_drop_numbering accepted a short capacity\n");
            result = 1;
        } else {
            if (!sentinel_intact(numbering, n_owner)) {
                printf("  FAIL: get_drop_numbering wrote into a rejected buffer\n");
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

    /* An oversized buffer is legal: a host may allocate once for the largest
     * geometry it expects and reuse it. The data lands in the leading ngrid
     * entries, with the capacity setting the stride of the vector array. */
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
            if (moist_check_error(error)) {
                show_error(error);
                printf("  FAIL: get_cavity_results rejected an oversized capacity\n");
                result = 1;
            } else if (out_ngrid != ngrid || out_nsph != nsph ||
                       !(area > 0.0) || !(big_a[0] > 0.0) ||
                       !(big_radii[nsph - 1] > 0.0) ||
                       big_a[ngrid] != 0.0 || big_radii[nsph] != 0.0) {
                printf("  FAIL: oversized capacity produced unexpected results\n");
                result = 1;
            }
        }
    }

    if (result == 0) {
        printf("  Array capacity validation behaved as specified!\n");
    }

cleanup:
    free(xyz); free(a); free(owner); free(converged); free(vradii); free(asph);
    free(wleb); free(r_iI0); free(f); free(rho); free(normal); free(numbering);
    free(amat0); free(xi); free(gaussian_f);
    free(big_xyz); free(big_a); free(big_owner); free(big_conv);
    free(big_radii); free(big_asph);

    moist_delete_cavity(&cav);
    moist_delete_structure(&mol);
    moist_delete_error(&error);
    return result;
}

/* Registry rather than a running sum: the exit status alone tells you how many
 * tests failed but not which, and the "Start test:" lines are interleaved with
 * everything the tests print. */
static const struct {
    const char* name;
    int (*fn)(void);
} test_registry[] = {
    {"header_and_version",                 test_header_and_version},
    {"version",                            test_version},
    {"uninitialized_error",                test_uninitialized_error},
    {"null_handle",                        test_null_handle},
    {"delete_resets_handle",               test_delete_resets_handle},
    {"drop_cavity",                        test_drop_cavity},
    {"custom_radii",                       test_custom_radii},
    {"h2o_cavity",                         test_h2o_cavity},
    {"cavity_gradient",                    test_cavity_gradient},
    {"isodensity_internal_cavity",         test_isodensity_internal_cavity},
    {"isodensity_callback_cavity",         test_isodensity_callback_cavity},
    {"isodensity_callback_third_derivative", test_isodensity_callback_third_derivative},
    {"isodensity_callback_failure",        test_isodensity_callback_failure},
    {"update_drop_cavity_keeps_params",    test_update_drop_cavity_keeps_params},
    {"capacity_validation",                test_capacity_validation},
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
