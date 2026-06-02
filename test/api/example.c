#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "moist.h"

static inline void
show_error(moist_error error)
{
    char message[512];
    const int message_len = (int)sizeof(message);
    moist_get_error(error, message, &message_len);
    printf("[Message] %s\n", message);
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

// int test_gems_model(void)
// {
//     printf("Start test: GEMS solvation model\n");
//     moist_error error = moist_new_error();

//     const int natoms = 3;
//     int numbers[3] = {8, 1, 1};  // O, H, H
//     double positions[9] = {
//         0.0000,  0.0000,  0.1173,
//         0.0000,  1.4309, -0.9370,
//         0.0000, -1.4309, -0.9370
//     };

//     moist_structure mol = moist_new_structure(error, natoms, numbers, positions, NULL, NULL);
//     if (moist_check_error(error)) {
//         show_error(error);
//         moist_delete_error(&error);
//         return 1;
//     }

//     const bool debug = false;
//     const int verbosity = 0;
//     const bool read_parameters = false;
//     const char* parameter_file = NULL;

//     moist_model model = moist_new_gems_solvation_model(error, "water",
//                                                        &debug, &verbosity,
//                                                        &read_parameters,
//                                                        parameter_file);
//     if (moist_check_error(error)) {
//         show_error(error);
//         moist_delete_structure(&mol);
//         moist_delete_error(&error);
//         return 1;
//     }

//     moist_update_solvation_model(error, model, mol);
//     if (moist_check_error(error)) {
//         show_error(error);
//         moist_delete_solvation_model(&model);
//         moist_delete_structure(&mol);
//         moist_delete_error(&error);
//         return 1;
//     }

//     double energy = 0.0;
//     moist_get_solvation_model_energy(error, model, &energy);
//     if (moist_check_error(error)) {
//         show_error(error);
//         moist_delete_solvation_model(&model);
//         moist_delete_structure(&mol);
//         moist_delete_error(&error);
//         return 1;
//     }

//     moist_delete_solvation_model(&model);
//     moist_delete_structure(&mol);
//     moist_delete_error(&error);

//     if (model != NULL || mol != NULL || error != NULL) {
//         printf("Error: handles not reset to NULL after deletion\n");
//         return 1;
//     }

//     return 0;
// }

// int test_gems_unknown_solvent_reports_error(void)
// {
//     printf("Start test: GEMS unknown solvent reports error\n");
//     moist_error error = moist_new_error();
//     const bool debug = false;
//     const int verbosity = 0;
//     const bool read_parameters = false;

//     moist_model model = moist_new_gems_solvation_model(error, "definitely-not-a-solvent",
//                                                        &debug, &verbosity,
//                                                        &read_parameters, NULL);
//     if (model != NULL) {
//         printf("Error: model handle should be NULL for unknown solvent\n");
//         moist_delete_solvation_model(&model);
//         moist_delete_error(&error);
//         return 1;
//     }
//     if (!moist_check_error(error)) {
//         printf("Error: expected moist_error for unknown solvent\n");
//         moist_delete_error(&error);
//         return 1;
//     }

//     moist_delete_error(&error);
//     return error == NULL ? 0 : 1;
// }

int test_drop_cavity(void)
{
    printf("Start test: DROP cavity build\n");
    moist_error error = moist_new_error();
    moist_radii radii_model = NULL;

    const int natoms = 1;
    int numbers[1] = {1};
    double positions[3] = {0.0, 0.0, 0.0};
    moist_structure mol =
        moist_new_structure(error, natoms, numbers, positions, NULL, NULL);

    if (moist_check_error(error)) {
        show_error(error);
        moist_delete_error(&error);
        return 1;
    }

    // Create explicit radii model and cavity handle (does not build yet)
    radii_model = moist_new_cpcm_radii(error);
    if (moist_check_error(error)) {
        show_error(error);
        moist_delete_structure(&mol);
        moist_delete_error(&error);
        return 1;
    }

    moist_cavity cav = moist_new_drop_cavity_with_radii(error, radii_model,
                                                       NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL);

    if (moist_check_error(error)) {
        show_error(error);
        moist_delete_radii(&radii_model);
        moist_delete_structure(&mol);
        moist_delete_error(&error);
        return 1;
    }

    // Now build the cavity with structure
    moist_update_cavity(error, cav, mol);

    if (moist_check_error(error)) {
        show_error(error);
        moist_delete_cavity(&cav);
        moist_delete_radii(&radii_model);
        moist_delete_structure(&mol);
        moist_delete_error(&error);
        return 1;
    }

    // First get generic cavity sizes to allocate arrays correctly
    double area = 0.0, volume = 0.0;
    int ngrid = 0, nsph = 0;
    moist_get_cavity_sizes(error, cav, &ngrid, &nsph);
    
    if (moist_check_error(error)) {
        show_error(error);
        moist_delete_cavity(&cav);
        moist_delete_radii(&radii_model);
        moist_delete_structure(&mol);
        moist_delete_error(&error);
        return 1;
    }

    // Use Variable Length Arrays (VLAs) - C99 feature
    // Now we have correct sizes from the cavity!
    double xyz[3 * ngrid];
    double a[ngrid];
    int owner[ngrid];
    bool converged[ngrid];
    double vradii[nsph];
    double asph[nsph];
    
    // DROP-specific arrays
    double normal[3 * ngrid];
    double wleb[ngrid], r_iI0[ngrid], xi[ngrid];
    double f[ngrid], rho[ngrid];
    int nmax = 0;

    // Get generic cavity results (Tier 1 - works for all cavity types)
    moist_get_cavity_results(error, cav,
                             &area, &volume, &ngrid, &nsph,
                             xyz, a, owner, converged, vradii, asph);

    int status = moist_check_error(error);
    
    if (status != 0) {
        show_error(error);
        moist_delete_cavity(&cav);
        moist_delete_radii(&radii_model);
        moist_delete_structure(&mol);
        moist_delete_error(&error);
        return 1;
    }
    
    // Get DROP-specific data (Tier 2 - only works for DROP cavities)
    moist_get_drop_specific(error, cav, &nmax,
                           normal, wleb, r_iI0,
                           f, rho);
    
    status = moist_check_error(error);
    
    if (status != 0) {
        show_error(error);
        moist_delete_cavity(&cav);
        moist_delete_radii(&radii_model);
        moist_delete_structure(&mol);
        moist_delete_error(&error);
        return 1;
    }
    
    // Assemble A-matrix and get xi values (DROP-specific)
    double amat0[ngrid * ngrid];
    moist_assemble_amat(error, cav, &ngrid, amat0, xi);
    
    status = moist_check_error(error);
    
    // Validate results before cleanup
    int result = status == 0 && area > 0.0 && ngrid > 0 &&
                   nmax >= ngrid && wleb[0] > 0.0 && a[0] > 0.0 && 
                   r_iI0[0] > 0.0 && xi[0] > 0.0 &&
                   f[0] > 0.0 && rho[0] >= 0.0 && 
                   converged[0] == true && asph[0] > 0.0 &&
                   vradii[0] > 0.0
               ? 0
               : 1;

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
    moist_error error = moist_new_error();
    moist_radii radii_model = NULL;

    const int natoms = 3;
    int numbers[3] = {8, 1, 1};  // O, H, H
    double positions[9] = {
        0.0000,  0.0000,  0.1173,
        0.0000,  1.4309, -0.9370,
        0.0000, -1.4309, -0.9370
    };

    moist_structure mol = moist_new_structure(error, natoms, numbers, positions, NULL, NULL);
    if (moist_check_error(error)) {
        show_error(error);
        moist_delete_error(&error);
        return 1;
    }

    radii_model = moist_new_custom_radii(error);
    if (moist_check_error(error)) {
        show_error(error);
        moist_delete_structure(&mol);
        moist_delete_error(&error);
        return 1;
    }

    // Element-specific custom radii.
    int z_list[2] = {1, 8};
    double r_elem[2] = {2.0, 3.2};
    moist_set_custom_radii_elements(error, radii_model, 2, z_list, r_elem);
    if (moist_check_error(error)) {
        show_error(error);
        moist_delete_radii(&radii_model);
        moist_delete_structure(&mol);
        moist_delete_error(&error);
        return 1;
    }

    moist_cavity cav = moist_new_drop_cavity_with_radii(error, radii_model,
                                                       NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL);
    if (moist_check_error(error)) {
        show_error(error);
        moist_delete_radii(&radii_model);
        moist_delete_structure(&mol);
        moist_delete_error(&error);
        return 1;
    }

    moist_update_cavity(error, cav, mol);
    if (moist_check_error(error)) {
        show_error(error);
        moist_delete_cavity(&cav);
        moist_delete_radii(&radii_model);
        moist_delete_structure(&mol);
        moist_delete_error(&error);
        return 1;
    }

    int ngrid = 0, nsph = 0;
    moist_get_cavity_sizes(error, cav, &ngrid, &nsph);
    if (moist_check_error(error) || nsph != natoms) {
        if (moist_check_error(error)) show_error(error);
        moist_delete_cavity(&cav);
        moist_delete_radii(&radii_model);
        moist_delete_structure(&mol);
        moist_delete_error(&error);
        return 1;
    }

    double* xyz = (double*)malloc(3 * ngrid * sizeof(double));
    double* a = (double*)malloc(ngrid * sizeof(double));
    int* owner = (int*)malloc(ngrid * sizeof(int));
    bool* converged = (bool*)malloc(ngrid * sizeof(bool));
    double* vradii = (double*)malloc(nsph * sizeof(double));
    double* asph = (double*)malloc(nsph * sizeof(double));
    double area = 0.0, volume = 0.0;
    if (!xyz || !a || !owner || !converged || !vradii || !asph) {
        free(xyz); free(a); free(owner); free(converged); free(vradii); free(asph);
        moist_delete_cavity(&cav);
        moist_delete_radii(&radii_model);
        moist_delete_structure(&mol);
        moist_delete_error(&error);
        return 1;
    }

    moist_get_cavity_results(error, cav, &area, &volume, &ngrid, &nsph,
                             xyz, a, owner, converged, vradii, asph);
    if (moist_check_error(error)) {
        show_error(error);
        free(xyz); free(a); free(owner); free(converged); free(vradii); free(asph);
        moist_delete_cavity(&cav);
        moist_delete_radii(&radii_model);
        moist_delete_structure(&mol);
        moist_delete_error(&error);
        return 1;
    }

    int element_ok = fabs(vradii[0] - r_elem[1]) < 1e-12 &&
                     fabs(vradii[1] - r_elem[0]) < 1e-12 &&
                     fabs(vradii[2] - r_elem[0]) < 1e-12;
    free(xyz); free(a); free(owner); free(converged); free(vradii); free(asph);
    moist_delete_cavity(&cav);

    // Atom-specific custom radii.
    double r_atom[3] = {1.7, 2.1, 2.4};
    moist_set_custom_radii_atoms(error, radii_model, natoms, r_atom);
    if (moist_check_error(error)) {
        show_error(error);
        moist_delete_radii(&radii_model);
        moist_delete_structure(&mol);
        moist_delete_error(&error);
        return 1;
    }

    cav = moist_new_drop_cavity_with_radii(error, radii_model,
                                          NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL);
    if (moist_check_error(error)) {
        show_error(error);
        moist_delete_radii(&radii_model);
        moist_delete_structure(&mol);
        moist_delete_error(&error);
        return 1;
    }

    moist_update_cavity(error, cav, mol);
    if (moist_check_error(error)) {
        show_error(error);
        moist_delete_cavity(&cav);
        moist_delete_radii(&radii_model);
        moist_delete_structure(&mol);
        moist_delete_error(&error);
        return 1;
    }

    moist_get_cavity_sizes(error, cav, &ngrid, &nsph);
    if (moist_check_error(error) || nsph != natoms) {
        if (moist_check_error(error)) show_error(error);
        moist_delete_cavity(&cav);
        moist_delete_radii(&radii_model);
        moist_delete_structure(&mol);
        moist_delete_error(&error);
        return 1;
    }

    xyz = (double*)malloc(3 * ngrid * sizeof(double));
    a = (double*)malloc(ngrid * sizeof(double));
    owner = (int*)malloc(ngrid * sizeof(int));
    converged = (bool*)malloc(ngrid * sizeof(bool));
    vradii = (double*)malloc(nsph * sizeof(double));
    asph = (double*)malloc(nsph * sizeof(double));
    if (!xyz || !a || !owner || !converged || !vradii || !asph) {
        free(xyz); free(a); free(owner); free(converged); free(vradii); free(asph);
        moist_delete_cavity(&cav);
        moist_delete_radii(&radii_model);
        moist_delete_structure(&mol);
        moist_delete_error(&error);
        return 1;
    }

    moist_get_cavity_results(error, cav, &area, &volume, &ngrid, &nsph,
                             xyz, a, owner, converged, vradii, asph);
    if (moist_check_error(error)) {
        show_error(error);
        free(xyz); free(a); free(owner); free(converged); free(vradii); free(asph);
        moist_delete_cavity(&cav);
        moist_delete_radii(&radii_model);
        moist_delete_structure(&mol);
        moist_delete_error(&error);
        return 1;
    }

    int atom_ok = fabs(vradii[0] - r_atom[0]) < 1e-12 &&
                  fabs(vradii[1] - r_atom[1]) < 1e-12 &&
                  fabs(vradii[2] - r_atom[2]) < 1e-12;
    free(xyz); free(a); free(owner); free(converged); free(vradii); free(asph);

    moist_delete_cavity(&cav);
    moist_delete_radii(&radii_model);
    moist_delete_structure(&mol);
    moist_delete_error(&error);

    return (element_ok && atom_ok) ? 0 : 1;
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
    moist_error error = moist_new_error();

    // Define H2O molecule (coordinates in Bohr)
    // Standard equilibrium geometry: r(OH) = 0.9572 Å = 1.8088 Bohr, angle = 104.52°
    const int natoms = 3;
    int numbers[3] = {8, 1, 1};  // O, H, H
    double positions[9] = {
        0.0000,  0.0000,  0.1173,   // O atom
        0.0000,  1.4309, -0.9370,   // H atom 1
        0.0000, -1.4309, -0.9370    // H atom 2
    };

    moist_structure mol = moist_new_structure(error, natoms, numbers, positions, NULL, NULL);
    if (moist_check_error(error)) {
        show_error(error);
        moist_delete_error(&error);
        return 1;
    }

    // Create DROP cavity handle (does not build yet)
    moist_cavity cav = moist_new_drop_cavity(error, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL);
    if (moist_check_error(error)) {
        show_error(error);
        moist_delete_structure(&mol);
        moist_delete_error(&error);
        return 1;
    }

    // Build cavity with structure
    printf("  Building cavity...\n");
    moist_update_cavity(error, cav, mol);
    if (moist_check_error(error)) {
        show_error(error);
        moist_delete_cavity(&cav);
        moist_delete_structure(&mol);
        moist_delete_error(&error);
        return 1;
    }
    printf("  Cavity built successfully\n");

    // Get cavity sizes for array allocation
    int ngrid = 0, nsph = 0;
    moist_get_cavity_sizes(error, cav, &ngrid, &nsph);
    if (moist_check_error(error)) {
        show_error(error);
        moist_delete_cavity(&cav);
        moist_delete_structure(&mol);
        moist_delete_error(&error);
        return 1;
    }

    printf("  Cavity grid points: %d\n", ngrid);
    printf("  Number of spheres: %d\n", nsph);

    // Allocate arrays dynamically (C++ compatible, no VLAs)
    double* xyz = (double*)malloc(3 * ngrid * sizeof(double));
    double* a = (double*)malloc(ngrid * sizeof(double));
    int* owner = (int*)malloc(ngrid * sizeof(int));
    bool* converged = (bool*)malloc(ngrid * sizeof(bool));
    double* vradii = (double*)malloc(nsph * sizeof(double));
    double* asph = (double*)malloc(nsph * sizeof(double));

    if (!xyz || !a || !owner || !converged || !vradii || !asph) {
        printf("Error: Memory allocation failed\n");
        free(xyz); free(a); free(owner); free(converged); free(vradii); free(asph);
        moist_delete_cavity(&cav);
        moist_delete_structure(&mol);
        moist_delete_error(&error);
        return 1;
    }

    // Get generic cavity results
    double area = 0.0, volume = 0.0;
    moist_get_cavity_results(error, cav, &area, &volume, &ngrid, &nsph,
                             xyz, a, owner, converged, vradii, asph);
    if (moist_check_error(error)) {
        show_error(error);
        free(xyz); free(a); free(owner); free(converged); free(vradii); free(asph);
        moist_delete_cavity(&cav);
        moist_delete_structure(&mol);
        moist_delete_error(&error);
        return 1;
    }

    printf("  Cavity surface area: %.4f Bohr²\n", area);
    printf("  Cavity volume: %.4f Bohr³\n", volume);

    // Allocate DROP-specific arrays dynamically
    double* normal = (double*)malloc(3 * ngrid * sizeof(double));
    double* wleb = (double*)malloc(ngrid * sizeof(double));
    double* r_iI0 = (double*)malloc(ngrid * sizeof(double));
    double* f = (double*)malloc(ngrid * sizeof(double));
    double* rho = (double*)malloc(ngrid * sizeof(double));
    int nmax = 0;

    if (!normal || !wleb || !r_iI0 || !f || !rho) {
        printf("Error: Memory allocation failed for DROP arrays\n");
        free(xyz); free(a); free(owner); free(converged); free(vradii); free(asph);
        free(normal); free(wleb); free(r_iI0); free(f); free(rho);
        moist_delete_cavity(&cav);
        moist_delete_structure(&mol);
        moist_delete_error(&error);
        return 1;
    }

    moist_get_drop_specific(error, cav, &nmax, normal, wleb, r_iI0,
                           f, rho);
    if (moist_check_error(error)) {
        show_error(error);
        free(xyz); free(a); free(owner); free(converged); free(vradii); free(asph);
        free(normal); free(wleb); free(r_iI0); free(f); free(rho);
        moist_delete_cavity(&cav);
        moist_delete_structure(&mol);
        moist_delete_error(&error);
        return 1;
    }

    printf("  Raw grid size (nmax): %d\n", nmax);

    // Validate results
    int result = area > 0.0 && volume > 0.0 && ngrid > 0 && nsph == 3 &&
                 nmax >= ngrid && wleb[0] > 0.0 && a[0] > 0.0 &&
                 vradii[0] > 0.0 && vradii[1] > 0.0 && vradii[2] > 0.0
             ? 0
             : 1;

    // Free allocated memory
    free(xyz); free(a); free(owner); free(converged); free(vradii); free(asph);
    free(normal); free(wleb); free(r_iI0); free(f); free(rho);

    // Cleanup handles
    moist_delete_cavity(&cav);
    moist_delete_structure(&mol);
    moist_delete_error(&error);

    if (cav != NULL || mol != NULL || error != NULL) {
        printf("Error: handles not reset to NULL after deletion\n");
        return 1;
    }

    printf("  H2O cavity computation completed successfully!\n");
    return result;
}

int test_cavity_gradient(void)
{
    printf("Start test: cavity gradient computation\n");
    moist_error error = moist_new_error();

    // Define H2O molecule (coordinates in Bohr)
    const int natoms = 3;
    int numbers[3] = {8, 1, 1};  // O, H, H
    double positions[9] = {
        0.0000,  0.0000,  0.1173,   // O atom
        0.0000,  1.4309, -0.9370,   // H atom 1
        0.0000, -1.4309, -0.9370    // H atom 2
    };

    moist_structure mol = moist_new_structure(error, natoms, numbers, positions, NULL, NULL);
    if (moist_check_error(error)) {
        show_error(error);
        moist_delete_error(&error);
        return 1;
    }

    // Create DROP cavity
    moist_cavity cav = moist_new_drop_cavity(error, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL);
    if (moist_check_error(error)) {
        show_error(error);
        moist_delete_structure(&mol);
        moist_delete_error(&error);
        return 1;
    }

    // Build cavity
    moist_update_cavity(error, cav, mol);
    if (moist_check_error(error)) {
        show_error(error);
        moist_delete_cavity(&cav);
        moist_delete_structure(&mol);
        moist_delete_error(&error);
        return 1;
    }

    // Get cavity sizes
    int ngrid = 0, nsph = 0;
    moist_get_cavity_sizes(error, cav, &ngrid, &nsph);
    if (moist_check_error(error)) {
        show_error(error);
        moist_delete_cavity(&cav);
        moist_delete_structure(&mol);
        moist_delete_error(&error);
        return 1;
    }

    printf("  Grid: ngrid=%d, nsph=%d\n", ngrid, nsph);

    // Compute cavity gradient w.r.t. nuclear coordinates
    printf("  Computing cavity gradient...\n");
    moist_compute_cavity_gradient(error, cav);
    if (moist_check_error(error)) {
        show_error(error);
        moist_delete_cavity(&cav);
        moist_delete_structure(&mol);
        moist_delete_error(&error);
        return 1;
    }

    // Allocate gradient arrays
    double* A_tot1_rA = (double*)malloc(3 * nsph * sizeof(double));              // (3, nsph)
    double* V_tot1_rA = (double*)malloc(3 * nsph * sizeof(double));              // (3, nsph)
    double* asph1_rA = (double*)malloc(3 * nsph * nsph * sizeof(double));        // (3, nsph, nsph)
    double* vsph1_rA = (double*)malloc(3 * nsph * nsph * sizeof(double));        // (3, nsph, nsph)
    double* xyz1_rA = (double*)malloc(3 * 3 * nsph * ngrid * sizeof(double));  // (3, 3, nsph, ngrid)
    double* r_iI1_rA = (double*)malloc(3 * nsph * ngrid * sizeof(double));       // (3, nsph, ngrid)
    double* rho1_rA = (double*)malloc(3 * nsph * ngrid * sizeof(double));        // (3, nsph, ngrid)

    if (!A_tot1_rA || !V_tot1_rA || !asph1_rA || !vsph1_rA || 
        !xyz1_rA || !r_iI1_rA || !rho1_rA) {
        printf("Error: Memory allocation failed for gradient arrays\n");
        free(A_tot1_rA); free(V_tot1_rA); free(asph1_rA); free(vsph1_rA);
        free(xyz1_rA); free(r_iI1_rA); free(rho1_rA);
        moist_delete_cavity(&cav);
        moist_delete_structure(&mol);
        moist_delete_error(&error);
        return 1;
    }

    // Get cavity gradient arrays
    moist_get_cavity_gradient(error, cav, &nsph, &ngrid,
                              A_tot1_rA, V_tot1_rA, 
                              asph1_rA, vsph1_rA,
                              xyz1_rA, r_iI1_rA, rho1_rA);
    if (moist_check_error(error)) {
        show_error(error);
        free(A_tot1_rA); free(V_tot1_rA); free(asph1_rA); free(vsph1_rA);
        free(xyz1_rA); free(r_iI1_rA); free(rho1_rA);
        moist_delete_cavity(&cav);
        moist_delete_structure(&mol);
        moist_delete_error(&error);
        return 1;
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
    double* Amat0 = (double*)malloc(ngrid * ngrid * sizeof(double));
    double* Amat1_rA = (double*)malloc(3 * nsph * ngrid * ngrid * sizeof(double));
    double* xi = (double*)malloc(ngrid * sizeof(double));

    if (!Amat0 || !Amat1_rA || !xi) {
        printf("Error: Memory allocation failed for A-matrix arrays\n");
        free(A_tot1_rA); free(V_tot1_rA); free(asph1_rA); free(vsph1_rA);
        free(xyz1_rA); free(r_iI1_rA); free(rho1_rA);
        free(Amat0); free(Amat1_rA); free(xi);
        moist_delete_cavity(&cav);
        moist_delete_structure(&mol);
        moist_delete_error(&error);
        return 1;
    }

    moist_get_amat_gradient(error, cav, &nsph, &ngrid, Amat0, Amat1_rA, xi);
    if (moist_check_error(error)) {
        show_error(error);
        free(A_tot1_rA); free(V_tot1_rA); free(asph1_rA); free(vsph1_rA);
        free(xyz1_rA); free(r_iI1_rA); free(rho1_rA);
        free(Amat0); free(Amat1_rA); free(xi);
        moist_delete_cavity(&cav);
        moist_delete_structure(&mol);
        moist_delete_error(&error);
        return 1;
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
    double* q1 = (double*)malloc(ngrid * sizeof(double));
    double* q2 = (double*)malloc(ngrid * sizeof(double));
    double* grad_contract = (double*)malloc(3 * nsph * sizeof(double));    // (3, nsph)
    double* grad_ref = (double*)malloc(3 * nsph * sizeof(double));         // (3, nsph)

    // Test API: contract_nuc_elec_qefield_rA
    double* surface_q = (double*)malloc(ngrid * sizeof(double));           // (ngrid)
    double* surface_q_scaled = (double*)malloc(ngrid * sizeof(double));    // (ngrid)
    double* qefield = (double*)malloc(3 * ngrid * sizeof(double));         // (3, ngrid)
    double* qefield_scaled = (double*)malloc(3 * ngrid * sizeof(double));  // (3, ngrid)
    double* za = (double*)malloc(nsph * sizeof(double));                    // (nsph)
    double* grad_ne_zero = (double*)malloc(3 * nsph * sizeof(double));      // (3, nsph)
    double* grad_ne = (double*)malloc(3 * nsph * sizeof(double));           // (3, nsph)
    double* grad_ne_scaled = (double*)malloc(3 * nsph * sizeof(double));    // (3, nsph)

    if (!q1 || !q2 || !grad_contract || !grad_ref ||
        !surface_q || !surface_q_scaled || !qefield || !qefield_scaled ||
        !za || !grad_ne_zero || !grad_ne || !grad_ne_scaled) {
        printf("Error: Memory allocation failed for contraction API test arrays\n");
        free(A_tot1_rA); free(V_tot1_rA); free(asph1_rA); free(vsph1_rA);
        free(xyz1_rA); free(r_iI1_rA); free(rho1_rA);
        free(Amat0); free(Amat1_rA); free(xi);
        free(q1); free(q2); free(grad_contract); free(grad_ref);
        free(surface_q); free(surface_q_scaled); free(qefield); free(qefield_scaled);
        free(za); free(grad_ne_zero); free(grad_ne); free(grad_ne_scaled);
        moist_delete_cavity(&cav);
        moist_delete_structure(&mol);
        moist_delete_error(&error);
        return 1;
    }

    for (int i = 0; i < ngrid; i++) {
        q1[i] = xi[i];
        q2[i] = 0.5 * xi[i] + 1.0e-3 * (double)(i + 1);
    }

    moist_contract_amat1_q1q2_rA(error, cav, q1, q2, grad_contract);
    if (moist_check_error(error)) {
        show_error(error);
        free(A_tot1_rA); free(V_tot1_rA); free(asph1_rA); free(vsph1_rA);
        free(xyz1_rA); free(r_iI1_rA); free(rho1_rA);
        free(Amat0); free(Amat1_rA); free(xi);
        free(q1); free(q2); free(grad_contract); free(grad_ref);
        free(surface_q); free(surface_q_scaled); free(qefield); free(qefield_scaled);
        free(za); free(grad_ne_zero); free(grad_ne); free(grad_ne_scaled);
        moist_delete_cavity(&cav);
        moist_delete_structure(&mol);
        moist_delete_error(&error);
        return 1;
    }

    // Symmetry check: for symmetric dA, q1^T(dA)q2 == q2^T(dA)q1
    moist_contract_amat1_q1q2_rA(error, cav, q2, q1, grad_ref);
    if (moist_check_error(error)) {
        show_error(error);
        free(A_tot1_rA); free(V_tot1_rA); free(asph1_rA); free(vsph1_rA);
        free(xyz1_rA); free(r_iI1_rA); free(rho1_rA);
        free(Amat0); free(Amat1_rA); free(xi);
        free(q1); free(q2); free(grad_contract); free(grad_ref);
        free(surface_q); free(surface_q_scaled); free(qefield); free(qefield_scaled);
        free(za); free(grad_ne_zero); free(grad_ne); free(grad_ne_scaled);
        moist_delete_cavity(&cav);
        moist_delete_structure(&mol);
        moist_delete_error(&error);
        return 1;
    }

    int contract_ok = 1;
    const double contract_tol = 1e-9;
    double contract_norm = 0.0;
    for (int iatom = 0; iatom < nsph && contract_ok; iatom++) {
        for (int iaxis = 0; iaxis < 3; iaxis++) {
            size_t idx = idx_f2(iaxis, iatom, 3);
            if (fabs(grad_contract[idx] - grad_ref[idx]) > contract_tol) {
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
    if (contract_ok && contract_norm < 1e-20) {
        printf("  Warning: contract_amat1_q1q2_rA norm is essentially zero\n");
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
        free(A_tot1_rA); free(V_tot1_rA); free(asph1_rA); free(vsph1_rA);
        free(xyz1_rA); free(r_iI1_rA); free(rho1_rA);
        free(Amat0); free(Amat1_rA); free(xi);
        free(q1); free(q2); free(grad_contract); free(grad_ref);
        free(surface_q); free(surface_q_scaled); free(qefield); free(qefield_scaled);
        free(za); free(grad_ne_zero); free(grad_ne); free(grad_ne_scaled);
        moist_delete_cavity(&cav);
        moist_delete_structure(&mol);
        moist_delete_error(&error);
        return 1;
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
        free(A_tot1_rA); free(V_tot1_rA); free(asph1_rA); free(vsph1_rA);
        free(xyz1_rA); free(r_iI1_rA); free(rho1_rA);
        free(Amat0); free(Amat1_rA); free(xi);
        free(q1); free(q2); free(grad_contract); free(grad_ref);
        free(surface_q); free(surface_q_scaled); free(qefield); free(qefield_scaled);
        free(za); free(grad_ne_zero); free(grad_ne); free(grad_ne_scaled);
        moist_delete_cavity(&cav);
        moist_delete_structure(&mol);
        moist_delete_error(&error);
        return 1;
    }

    moist_contract_nuc_elec_qefield_rA(error, cav, surface_q_scaled, qefield_scaled, za, grad_ne_scaled);
    if (moist_check_error(error)) {
        show_error(error);
        free(A_tot1_rA); free(V_tot1_rA); free(asph1_rA); free(vsph1_rA);
        free(xyz1_rA); free(r_iI1_rA); free(rho1_rA);
        free(Amat0); free(Amat1_rA); free(xi);
        free(q1); free(q2); free(grad_contract); free(grad_ref);
        free(surface_q); free(surface_q_scaled); free(qefield); free(qefield_scaled);
        free(za); free(grad_ne_zero); free(grad_ne); free(grad_ne_scaled);
        moist_delete_cavity(&cav);
        moist_delete_structure(&mol);
        moist_delete_error(&error);
        return 1;
    }

    const double homo_tol = 1e-9;
    for (int iatom = 0; iatom < nsph && nuc_elec_ok; iatom++) {
        for (int iaxis = 0; iaxis < 3; iaxis++) {
            size_t idx = idx_f2(iaxis, iatom, 3);
            if (!isfinite(grad_ne[idx]) || !isfinite(grad_ne_scaled[idx])) {
                printf("  Error: contract_nuc_elec_qefield_rA produced non-finite values\n");
                nuc_elec_ok = 0;
                break;
            }
            if (fabs(grad_ne_scaled[idx] - 2.0 * grad_ne[idx]) > homo_tol) {
                printf("  Error: contract_nuc_elec_qefield_rA homogeneity check failed\n");
                nuc_elec_ok = 0;
                break;
            }
        }
    }

    // Validate results - check that gradients are computed and non-trivial
    int result = 0;
    // A_tot1_rA should have non-zero values (area depends on nuclear coords)
    double grad_sum = 0.0;
    for (int i = 0; i < 3 * nsph; i++) {
        grad_sum += A_tot1_rA[i] * A_tot1_rA[i];
    }
    if (grad_sum < 1e-20) {
        printf("  Warning: Area gradient is essentially zero\n");
        result = 1;
    }
    if (!layout_ok) {
        result = 1;
    }
    if (!contract_ok || !nuc_elec_ok) {
        result = 1;
    }

    // Free all arrays
    free(A_tot1_rA); free(V_tot1_rA); free(asph1_rA); free(vsph1_rA);
    free(xyz1_rA); free(r_iI1_rA); free(rho1_rA);
    free(Amat0); free(Amat1_rA); free(xi);
    free(q1); free(q2); free(grad_contract); free(grad_ref);
    free(surface_q); free(surface_q_scaled); free(qefield); free(qefield_scaled);
    free(za); free(grad_ne_zero); free(grad_ne); free(grad_ne_scaled);

    // Cleanup handles
    moist_delete_cavity(&cav);
    moist_delete_structure(&mol);
    moist_delete_error(&error);

    printf("  Gradient computation completed successfully!\n");
    return result;
}

int main(void)
{
    int stat = 0;
    stat += test_header_and_version();
    stat += test_version();
    stat += test_uninitialized_error();
    stat += test_null_handle();
    stat += test_delete_resets_handle();
    // stat += test_gems_model();
    // stat += test_gems_unknown_solvent_reports_error();
    stat += test_drop_cavity();
    stat += test_custom_radii();
    stat += test_h2o_cavity();
    stat += test_cavity_gradient();

    return stat == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}
