/* C shim exposing the vendored libmsym point-group detection to Fortran
 *
 * libmsym hands out pointers into its own context (operations, subgroups),
 * enums and fixed-size name buffers; this file flattens that into plain arrays
 * so moist_math_symmetry_msym (msym.f90) binds these entry points and nothing else
 *
 * - every element gets unit mass, the caller's species id as nuclear charge
 *   and a fixed name: libmsym calls two atoms equivalent when all three match,
 *   so the species id alone decides which atoms may be permuted, and the
 *   center libmsym works about is the plain centroid (no isotope dependence)
 * - operation matrices are written column-major (Fortran order), each acting
 *   on coordinates relative to `center`
 * - subgroup operations are returned as indices into the operation list,
 *   found by matrix comparison, not by pointer identity into libmsym's arrays
 * - detection runs at libmsym's default thresholds except `equivalence`,
 *   which the caller passes: the default merges unrelated orbits of a large
 *   molecule into one set, which libmsym then rejects or reads as a lower group
 * - every entry point returns 0 on success and a nonzero status otherwise;
 *   moist_msym_error_message() renders the status and libmsym's last details
 *
 * - libmsym's floating-point exceptions stay inside: it divides zero by zero
 *   on purpose (vequal), which a caller trapping exceptions (gfortran
 *   -ffpe-trap) would die on; every entry point computing in libmsym holds
 *   the caller's environment, then drops the raised flags and restores it
 *
 * libmsym keeps its error details in a static buffer, so calls must not run
 * concurrently; moist calls this once per cavity update, outside parallel
 * regions
 */
#include <fenv.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "msym.h"
#include "msym_error.h"
#include "symop.h"

/* Largest elementwise difference for two operation matrices to match */
#define MOIST_MSYM_MATCH 1.0e-6

/* Size of libmsym's error details buffer (private to its msym_error.c) */
#define MOIST_MSYM_DETAILS 1024

/* Status codes beyond libmsym's own (all negative there) */
#define MOIST_MSYM_ALLOC 1
#define MOIST_MSYM_BOUNDS 2
#define MOIST_MSYM_NOMATCH 3

static void copy_name(char *dst, const char *src) {
  /* Fixed 8-byte, NUL-padded name for Fortran */
  memset(dst, 0, 8);
  strncpy(dst, src, 7);
}

/* Caller's floating-point environment while libmsym computes */
typedef struct {
  fenv_t env;
  int held;
} fp_guard;

static void fp_hold(fp_guard *guard) {
  /* Saves the environment, clears the flags and masks the traps; nonzero
   * where non-stop mode does not exist, with nothing to restore */
  guard->held = feholdexcept(&guard->env) == 0;
}

static void fp_restore(const fp_guard *guard) {
  if (!guard->held) return;
  /* Flags first: restoring the traps over a raised flag must not fire one */
  feclearexcept(FE_ALL_EXCEPT);
  fesetenv(&guard->env);
}

static void column_major(double m[3][3], double *out) {
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) out[i + 3 * j] = m[i][j];
}

int moist_msym_create(void **handle) {
  msym_context ctx = msymCreateContext();
  *handle = (void *)ctx;
  return ctx == NULL ? MOIST_MSYM_ALLOC : 0;
}

void moist_msym_release(void *handle) {
  if (handle != NULL) msymReleaseContext((msym_context)handle);
}

static int find_symmetry(void *handle, int nat, const double *xyz,
                         const int *species, double equivalence) {
  msym_context ctx = (msym_context)handle;
  /* Default thresholds with the caller's equivalence test; not positive
   * keeps libmsym's own */
  if (equivalence > 0.0) {
    msym_thresholds_t thresholds = *msymGetDefaultThresholds();
    thresholds.equivalence = equivalence;
    msym_error_t tret = msymSetThresholds(ctx, &thresholds);
    if (tret != MSYM_SUCCESS) return (int)tret;
  }
  msym_element_t *elements =
      calloc(nat > 0 ? (size_t)nat : 1, sizeof(msym_element_t));
  if (elements == NULL) return MOIST_MSYM_ALLOC;
  for (int i = 0; i < nat; i++) {
    elements[i].id = NULL;
    elements[i].m = 1.0;
    elements[i].n = species[i];
    elements[i].name[0] = 'X';
    for (int k = 0; k < 3; k++) elements[i].v[k] = xyz[3 * i + k];
  }
  /* libmsym copies the elements */
  msym_error_t ret = msymSetElements(ctx, nat, elements);
  free(elements);
  if (ret != MSYM_SUCCESS) return (int)ret;
  ret = msymFindSymmetry(ctx);
  /* libmsym cannot orient a group without any axis or plane and reports C1
   * as this failure; set C1 itself so the caller sees the identity */
  if (ret == MSYM_POINT_GROUP_ERROR) {
    /* Reading the details clears libmsym's buffer: keep a copy and put it
     * back for moist_msym_error_message() when this is another failure */
    char details[MOIST_MSYM_DETAILS];
    snprintf(details, sizeof(details), "%s", msymGetErrorDetails());
    if (strcmp(details, "Point group has no primary axis for reorientation") == 0)
      ret = msymSetPointGroupByName(ctx, "C1");
    else
      msymSetErrorDetails("%s", details);
  }
  return (int)ret;
}

int moist_msym_find(void *handle, int nat, const double *xyz,
                    const int *species, double equivalence) {
  fp_guard guard;
  fp_hold(&guard);
  int status = find_symmetry(handle, nat, xyz, species, equivalence);
  fp_restore(&guard);
  return status;
}

static int group_summary(void *handle, char *name, int *type, int *n,
                         int *order) {
  msym_context ctx = (msym_context)handle;
  char buf[8] = {0};
  msym_point_group_type_t t;
  const msym_symmetry_operation_t *sops = NULL;
  msym_error_t ret;
  if ((ret = msymGetPointGroupName(ctx, sizeof(buf), buf)) != MSYM_SUCCESS)
    return (int)ret;
  if ((ret = msymGetPointGroupType(ctx, &t, n)) != MSYM_SUCCESS) return (int)ret;
  if ((ret = msymGetSymmetryOperations(ctx, order, &sops)) != MSYM_SUCCESS)
    return (int)ret;
  copy_name(name, buf);
  *type = (int)t;
  return 0;
}

int moist_msym_group(void *handle, char *name, int *type, int *n, int *order) {
  fp_guard guard;
  fp_hold(&guard);
  int status = group_summary(handle, name, type, n, order);
  fp_restore(&guard);
  return status;
}

static int operations(void *handle, int order, double *mats, int *info,
                      double *center) {
  msym_context ctx = (msym_context)handle;
  const msym_symmetry_operation_t *sops = NULL;
  int sopsl = 0;
  msym_error_t ret;
  if ((ret = msymGetSymmetryOperations(ctx, &sopsl, &sops)) != MSYM_SUCCESS)
    return (int)ret;
  if (sopsl != order) return MOIST_MSYM_BOUNDS;
  if ((ret = msymGetCenterOfMass(ctx, center)) != MSYM_SUCCESS) return (int)ret;
  for (int k = 0; k < sopsl; k++) {
    double m[3][3];
    msym_symmetry_operation_t sop = sops[k];
    symmetryOperationMatrix(&sop, m);
    column_major(m, mats + 9 * k);
    info[4 * k + 0] = (int)sop.type;
    info[4 * k + 1] = sop.order;
    info[4 * k + 2] = sop.power;
    info[4 * k + 3] = (int)sop.orientation;
  }
  return 0;
}

int moist_msym_operations(void *handle, int order, double *mats, int *info,
                          double *center) {
  fp_guard guard;
  fp_hold(&guard);
  int status = operations(handle, order, mats, info, center);
  fp_restore(&guard);
  return status;
}

static int operation_axis(void *handle, int iop, double *axis) {
  /* Rotation axis / plane normal of operation iop (0-based) */
  msym_context ctx = (msym_context)handle;
  const msym_symmetry_operation_t *sops = NULL;
  int sopsl = 0;
  msym_error_t ret;
  if ((ret = msymGetSymmetryOperations(ctx, &sopsl, &sops)) != MSYM_SUCCESS)
    return (int)ret;
  if (iop < 0 || iop >= sopsl) return MOIST_MSYM_BOUNDS;
  for (int k = 0; k < 3; k++) axis[k] = sops[iop].v[k];
  return 0;
}

int moist_msym_axis(void *handle, int iop, double *axis) {
  fp_guard guard;
  fp_hold(&guard);
  int status = operation_axis(handle, iop, axis);
  fp_restore(&guard);
  return status;
}

static int subgroup_count(void *handle, int *nsub) {
  msym_context ctx = (msym_context)handle;
  const msym_subgroup_t *sg = NULL;
  msym_error_t ret = msymGetSubgroups(ctx, nsub, &sg);
  if (ret != MSYM_SUCCESS) *nsub = 0;
  return (int)ret;
}

int moist_msym_subgroup_count(void *handle, int *nsub) {
  fp_guard guard;
  fp_hold(&guard);
  int status = subgroup_count(handle, nsub);
  fp_restore(&guard);
  return status;
}

static int subgroup(void *handle, int isub, char *name, int *order,
                    int maxorder, int *ops) {
  msym_context ctx = (msym_context)handle;
  const msym_subgroup_t *sg = NULL;
  const msym_symmetry_operation_t *sops = NULL;
  int sgl = 0, sopsl = 0;
  msym_error_t ret;
  if ((ret = msymGetSubgroups(ctx, &sgl, &sg)) != MSYM_SUCCESS) return (int)ret;
  if ((ret = msymGetSymmetryOperations(ctx, &sopsl, &sops)) != MSYM_SUCCESS)
    return (int)ret;
  if (isub < 0 || isub >= sgl) return MOIST_MSYM_BOUNDS;
  const msym_subgroup_t *s = &sg[isub];
  if (s->order > maxorder) return MOIST_MSYM_BOUNDS;
  copy_name(name, s->name);
  *order = s->order;
  for (int i = 0; i < s->order; i++) {
    double a[3][3];
    msym_symmetry_operation_t sop = *s->sops[i];
    symmetryOperationMatrix(&sop, a);
    int found = -1;
    for (int k = 0; k < sopsl && found < 0; k++) {
      double b[3][3], diff = 0.0;
      msym_symmetry_operation_t gop = sops[k];
      symmetryOperationMatrix(&gop, b);
      for (int p = 0; p < 3; p++)
        for (int q = 0; q < 3; q++) diff = fmax(diff, fabs(a[p][q] - b[p][q]));
      if (diff < MOIST_MSYM_MATCH) found = k;
    }
    if (found < 0) return MOIST_MSYM_NOMATCH;
    ops[i] = found;
  }
  return 0;
}

int moist_msym_subgroup(void *handle, int isub, char *name, int *order,
                        int maxorder, int *ops) {
  fp_guard guard;
  fp_hold(&guard);
  int status = subgroup(handle, isub, name, order, maxorder, ops);
  fp_restore(&guard);
  return status;
}

void moist_msym_error_message(int status, char *buf, int len) {
  if (len <= 0) return;
  const char *what;
  const char *details = "";
  switch (status) {
    case 0: what = "success"; break;
    case MOIST_MSYM_ALLOC: what = "allocation failed"; break;
    case MOIST_MSYM_BOUNDS: what = "index or size out of bounds"; break;
    case MOIST_MSYM_NOMATCH:
      what = "subgroup operation not found in the point group";
      break;
    default:
      what = msymErrorString((msym_error_t)status);
      details = msymGetErrorDetails();
      break;
  }
  snprintf(buf, (size_t)len, "%s%s%s", what, details[0] ? ": " : "", details);
}
