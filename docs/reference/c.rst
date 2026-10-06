C API
=====

``include/moist.h`` supports C11 and C++11. Names follow ``moist_<verb>_<object>[_<detail>]``.

Construction
------------

Construct a run context, radii and a level-set function (LSF), then a cavity
and model. The model requires the context; a cavity or component without one
runs on its model's.
DROP requires an LSF; iSwiG does not. NULL radii selects CPCM radii.
NULL options selects defaults. To override settings, initialize the typed
options with the caller's ``sizeof`` before assigning fields:

.. code-block:: c

   moist_error error = moist_new_error();
   moist_context context = moist_new_context(error, 4, 0, false);  /* nthreads, verbosity, debug */
   moist_drop_options options;
   moist_init_drop_options(error, &options, sizeof options);
   /* Check moist_check_error(error) after each fallible call. */
   options.nleb = 194;

   moist_lsf lsf = moist_new_svdw_lsf(error, NULL);
   moist_cavity cavity = moist_new_drop_cavity(error, NULL, lsf, NULL, &options);  /* NULL: the model's context */
   moist_delete(lsf);  /* cavity owns a copy */
   moist_model model = moist_new_model(error, context, cavity);
   moist_delete(cavity);  /* model owns a copy */
   moist_delete(context);  /* model retains the shared context */

   moist_component pcm = moist_new_cpcm_component(error, NULL, 78.4, NULL);
   moist_add_model_component(error, model, pcm);  /* model owns a copy */
   moist_delete(pcm);
   moist_structure mol = moist_new_structure(error, natoms, numbers, positions,
                                             NULL, NULL);  /* NULL: molecular */
   moist_update_model(error, model, mol);
   moist_delete(mol);  /* the cavity and components copy the structure */

Add every component before the first ``moist_update_model``; later additions
are rejected.

SvdW, CFC, isodensity, DROP, iSwiG and PCM settings have separate
options types. ``struct_size`` is set by the initializer; do not modify it.
The complete 1.0 layout is the minimum size; existing fields keep their offsets
and meanings, and future libraries default fields an older struct lacks.
Larger structs are accepted;
unknown trailing fields are ignored by constructors and untouched by initializers.

``moist_new_context(error, nthreads, verbosity, debug)`` mirrors the Fortran ``new_context(ctx, nthreads=, verbosity=, debug=)``; C has no optional arguments, so pass ``0, 0, false`` for the defaults.
The thread count is fixed for the context's lifetime: a positive count is used as given, ``0`` takes the calling thread's ``omp_get_max_threads()`` at construction, and a negative count is an error.
``verbosity`` and ``debug`` are the only logging settings; the options structs carry none.

Every cavity, component and model constructor takes the context as its second argument.
``moist_new_model`` refuses a NULL context; a cavity or component accepts NULL.
A cavity or component copied into a model keeps its own context, or runs on the model's when it has none.
Without a context, updating a standalone cavity is an error.
A cavity borrowed with ``moist_get_model_cavity`` retains the context its copy runs on until it is deleted.

Every handle retains its context, and a model retains those of its parts, so the public context handle can be deleted.
``moist_get_context_num_threads`` reports the count the context was constructed with.
The count sizes MOIST's own OpenMP regions and FFT workers; MOIST never changes the host's OpenMP runtime. BLAS/LAPACK threading is the host's to configure (``OMP_NUM_THREADS``, ``MKL_NUM_THREADS``, ``OPENBLAS_NUM_THREADS`` or the vendor's API).
Without OpenMP, MOIST's kernels are serial; the math backends may still use threads.

Ownership
---------

Each ``new_*`` handle has one owner. Constructors copy their configuration;
assigning a handle creates an alias, not another owner. Delete couplings before
their model. ``moist_get_model_cavity`` returns a borrowed cavity: it must not
outlive the model, ``moist_update_cavity`` rejects it, and
``moist_delete_cavity`` releases only the handle.
Responses own their results. Callback functions and contexts must outlive every
cavity or model using them.

``moist_delete(handle)`` dispatches to the typed deletion function in C and C++.
Deletion sets the supplied handle to NULL; deleting an already-NULL handle is safe.
Requests, outputs, response items and fields are addressed by name
(NUL-terminated strings of at most ``MOIST_NAME_MAX`` or ``MOIST_FIELD_NAME_MAX``
characters); there are no numeric tags. Requests and response items have no
handles or indices: a cursor makes one current at a time.

Serialize calls sharing a context, handles, their borrowed aliases, or their parent model.
Independent object graphs may run on separate threads. Isodensity callbacks may
run concurrently with the same context; they must be reentrant and protect shared
scratch data.

Errors and buffers
------------------

Fallible calls take a required ``moist_error`` first and replace its previous diagnostic.
Check ``moist_check_error`` before using results: ``moist_success`` means success,
``moist_failure`` means a diagnostic, and ``moist_invalid_error`` means a NULL
error handle. A NULL required pointer is an error; optional pointers have
documented meanings. Numeric outputs remain unchanged on failure.
Banner, version, and field-description getters clear writable buffers on failure.
Rejected coupling answers follow the invalidation rules in :doc:`coupling`.

``moist_get_error`` takes a buffer and a pointer to its positive
``int`` capacity. It NUL-terminates and truncates; no error produces an empty string.

``moist_get_banner``, ``moist_get_version_string`` and the
``_field_about`` getters take capacity by value and return
the full text length through a required ``size_t *length``, excluding the terminator.
Pass NULL and zero to query the length, then allocate ``length + 1`` bytes.
Truncation is successful and detectable as ``length >= capacity``. Positive-capacity buffers are always
NUL-terminated; errors leave ``length`` unchanged. The host prints banner text
to its own stream. Capacities above ``SIZE_MAX/2`` are rejected.

Arrays use flat C row-major order, with the last axis contiguous.
Dimensions reverse the native Fortran dimensions without rearranging the buffer.
For example, positions and gradients are ``[natoms][3]``; grid vectors are
``[ngrid][3]``. The entries
of the host loop -- the coupling and response entries below and the cavity
and model field getters -- take no size: moist reads or writes exactly the documented
shape, and the host allocates it. ``moist_get_cavity_results``,
``moist_get_cavity_gaussian``, ``moist_assemble_amat``,
``moist_get_model_gradient`` and the gradient-tensor getters take ``int``
capacities, which may exceed the logical size; padding remains untouched.
The ``moist_contract_*`` entries use the cavity's own ``ngrid`` and ``nsph``.
The header defines shapes and strides.
See :doc:`coupling` for the host exchange protocol.

.. _coupling-c:

Host coupling
-------------

Staging (``moist_prepare_model_energy``, ``moist_prepare_model_response``,
``moist_prepare_model_gradient``) and evaluation (``moist_get_model_energy``,
``moist_get_model_response``, ``moist_get_model_gradient``) take the model and
the coupling; :doc:`coupling` describes the shared protocol. Requests and response
items are reached through these entries (declaration macros omitted):

.. code-block:: c

   bool moist_next_coupling_request(moist_error error, moist_coupling coupling);
   void moist_get_coupling_request_name(moist_error error, moist_coupling coupling,
                                        char* name);
   void moist_get_coupling_request_missing(moist_error error, moist_coupling coupling,
                                           const char* output, bool* missing);
   void moist_answer_coupling_request(moist_error error, moist_coupling coupling,
                                      const char* output, const double* values);
   void moist_get_coupling_request_field_count(moist_error error, moist_coupling coupling,
                                               int* nfield);
   void moist_get_coupling_request_field_info(moist_error error, moist_coupling coupling,
                                              int index,
                                              char* name /* [MOIST_FIELD_NAME_MAX + 1] */,
                                              int* dtype, int* rank,
                                              int* dims /* [MOIST_FIELD_MAX_RANK] */,
                                              int* count);
   void moist_get_coupling_request_field_about(moist_error error, moist_coupling coupling,
                                               const char* name, char* about,
                                               size_t capacity, size_t* length);
   void moist_get_coupling_request_field_real(moist_error error, moist_coupling coupling,
                                              const char* name, double* values);
   bool moist_next_response_item(moist_error error, moist_response response);
   void moist_get_response_item_name(moist_error error, moist_response response,
                                     char* name);
   void moist_get_response_field_count(moist_error error, moist_response response,
                                       int* nfield);
   void moist_get_response_field_info(moist_error error, moist_response response,
                                      int index,
                                      char* name /* [MOIST_FIELD_NAME_MAX + 1] */,
                                      int* dtype, int* rank,
                                      int* dims /* [MOIST_FIELD_MAX_RANK] */,
                                      int* count);
   void moist_get_response_field_about(moist_error error, moist_response response,
                                       const char* name, char* about,
                                       size_t capacity, size_t* length);
   void moist_get_response_field_real(moist_error error, moist_response response,
                                      const char* name, double* values);

The two cursor entries, ``moist_next_coupling_request`` and
``moist_next_response_item``, are the only ones here that return a value. They
return false at the end of a pass and on any failure, so an error cannot keep a
host loop running; check the error after the loop. A NULL error handle gives
false without a diagnostic; a NULL or invalid coupling and a NULL response give
false with one. The name, missing, answer and request field entries act on the
current request, the item name and response field entries on the current
response item; each reports an error when none is current.
``moist_get_coupling_request_missing`` writes false for an output the request
does not have; ``moist_answer_coupling_request`` rejects that name.

The field entries read request inputs and response arrays as :doc:`fields`,
with the conventions of the model field getters; only real fields exist.
``gaussian_moments`` declares ``width[ngrid]``; the other requests declare no
field. ``_count`` and ``_info`` enumerate whatever the current item declares,
so ``count`` sizes a buffer for any item. The items of the
:ref:`response table <coupling-response>` declare ``w_phi[ngrid]`` of
``potential_adjoint``; ``w_rho[ngrid]``, ``w_grad_rho[ngrid][3]`` and
``w_hess_rho[ngrid][3][3]`` of ``density``; ``w_overlap[ngrid]`` and
``w_normal_deriv[ngrid]`` of ``gaussian_amplitude``. ``w_hess_rho`` is
``[point][b][a]`` for native ``(a,b,point)`` and need not be symmetric;
``_real`` writes exactly the ``count`` elements ``_info`` reports. No current
item or request, a field it does not declare and a NULL buffer are errors
reported by name.

Evaluation example
~~~~~~~~~~~~~~~~~~

This example evaluates energy and the Fock contribution for PCM on a cavity
fixed with respect to the density. The model has already been updated;
``host_potential`` and ``host_fock`` stand for the host's integral routines.
Read grid inputs as model :doc:`fields` and size response buffers from the
current item's fields. Allocate name buffers with ``MOIST_NAME_MAX + 1``
characters.

``moist_get_coupling_request_missing`` reads the current state on every call,
including answers submitted earlier in the same pass.

.. code-block:: c

   moist_coupling cpl = moist_new_coupling(err, model);
   moist_response resp = moist_new_response(err);
   int ngrid;
   moist_get_model_field_int(err, model, "ngrid", &ngrid);
   /* Check err after each call. Allocate xyz[3*ngrid], xi[ngrid], phi[ngrid]. */
   moist_get_model_field_real(err, model, "xyz", xyz);
   moist_get_model_field_real(err, model, "xi0", xi);

   moist_prepare_model_energy(err, model, cpl);
   while (moist_next_coupling_request(err, cpl)) {
       char name[MOIST_NAME_MAX + 1];
       bool missing;
       moist_get_coupling_request_name(err, cpl, name);
       if (strcmp(name, "gaussian_potential") == 0) {
           moist_get_coupling_request_missing(err, cpl, "phi", &missing);
           if (missing) {
               host_potential(ngrid, xyz, xi, phi);
               moist_answer_coupling_request(err, cpl, "phi", phi);
           }
       } else {
           fprintf(stderr, "unsupported request %s\n", name);
           exit(EXIT_FAILURE);
       }
   }
   /* Check err: next also returns false on failure. get_* fails by name on anything missing. */
   double energy = 0.0; /* Or start from an existing contribution. */
   moist_get_model_energy(err, model, cpl, &energy);

   moist_prepare_model_response(err, model, cpl);
   /* Walk the requests again as above: a density-dependent cavity also needs dphi_dr, dphi_dxi. */
   moist_get_model_response(err, model, cpl, resp);
   while (moist_next_response_item(err, resp)) {
       char item[MOIST_NAME_MAX + 1];
       moist_get_response_item_name(err, resp, item);
       if (strcmp(item, "potential_adjoint") == 0) {
           int nfield;
           moist_get_response_field_count(err, resp, &nfield);
           for (int i = 0; i < nfield; ++i) {
               char field[MOIST_FIELD_NAME_MAX + 1];
               int dtype, rank, dims[MOIST_FIELD_MAX_RANK], count;
               moist_get_response_field_info(err, resp, i, field, &dtype, &rank, dims, &count);
               if (strcmp(field, "w_phi") != 0) continue;
               double* w_phi = malloc((size_t)count * sizeof *w_phi);  /* count == ngrid */
               moist_get_response_field_real(err, resp, "w_phi", w_phi);
               host_fock(ngrid, xyz, xi, w_phi);
               free(w_phi);
           }
       } else {
           fprintf(stderr, "unsupported response item %s\n", item);
           exit(EXIT_FAILURE);
       }
   }
   /* Check err: next also returns false on failure. */

   moist_delete_response(&resp);
   moist_delete_coupling(&cpl);

``moist_get_model_gradient`` adds into ``gradient[nat_cap][3]`` and fills the
response with the host part of the gradient phase. For density-dependent
cavities, retain the response-phase density weights as described in
:doc:`coupling`. See ``test/api/example.c`` for executable error-handling
and lifetime examples.

Component energies
~~~~~~~~~~~~~~~~~~

Continuum components use 0-based insertion indices; count and names need no update.
Their real scalar ``energy`` :doc:`field <fields>` stores the latest ``moist_get_model_energy`` contribution after successful evaluation and clears on ``moist_update_model``. Other model families and out-of-range indices are errors.

.. code-block:: c

   int ncomponents;
   moist_get_model_component_count(err, model, &ncomponents);
   for (int i = 0; i < ncomponents; ++i) {
       char name[MOIST_NAME_MAX + 1];
       size_t length;
       double energy;
       moist_get_model_component_name(err, model, i, name, sizeof name, &length);
       moist_get_model_component_field_real(err, model, i, "energy", &energy);
       printf("%s %.10f\n", name, energy);  /* Check err after each call. */
   }

``moist_get_model_component_field_*`` provides ``_count``, ``_info``, ``_about`` and ``_real`` as for model fields. ``moist_get_model_component_description`` copies a one-line description.

Settings printout
~~~~~~~~~~~~~~~~~

``moist_get_model_parameters_text`` returns continuum settings without an update:
the cavity, then components numbered from 1, with descriptions, solvent inputs and registered settings.
Use ``moist_get_banner`` buffer rules to retrieve and print the text.
Other model families are errors.

Contracts
---------

``moist_get_model_energy`` and ``moist_get_model_gradient`` **accumulate** into
caller-owned buffers; initialize them before the first call. Failure leaves the
accumulator unchanged. If a component fails during energy evaluation, only
earlier components' energies remain available.

The full derivative tensors of ``moist_get_cavity_gradient`` and
``moist_get_amat_gradient`` are diagnostic facilities. Normal host integration
uses model gradients and response weights, keeping cavity-motion contractions
inside MOIST. Their dimensions reverse the native axes exactly; consult the
header for each axis's meaning. The ``r_iI1_rA`` and ``rho1_rA`` outputs of
``moist_get_cavity_gradient`` accept ``NULL`` and exist only for a cavity
built with ``do_fine``; asking for one that was not computed is an error.

For internal isodensity models, ``moist_set_model_isodensity_density`` installs
new density data on the owned cavity and invalidates model results and every
coupling of the model. Follow it with ``moist_update_model``. The standalone
``moist_set_isodensity_density`` is restricted to independently owned cavities.
