.. _api-context:

Run context
===========

The context holds the settings that every part of a calculation shares: the thread count, the output level and the timer.
A model requires one when it is constructed.
A cavity or component may be given its own, which it keeps inside the model; one without runs on the model's.
A cavity used on its own, outside a model, needs a context.
The context must outlive everything built on it: the C and Python objects retain it, so the host may drop its own handle early, while a Fortran host keeps its ``target`` variable alive.

Threads
-------

The thread count is fixed when the context is created: a positive count is used as given and ``0`` takes the calling thread's OpenMP setting at that moment.
It sizes MOIST's own OpenMP regions and FFT workers.
MOIST never changes the host's OpenMP runtime; BLAS/LAPACK threading remains the host's to configure (``OMP_NUM_THREADS``, ``MKL_NUM_THREADS``, ``OPENBLAS_NUM_THREADS`` or the vendor's API).
Without OpenMP, MOIST's own kernels are serial; the math library may still use threads.

Calls sharing a context must not run at the same time.
Independent calculations, each with its own context, may run on separate threads.

.. _verbosity:

Output
------

The verbosity of the context sets how much MOIST prints.
A host calls MOIST in every SCF iteration and geometry step, so each level also bounds how often it prints:

== ============== ===========================================================
0  silent         nothing
1  minimal        header block and results, once per run
2  normal         settings relevant to a calculation, a summary per geometry
3  detailed       every setting and per-geometry tables
4  solver summary one block per inner solve, e.g. the DROP projection
5  solver trace   per-iteration and per-point lines of inner solves
6  debug          everything; the debug stream itself needs ``debug``
== ============== ===========================================================

Level 2 is the default of the Fortran, C and Python context and suits a production calculation.
It prints nothing inside the host's SCF iterations, only once per run and once per geometry.
It shows what is needed to report and reproduce a result, such as the radii, the solvent potential of a MOZ model or the PCM solver, and leaves out numerical internals such as the DROP switching parameters.

Messages go to standard output; a Fortran host can redirect them to a unit or a file.
Warnings go to the same output as the messages at every verbosity, including 0; ``warnings_to_stderr`` also copies them to standard error.
The ``debug`` flag raises the verbosity to 6 and enables an additional debug stream; verbosity 6 alone does not enable it.
Fine-grained profiling timers are a separate switch, independent of the verbosity (Fortran ``do_profile``, the command line's ``--profile``).
The command-line program defaults to verbosity 4.

MOIST prints its settings and timings only when the host asks:

- The settings printout lists the cavity and then each component: nothing at verbosity 0, the normal settings at 1 and 2, every setting from 3.
- The radii table of the current structure has one row per element, or a short summary for custom radii: nothing at verbosity 0, the table from 1.
  It needs an update and is ready before a deferred build: Fortran ``model%print_radii()`` prints nothing before the update, C ``moist_get_model_radii_text`` and Python ``model.radii_text()`` fail.
- The timer prints nothing at verbosity 0 and the top-level phases at 1, adds one level of the tree per step up to 3, and prints the whole tree from 4.
- The banner is returned as text for the host to print: Fortran ``moist_banner_text(style)``, C ``moist_get_banner``, Python ``get_banner(style)``.

.. tab-set::

   .. tab-item:: Fortran
      :sync: fortran

      .. code-block:: fortran

         use moist_context, only : moist_context_type, new_context

         type(moist_context_type), target :: ctx

         call new_context(ctx, nthreads=0, verbosity=2)
         ! ... construct the model on ctx ...
         call model%print_parameters()
         ! ... run the calculation ...
         call ctx%timer%print()  ! or ctx%timer%print_simple()

   .. tab-item:: C
      :sync: c

      .. code-block:: c

         moist_context context = moist_new_context(error, 0, 2, false, false);  /* nthreads, verbosity, debug, warnings_to_stderr */
         /* ... construct the model on context; keep this handle to print the timings ... */
         size_t length;
         moist_get_model_parameters_text(error, model, NULL, 0, &length);
         char* text = malloc(length + 1);
         moist_get_model_parameters_text(error, model, text, length + 1, &length);
         fputs(text, stdout);
         free(text);
         /* ... run the calculation ... */
         moist_print_context_timer(error, context);  /* or moist_print_context_timer_simple */

   .. tab-item:: Python
      :sync: python

      .. code-block:: python

         from moist import Context

         context = Context(nthreads=0, verbosity=2)
         # ... construct the model on context ...
         model.print_parameters()
         # ... run the calculation ...
         context.print_timer()  # or context.print_timer(simple=True)

The library default is verbosity 2 in every language; each language page lists the remaining context options.

Host timers
-----------

Host timers collect coupling and response work under a separate ``Host`` node in the timing report.
PySCF records this work automatically.

- Time only host work; MOIST calls inside a host timer are counted in both groups.
- Names contain 1 to 32 UTF-8 bytes, no control characters and no leading or trailing spaces. In Fortran, use ``trim(buffer)`` for padded names.
- Stop nested timers innermost first. Repeated names under the same parent accumulate.
- Stop timers on failure paths too. In C and Fortran, use a separate error variable for stopping to preserve the original failure.
- Serialize timer calls with other calls on the context; never call them from a callback.

Start and stop timers through the context:

- Fortran: ``ctx%start_timer(name, error)`` and ``ctx%stop_timer(name, diagnostic)``.
- C/C++: ``moist_start_context_timer(error, context, name)`` and ``moist_stop_context_timer(diagnostic, context, name)``.
- Python: prefer the scoped form, which stops the timer even when the block raises:

  .. code-block:: python

     with context.timer("Fock build"):
         ...  # contract the response

Print the report after closing the timers. Open timers appear as ``running`` and keep measuring.
