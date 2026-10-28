.. _moz-checkpoints:

Checkpoints and initial guesses
===============================

The 1D VV, 1D UV and 3D UV solves can start from a stored checkpoint and save theirs result while and after they are solved.

.. code-block:: fortran

   use moist_model_moz_checkpoint_file, only: moz_file_checkpoint_type, new_file_checkpoint
   type(moz_file_checkpoint_type) :: vv_store, uv_store

   call new_file_checkpoint(vv_store, "water-kh.mbin")
   call solvent%set_checkpoint(vv_store)          ! before the model copies the solvent
   call new_moz_3d_model(model, ctx, grid, solvent, error, parameters)
   call new_file_checkpoint(uv_store, "chloride.mbin", interval=10)
   call model%set_checkpoint(uv_store)

Without a store nothing is read or written. The UV models copy the solvent, and the copy
carries its store, so the bulk solve inside ``prepare_solvent`` uses it.

Acceptance rule
---------------

A loaded checkpoint is converged exactly when its residual under the current settings (closure, potential, ``w(k)``, temperature and bulk theory) is below the tolerance. 
Any other checkpoint that fits the discretization starts the solve; a failed warm attempt retries the cold start once, as for in-memory reuse.
The metadata in a state (closure, theory, temperature, residual, convergence flag) is provenance only.

A state fits a solve when its kind, grid shape, grid spacing, solvent sites and solute sites (atomic numbers in order) agree and its fields have the shape the solve needs.
Positions and charges may differ, so a field of a displaced or recharged solute, of another closure, or of an interrupted run is a valid starting point.
The Cartesian grid origin follows the solute and is not compared.

States hold physical fields, ``t = t_s + beta u_lr``, ``c = c_s - beta u_lr`` and ``h = c_s + t_s``.
Loading converts ``t`` with the current long-range potential, so a state stays usable when the Ng split or ``alpha`` changes.

When states are loaded and saved
--------------------------------

- **1D VV**: every solve loads a guess.
- **1D and 3D UV**: a solve loads from the checkpoint when no in-memory field is reused, i.e. on the first solve after construction or ``update`` and when ``convergence.reuse_initial_guess = false``.
- A store with ``interval > 0`` receives the accepted iterate every ``interval`` iterations, marked unconverged. When the solve iterated and ended without an error, its final field is saved, converged or not. A field accepted without iterations is not saved again.
- Failed saves and unreadable or unsuitable stored states produce a warning; the solve proceeds.

Stores
------

``moz_checkpoint_type`` (``moist_model_moz_checkpoint``) is abstract with ``save``, ``load`` and ``describe``; its ``header`` describes the states of the current solve and is set by the owner.
``moz_file_checkpoint_type`` (``moist_model_moz_checkpoint_file``) saves by replacing ``file`` and loads from a list of candidate files, taking the first existing, readable file that fits.

A save writes ``<file>.tmp`` and renames it over ``file``, so an interrupted or failed save keeps the previous checkpoint.
A load checks the metadata and the extents of ``/tr`` before it reads ``/tr``.
Only t(r) is stored: c(r) and h(r) follow from it by the closure.
Malformed containers are refused when opened, so a corrupted file costs a warning, never the run.

File layout
-----------

Files are binary containers (``moist_utils_binio``, ``.mbin``), format version 1.

.. list-table::
   :header-rows: 1

   * - Path
     - Content
   * - ``/checkpoint/version``, ``/checkpoint/kind``
     - format version; ``1D VV``, ``1D UV`` or ``3D UV``
   * - ``/grid/shape``, ``/grid/spacing``, ``/grid/origin``
     - ``[npts]``, ``[nx, ny, nz]`` or ``[ngrid]``; spacing, 0 without one;
       Cartesian box corner (optional)
   * - ``/sites/solvent``, ``/sites/solute``, ``/solute/coordinates``
     - atomic numbers; solute coordinates ``(3, nu)`` (UV only)
   * - ``/metadata/*``
     - ``closure``, ``theory``, ``temperature``, ``residual`` (last measured),
       ``lambda``, ``iterations``, ``converged``
   * - ``/tr``
     - physical indirect correlation: VV ``(npts, npair)`` packed, 1D UV ``(npts, nu nv)``,
       3D UV ``(ngrid, nv)``
