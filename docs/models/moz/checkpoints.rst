.. _moz-checkpoints:

Checkpoints and initial guesses
===============================

The 1D VV, 1D UV and 3D UV solves can start from a stored checkpoint and save their results while and after they are solved.

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
Any other start from a checkpoint begins the solve; a non-finite or infeasible one is ignored with a warning and the solve starts cold, and a failed warm attempt retries the cold start once, as for in-memory reuse.
The metadata in a state (closure, theory, residual, convergence flag) is provenance only.

Which states seed a solve
-------------------------

A state seeds a solve of the same kind with the same solvent sites, the same solute sites (atomic numbers in order) and the same number of field columns.
Positions and charges may differ, so a field of a displaced or recharged solute, of another closure, or of an interrupted run is a valid starting point.
The grid may differ too, and the source message names the method that moves the stored field onto the grid of the solve (verbosity 2 for 1D VV, 4 for 1D and 3D UV).
A state that ``identity`` applies to is moved by it alone; for any other state the points that coincide with stored ones are copied and the others are interpolated from the stored grid, and the method is ``copy`` when every point coincides.

- ``identity``: molecular grids with the same atomic grids (shell radii and angular point counts per atom) on the built-in Lebedev angular grids.
  Each point is named by its position among the unpruned points (atom, radial shell, angular point), which moves with its atom.
  Pruning depends on the structure, so the stored values go to the points retained at both structures.
- ``copy``: a current point within 1e-9 bohr of exactly one stored point, once the stored solute is fitted onto the current one, takes its value: the same grid, any other grid restarting on itself, the same box around a translated solute, the shared points of a box larger by an even number of points along each axis with the same spacing, or of a longer radial grid.
  Points at two stored points that close are not copied.
- ``fourier`` and ``trilinear``: a stored Cartesian grid seeds any other 3D grid, of another size, spacing or type, or around a rotated solute.
  The stored solute is fitted onto the current one by a rotation and a translation (one atom only translates; collinear atoms are not turned about their axis), and the stored field is evaluated at the current points: on a Cartesian grid through the Fourier modes of the stored box, on a molecular grid by trilinear interpolation between the eight surrounding stored points.
  Only the short-range part is interpolated: the long-range part of the stored state is computed at the new points and added back, so this needs the long-range data of the stored state (``/longrange``).
  Points that coincide with stored ones are copied, the others interpolated.
  A state whose long-range data is unknown or incomplete seeds only the coinciding points (``copy``).
- ``points``: a stored molecular grid, or any grid stored with its points, seeds any other 3D grid by moving least squares.
  After the same rigid fit of the solute, each current point fits an affine function to the nearest stored points (16, else 32 or 64), weighted by their distance; the stencil is accepted when the neighbourhood spans three dimensions and the stencil does not amplify the stored values by more than a factor of 4.
  Points beyond four local spacings of their nearest stored point (the distance to its nearest distinct stored point), and points whose neighbourhood is flat, are not covered.
  Within that reach a point outside the stored points is extrapolated linearly, as long as the stencil does not amplify the stored values by more than a factor of 4.
  Like ``fourier``, this moves the short-range part and needs the long-range data of the stored state.
- ``spline``: radial grids with other nodes; ``t`` is splined through the stored nodes.

Points the stored grid does not cover (outside a smaller box, beyond the stored radial range or the reach of the stored points, or not retained by a molecular state) start cold, from ``c_s = 0`` in 1D and from the selected cold guess (``cs0`` or ``ts0``) in 3D, and the start is checked like any warm guess.
When some points are not covered, the source message gives the number of covered points.
Point positions name the same directions only on the built-in Lebedev angular grids: a molecular state from a grid with another angular generator, or without this record, is moved by ``copy`` and ``points`` instead of ``identity``.

States hold physical fields, ``t = t_s + beta u_lr``, ``c = c_s - beta u_lr`` and ``h = c_s + t_s``.
Loading converts ``t`` with the current long-range potential, so a state stays usable when the Ng split or ``alpha`` changes.

When states are loaded and saved
--------------------------------

- **1D VV**: every solve loads a guess.
- **1D and 3D UV**: a solve loads from the checkpoint when no in-memory field is reused, i.e. on the first solve after construction or ``update``, after a failed or unconverged solve, after a change of the grids or of the bulk solution, and when ``convergence.reuse_initial_guess = false``.
- A store with ``interval > 0`` receives the accepted iterate every ``interval`` iterations on the requested closure, marked unconverged. Fields on the lower closures of the cold-start ramp are not saved. When the solve iterated and ended without an error on the requested closure, its final field is saved, converged or not. A field accepted without iterations is not saved again.
- Failed saves and unreadable or unsuitable stored states produce a warning; the solve proceeds.

Stores
------

``moz_checkpoint_type`` (``moist_model_moz_checkpoint``) is abstract with ``save``, ``load`` and ``describe``; its ``header`` describes the states of the current solve and is set by the owner.
``moz_file_checkpoint_type`` (``moist_model_moz_checkpoint_file``) saves by replacing ``file`` and loads from a list of candidate files, taking the first existing, readable file that can seed the solve; a candidate of the right sites that no method moves onto the grid, or that covers none of its points, is skipped with a warning.
``moz_memory_checkpoint_type`` (``moist_model_moz_checkpoint_memory``) keeps the state in a variable of the caller, the slot.

.. code-block:: fortran

   use moist_model_moz_checkpoint, only: moz_checkpoint_state_type
   use moist_model_moz_checkpoint_memory, only: moz_memory_checkpoint_type, new_memory_checkpoint
   type(moz_checkpoint_state_type), target :: slot
   type(moz_memory_checkpoint_type) :: store

   call new_memory_checkpoint(store, slot)
   call model%set_checkpoint(store)
   do istep = 1, nstep
      call model%update(mol(istep), error)   ! each solve starts from the field in slot
      ...                                    ! and replaces it when it iterates; slot%tr is t(r)
   end do

Models copy the store they are given, and every copy saves into and loads from the same slot.
Each save replaces the slot, so it holds the last saved field with the grid, sites and provenance of the solve that saved it.
A field accepted without iterations is not saved: the slot then still fits the current grid, but its origin and solute coordinates are those of the structure that saved it.
Iterations on the lower closures of the cold-start ramp also leave the slot unchanged.
The caller may also fill the slot before a solve: it then needs ``tr``, the kind, the site lists and the layout of the grid that holds ``tr`` (``get_layout`` of that grid). The field shape is taken from ``tr``.
On a 3D grid every method except ``identity`` also needs the solute coordinates, on the same grid too.
Interpolated points also need the long-range data: ``tail`` ``zero``, or ``known`` with ``alpha``, the temperature and the solute and solvent charges; without it only coinciding points are copied.
A slot without ``tr`` offers nothing.
The slot must outlive the models that hold the store, and the solvent and the UV model need separate slots.

A save creates a temporary file ``<file>.<n>.tmp`` of its own, which no other save can claim, writes it and then replaces ``file`` in one step, also on Windows.
An interrupted write or a failed replacement keeps the previous checkpoint.
Jobs saving the same file at the same time do not write into each other's data; the last replacement is kept.
An interrupted save can leave its temporary file behind; it is never read and can be deleted.
A load checks the metadata and the extents of ``/tr`` before it reads ``/tr``.
Only t(r) is stored: c(r) and h(r) follow from it by the closure.
Malformed containers are refused when opened, so a corrupted file costs a warning, never the run.

Moving a field between grids
----------------------------

``moz_blowup`` (``moist_model_moz_blowup``) moves a field from one grid onto another and keeps the quantity, its units and its columns:

.. code-block:: fortran

   use moist_model_moz_blowup, only: moz_blowup
   logical, allocatable :: covered(:)

   call moz_blowup(grid_old, t_old, grid_new, t_new, error, covered=covered)

The grids are updated 3D grids or radial grids, or their layouts from ``get_layout``.
Each point of the new grid takes the first method that applies: ``copy`` from a point of the old grid within 1e-9 bohr; inside the box of a Cartesian old grid ``fourier`` onto a Cartesian new grid (the field is taken as periodic in that box) and ``trilinear`` onto other grids; ``points`` (moving least squares, as above) from the points of any other 3D grid; ``spline`` up to the last node of a radial old grid.
Points no method covers keep their values in ``t_new`` and are left out of ``covered``; without ``covered`` they are an error and ``t_new`` is unchanged.
Optional arguments: ``map``, a rigid fit (``fit_rigid`` in ``moist_math_grid_transfer``) under which the new points are evaluated in the frame of the old grid; ``options``, to allow only copies, to choose ``fourier`` or ``trilinear`` for every new grid (``cartesian``), to set the threads, or to limit the memory of stored ``points`` stencils (``stencil_memory``, 0.5 GB by default, briefly up to twice that while the stencils are collected; beyond it they are computed again when applied); ``report``, the number of points each method served.

For the correlations of a solve, move ``t`` and rebuild ``c`` and ``h`` with the closure on the new grid.
A restart moves the long-range part of ``t`` separately and fits the solute rigidly; ``moz_checkpoint_guess`` below does both.
To move several fields between the same grids, prepare a ``grid_transfer_type`` once and apply it to each field; prepare again after a grid changes.

Moving a state without a store
------------------------------

``moz_checkpoint_guess`` (``moist_model_moz_checkpoint``) applies the same methods to a state the caller keeps:

.. code-block:: fortran

   call moz_checkpoint_guess(state, header, ctx, beta, ur_lr, guess, matched, method, reason)

``header`` describes the target solve (kind, site lists, ``grid`` from ``get_layout`` of the target grid, solute coordinates), ``beta`` and ``ur_lr`` are its inverse temperature and long-range potential.
``guess`` is the short-range ``t_s`` on the target points and ``matched`` marks the covered ones; ``method`` names the method used.
Without a method ``guess`` stays unallocated and ``reason`` says why.

File layout
-----------

Files are binary containers (``moist_utils_binio``, ``.mbin``), format version 2; files of another version are skipped with a warning.

.. list-table::
   :header-rows: 1

   * - Path
     - Content
   * - ``/checkpoint/version``, ``/checkpoint/kind``
     - format version; ``1D VV``, ``1D UV`` or ``3D UV``
   * - ``/grid/kind``, ``/grid/shape``, ``/grid/spacing``
     - ``cartesian``, ``radial`` or ``points``; ``[nx, ny, nz]``, ``[npts]``, the unpruned
       point count of a molecular grid, or ``[ngrid]``; uniform spacing, 0 without one
   * - ``/grid/origin``
     - Cartesian grids: box corner
   * - ``/grid/r``, ``/grid/weights``
     - radial grids: nodes and their weights
   * - ``/grid/xyz``, ``/grid/weights``, ``/grid/owner``
     - point grids: coordinates ``(3, npoint)``, quadrature weights, owning atom (0 for none)
   * - ``/grid/kspace/modes``, ``/grid/kspace/spacing``, ``/grid/kspace/ref``
     - point grids with a reciprocal grid: mode counts, mode spacings and phase reference
   * - ``/grid/points``, ``/grid/shells/offset``, ``/grid/shells/radius``, ``/grid/shells/points``
     - molecular grids only: position of each point among the unpruned points;
       first shell of each solute atom, shell radii and angular points per shell
   * - ``/grid/angular``
     - molecular grids only: ``lebedev`` for the built-in angular grids
   * - ``/sites/solvent``, ``/sites/solute``, ``/solute/coordinates``
     - atomic numbers; solute coordinates ``(3, nu)`` (UV only)
   * - ``/longrange/status``
     - 3D UV: ``known`` when the stored charges reproduce the long-range potential of the solve,
       ``zero`` without one, else ``unknown``
   * - ``/longrange/alpha``, ``/charges/solute``, ``/charges/solvent``
     - known long-range potential: screening parameter, charges of the solute atoms and of the
       solvent sites
   * - ``/metadata/*``
     - ``closure``, ``theory``, ``temperature``, ``residual`` (last measured),
       ``lambda``, ``iterations``, ``converged``
   * - ``/tr``
     - physical indirect correlation: VV ``(npts, npair)`` packed, 1D UV ``(npts, nu nv)``,
       3D UV ``(ngrid, nv)``
