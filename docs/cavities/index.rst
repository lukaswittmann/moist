Cavities
========

This section describes the cavity constructions available in MOIST and the settings that control their shape and surface discretization.

Overview
--------

For most applications, one of three cavities is the right choice:

- :doc:`iSwiG <iswig>` for a van der Waals cavity from atomic spheres;
- :doc:`SvdW-DROP <svdw>` for a smooth molecular surface without crevices;
- :doc:`ρ-DROP <isodensity>` for a cavity that follows the electron density of the QM host.

The remaining constructions serve special purposes, such as reference values for areas and volumes.

.. table::
   :class: matrix

   +---------------------------------+-------------+--------------------+------------------+--------------------------------+--------------+----------+---------+
   | Construction                    | Surfaces                                                                             | Output       | Derivatives        |
   |                                 +-------------+--------------------+------------------+--------------------------------+              +----------+---------+
   |                                 | vdW spheres | :doc:`SvdW <svdw>` | :doc:`CFC <cfc>` | :doc:`Isodensity <isodensity>` |              | Gradient | Hessian |
   +=================================+=============+====================+==================+================================+==============+==========+=========+
   | :doc:`iSwiG <iswig>`            | |yes|       | |no|               | |no|             | |no|                           | surface grid | |yes|    | |soon|  |
   +---------------------------------+-------------+--------------------+------------------+--------------------------------+--------------+----------+---------+
   | :doc:`DROP <drop>`              | |no|        | |yes|              | |yes|            | |yes|                          | surface grid | |yes|    | |soon|  |
   +---------------------------------+-------------+--------------------+------------------+--------------------------------+--------------+----------+---------+
   | :doc:`Marching tetrahedra <mt>` | |no|        | |yes|              | |yes|            | |part|                         | surface grid | |no|     | |no|    |
   +---------------------------------+-------------+--------------------+------------------+--------------------------------+--------------+----------+---------+
   | :doc:`Marching cubes <mc>`      | |no|        | |yes|              | |yes|            | |part|                         | area, volume | |no|     | |no|    |
   +---------------------------------+-------------+--------------------+------------------+--------------------------------+--------------+----------+---------+
   | :doc:`NumSA <numsa>`            | |yes|       | |no|               | |no|             | |no|                           | area         | |part|   | |no|    |
   +---------------------------------+-------------+--------------------+------------------+--------------------------------+--------------+----------+---------+

|yes| supported, |part| restricted, |no| not supported, |soon| in development.
Marching tetrahedra and marching cubes accept an isodensity level set for a frozen density only:
they provide no density response, so they cannot be used in an SCF.
NumSA provides the gradient of the total area only.

Common Interface
----------------

All cavity implementations extend the abstract ``cavity_type`` in ``src/moist/cavity/type.f90``.
It carries the shared allocatable quantities that concrete cavities fill where applicable:

- atomic sphere centers, radii, and per-sphere areas;
- grid point positions, owners, areas, outward unit normals, Gaussian widths,
  switching function, and volume elements;
- total cavity area (bohr\ :sup:`2`) and enclosed volume (bohr\ :sup:`3`);
- optional nuclear derivatives of positions, widths, switching function, and
  volume elements.

Concrete cavity types must implement two deferred procedures:

``update_surface(mol, error)``
   Prepare the cavity for the supplied structure, called through ``update``.
   Cavities needing host inputs finish through ``build_surface`` (:ref:`coupling-build`).

``get_gradient(error)``
   Prepare or validate the implementation-specific nuclear derivatives.

The base type also provides optional response hooks; ``get_surface_response``
maps model surface weights to host-response contributions.

Typical Use
-----------

1. Construct a concrete cavity and configure its radii and discretization.
2. Call ``update`` after geometry changes and follow the :ref:`coupling-build`.
3. A ``model_continuum_type`` owns one cavity, updates it before its
   components, and passes the same cavity to each component.
4. In the response phase, components accumulate model surface weights and
   ``get_surface_response`` maps them to the response items the host receives
   (the ``density`` item, for a cavity whose surface follows the density).
5. Call ``get_gradient`` before requesting supported nuclear derivatives.

The shared ``print``, ``write_xyz_debug``, ``write_csv_debug``, and ``write_pqr_debug`` procedures provide diagnostics and grid export.
``find_disconnected_cavities`` from ``moist_cavity_diagnostic`` is a separate diagnostic: it returns the point count of each connected island, largest first, with a connectivity radius of ``threshold`` times the mean nearest-neighbour spacing (4.0 when omitted). No cavity runs it during ``update``.

.. _cavity-radii:

Radius Models
-------------

Every cavity needs a radius model. ``CPCM`` radii are the default; ``SMD``,
``COSMO``, ``Bondi``, ``D3``, ``Gauss`` and custom per-atom or per-element radii
are also available. ``Gauss`` uses Gaussian charge scheme radii, including
their 1.2 scale factor. The cavity constructor copies the model.

.. tab-set::

   .. tab-item:: Fortran
      :sync: fortran

      .. code-block:: fortran

         use moist_radii, only : radius_type_static, new_cpcm_radii

         type(radius_type_static) :: radii

         call new_cpcm_radii(radii)

   .. tab-item:: C
      :sync: c

      .. code-block:: c

         moist_radii radii = moist_new_cpcm_radii(error);
         /* cavity constructors copy it; moist_delete(radii) afterwards.
            Passing NULL radii selects CPCM radii. */

   .. tab-item:: Python
      :sync: python

      .. code-block:: python

         from moist import CPCMRadii

         radii = CPCMRadii()

Per-atom and per-element radii are supported too, as ``CustomRadii`` in Python,
``moist_new_custom_radii`` plus its setters in C, and ``new_radii_custom_atoms`` /
``new_radii_custom_elements`` in Fortran; see :doc:`/reference/fortran`.

.. _cavity-symmetry:

Symmetry and Grid Orientation
-----------------------------

A Lebedev grid has the symmetry of a cube, not of a sphere, so areas, energies and gradients depend on how the grid of each atomic sphere is oriented.
In the *lab frame* every grid follows the input coordinate axes: a rotated molecule gives slightly different results, and a symmetric molecule a slightly asymmetric surface.
A symmetry model orients the grids by the molecule instead:

Molecular orientation
   Every grid follows a molecular frame computed from the atom positions.
   Results are invariant under rotation, translation and renumbering of the atoms, except for tied frames (see below).

Point group
   The grids of symmetry-equivalent atoms are mapped onto each other by the symmetry operations, so the surface carries the point group exactly (to round-off).
   A group always uses the molecular orientation.

The model serves the Lebedev-based cavities: :doc:`iSwiG <iswig>`, :doc:`DROP <drop>` with any level set, and :doc:`NumSA <numsa>`.
For :doc:`ρ-DROP <isodensity>` the surface is only as symmetric as the density the host supplies.

Symmetry is opt-in: a cavity constructed without a symmetry model keeps the lab frame, and its results are unchanged.
The command line always builds a symmetry model from its flags.

Molecular Frame
~~~~~~~~~~~~~~~

The frame is built from the orientation (via the shape) of the molecule: the second moments of the atom positions about their centroid (an inertia tensor with every atom weighted equally).

- Asymmetric molecule (three different eigenvalues): the three eigenvectors are the axes.
- Symmetric top (two equal eigenvalues, such as ammonia or benzene): the eigenvectors give the main axis only.
  The angle about it is taken from how the atoms are arranged around the axis.
- Spherical top (three equal eigenvalues, such as methane): the axes are taken from higher moments of the positions.
- Linear molecule: only the axis is defined.
  The angle about it follows the coordinate (lab frame) axes.

Where the symmetry leaves several equivalent choices (the three equivalent angles of ammonia), they give the same grid, or grids that a symmetry operation of the molecule maps onto each other, so the results are the same.

.. warning::

   Each of these choices is switched at a threshold (see the settings).
   Near a threshold a small displacement can switch a choice, and areas, energies and gradients jump with the grid.
   This can disturb geometry optimizations, finite differences and scans.
   The cases, each with its :ref:`field <fields>`:

   - nearly degenerate eigenvalues of the frame tensor: a nearly symmetric or nearly spherical top (``symmetry_near_gap``);
   - a nearly linear molecule, atoms within about 0.003 to 0.1 bohr RMS of one line (``symmetry_near_linear``);
   - a weak in-plane harmonic, or two harmonics or equivalent angles of nearly the same quality (``symmetry_near_harmonic``).

   A measure between 0.1 and 100 times its threshold counts as near: the model prints a warning, adds a *Near a switch* line to its report and sets the field.
   An exactly symmetric molecule gives no warning and, unless its point group is cubic (see below), is far from the thresholds.

   What to do: symmetrize the geometry and give its point group, use the lab frame (which has no switches), or move the threshold away from the geometries of interest.

Tied frames
   The frame uses the atom positions only, not the elements or radii.
   In rare cases the positions are more symmetric than the molecule, for example a regular five-membered ring with one heteroatom, or an icosahedral cage with one substituted corner.
   Several frames are then equally good, the one used depends on the input orientation, and so do the results, as in the lab frame (0.2 to 1 % of the area in these examples).
   The model prints a warning, adds a *Tied frames* line to its report and sets ``symmetry_tied_frames``.
   A point group does not help in these examples; common symmetric molecules such as benzene or methane are not affected.

Other jumps and orientation dependence
   - With ``group = "auto"``, a geometry that drifts out of the ``tolerance`` of its group gets a smaller group, and the results jump.
     A Schoenflies symbol keeps the group fixed and stops with an error instead.
   - A group does not remove the switches of the molecular frame described above.
   - With a subgroup of the molecule's point group (C\ :sub:`6h` for benzene), the grids can be placed in several equivalent ways, and the input orientation decides which.
     The surfaces are symmetry copies of each other, so results differ only in an environment that does not share the symmetry of the molecule.
   - A slightly bent molecule that is still taken as linear depends slightly on its orientation.

Cubic point groups
   A molecule with a cubic point group (methane, SF\ :sub:`6`) sits exactly on a tie: its axis is one of several equivalent maxima of a higher moment, and a displacement that lowers the symmetry picks one.
   Energies and gradients are continuous there and the analytic gradient is exact, but central finite differences converge only linearly in the step, and second derivatives depend on the direction of the displacement.
   Finite differences must also stay below the eigenvalue switch, which a hydrogen of methane reaches at about 2e-4 bohr.
   With a centre of inversion (SF\ :sub:`6`, a cube) the axis switches from the fourth to the third moment at a distortion of about 1e-5 bohr, which turns the grids by a finite angle and changes areas by a few percent; this switch has no setting and no warning.
   Give the point group with a ``tolerance`` that covers the noise of the geometry, or use the lab frame.

Settings
~~~~~~~~

The settings form a top-level ``symmetry`` table, beside (not inside) the settings of a cavity.
In Fortran they are the fields of ``moist_cavity_grid_model_symmetry_parameters_type``, which reads and writes the table with ``read_file`` and ``write_file`` (JSON or TOML).

.. code-block:: toml

   [symmetry]
   orientation = "molecular"
   group = "auto"
   tolerance = 1.0e-8

``orientation`` (string, ``"auto"``, ``"lab"`` or ``"molecular"``, default ``"auto"``)
   ``"auto"`` resolves to ``"molecular"`` with a group and to ``"lab"`` without one.
   ``"lab"`` together with a group other than ``"off"`` is an error.

``group`` (string, default ``"off"``)
   ``"off"``: no point group.
   ``"auto"``: the point group is detected for every geometry.
   A Schoenflies symbol such as ``"C2v"`` (case-insensitive): this group, which must be the detected group or one of its subgroups.
   ``"C1"`` gives the molecular orientation without a point group.
   For a linear molecule the symbol may name the linear group (``"D0h"`` or ``"Dinfh"``, ``"C0v"`` or ``"Cinfv"``) or a subgroup of D\ :sub:`4h` or C\ :sub:`4v`, which stand in for it.

``tolerance`` (real, default ``1.0e-8``)
   Largest distance in bohr by which a symmetry operation may move an atom away from an equivalent atom.
   Atoms are equivalent only with the same element and the same sphere radius.
   The geometry is never symmetrized.

``gap_threshold`` (real, default ``1.0e-4``)
   Gap between two neighbouring eigenvalues of the frame tensor, relative to their sum, above which the two count as different.
   The frame tensor is the second-moment tensor, joined by the fourth moments where a pair is not separated.

``harmonic_threshold`` (real, default ``1.0e-4``)
   Size of an in-plane harmonic, relative to its largest possible value, above which it can fix the angle about the axis.

``linear_threshold`` (real, default ``1.0e-4``)
   Mean square distance of the atoms from the long axis, in bohr\ :sup:`2`, below which the molecule is taken as linear (0.01 bohr RMS by default).
   It is an absolute measure, so a long molecule of finite thickness is never taken as a line.

Point Group
~~~~~~~~~~~

The point group is detected with libmsym, and every operation is then checked against ``tolerance``.
libmsym's own equivalence threshold follows ``tolerance`` (ten times its value, between 1e-10 and 5e-4).
A failed detection is repeated with tighter thresholds, and operations that miss the atoms are first fitted to them by one rigid rotation; the reported deviation is that of the fitted operations.

- With ``"auto"``, a detected group that the geometry misses by more than ``tolerance`` is replaced by its largest subgroup that passes, down to C\ :sub:`1`; a warning names the deviation.
- With a Schoenflies symbol, a group that is not a subgroup of the detected group, or that the geometry misses by more than ``tolerance``, stops the update with an error.

A grid can only follow the operations that keep its atom in place (the site group) if the site group is a subgroup of O\ :sub:`h`, the symmetry of a cube.
The group used is therefore the largest subgroup of the detected or requested group in which every site group qualifies.
Linear molecules use the largest qualifying subgroup of D\ :sub:`4h` or C\ :sub:`4v`.

.. list-table::
   :header-rows: 1
   :widths: 20 15 15 50

   * - Molecule
     - Detected
     - Used
     - Reason
   * - water
     - C\ :sub:`2v`
     - C\ :sub:`2v`
     - every site group qualifies
   * - BF\ :sub:`3`
     - D\ :sub:`3h`
     - C\ :sub:`3v`
     - the site group of boron, D\ :sub:`3h`, has an S\ :sub:`3` axis
   * - CO\ :sub:`2`
     - D\ :sub:`∞h`
     - D\ :sub:`4h`
     - linear

Usage
~~~~~

.. tab-set::

   .. tab-item:: Fortran
      :sync: fortran

      .. code-block:: fortran

         use moist_cavity_grid, only : grid_model_type, new_grid_model, &
            & moist_cavity_grid_model_symmetry_parameters_type

         type(moist_cavity_grid_model_symmetry_parameters_type) :: sym_param
         class(grid_model_type), allocatable :: grid_model

         sym_param%group = "auto"
         call new_grid_model(grid_model, error, symmetry=sym_param)
         if (allocated(error)) error stop error%message

         call new_cavity_iswig(cavity, radius_model=radii, error=error, ctx=ctx, &
            & grid_model=grid_model)
         if (allocated(error)) error stop error%message

      ``new_cavity_drop`` and ``new_cavity_numsa`` take the same optional ``grid_model`` argument.
      The cavity keeps a copy of the model and updates it in every ``update``; ``print_radii`` prints the report after the radii table.

   .. tab-item:: Command line
      :sync: cli

      .. code-block:: console

         moist cavity iswig water.xyz --sym auto
         moist cavity drop svdw water.xyz --orientation molecular
         moist cavity numsa water.xyz --sym C2v --sym-tolerance 1e-6

      ``--sym GROUP``, ``--orientation ORIENT`` and ``--sym-tolerance BOHR`` set ``group``, ``orientation`` and ``tolerance``.
      They are available for ``numsa``, ``iswig`` and ``drop svdw``, ``drop cfc`` and ``drop sses``.

.. _cavity-symmetry-report:

Report
~~~~~~

The report follows the radii table, here for water with ``--sym auto``:

.. code-block:: none

   Symmetry:
     Orientation               ... molecular (implied by group)
     Group requested           ... auto
     Group detected            ... C2v (deviation 1.75E-16 bohr)
     Group used                ... C2v (order 4)
     Max. deviation            ... 1.75E-16 bohr (tolerance 1.00E-08 bohr)
     Molecular frame           ... eigen (eigen-gap 6.67E-01)
     Orbits                    ... 2
     Grid rotations            ... 3

   Symmetry-equivalent atoms:
    ------ ---- ---------- ----- ----- ---------- --------------------
     Orbit  Sym   R (bohr)  Mult  Site  Rotations                Atoms
    ------ ---- ---------- ----- ----- ---------- --------------------
         1    O     3.2503     1     4          1                    1
         2    H     2.4566     2     2        2,3                  2,3

- *Group detected*: with the largest atom deviation under its operations; linear groups print as ``D0h`` and ``C0v``.
- *Group used*: a *Reason* line explains a group that differs from the detected or requested one.
- *Molecular frame*: the level that fixed the axes (``eigen``, ``plane``, ``sphere``, ``axial``) and the smallest relative gap of two neighbouring eigenvalues of the frame tensor.
- *Orbits*: sets of symmetry-equivalent atoms; *Grid rotations*: distinct grid orientations.
- The table has one row per orbit: number of atoms (*Mult*), order of the site group (*Site*), grid rotations and atoms.

The work appears as ``Symmetry`` in the timings.
A cavity with a symmetry model lists the numeric results among its :ref:`fields <fields>`:

- ``symmetry_orientation``: ``0`` lab, ``1`` molecular.
- ``symmetry_frame_level``: ``0`` none, ``1`` eigen, ``2`` plane, ``3`` sphere, ``4`` axial.
- ``symmetry_near_gap``, ``symmetry_near_linear``, ``symmetry_near_harmonic``, ``symmetry_tied_frames``: one truth value each.
- ``symmetry_order``, ``symmetry_deviation``, ``symmetry_frame_gap``, ``symmetry_nrot``, ``symmetry_rot_index``, ``symmetry_orbit``, ``symmetry_rotation``: order of the group used, its largest atom deviation, the eigen-gap, and the grid rotations and orbits of the spheres.

Gradients
~~~~~~~~~

The grids turn with the molecule, so moving any atom turns the grids of all atoms.
Analytic nuclear derivatives include this orientation term.
The reverse-mode gradient (``get_surface_gradient``) collects one torque per sphere in its existing sweep, at a few percent extra cost.
The forward-mode arrays (``get_gradient``) become dense in the atoms, at a cost proportional to points times atoms.

With a group of more than one operation, the orientation term is the derivative of the cavity whose group moves rigidly with the molecule.
Along the totally symmetric part of a displacement the grids follow the molecular frame within what each site allows; of the rest, only the rigid rotation counts.
The derivatives are therefore exact along rigid motions and along displacements that keep the group, for every quantity.
Along a symmetry-breaking distortion the group itself is lost, and the derivative leaves out the orientation response to it.

Limitations
~~~~~~~~~~~

- Results can change abruptly where a choice of the molecular frame or of the group is switched (see the warning above).
- Linear molecules carry D\ :sub:`4h` or C\ :sub:`4v`, not the full linear group.
- DROP: the octree seeding of the projection (levels 3 and 4, and the retry of level 1) still follows the coordinate axes.
  Results are rotation invariant whenever the seeding finds every projection branch.
- The C and Python interfaces do not expose the symmetry model yet.

Available Cavities
------------------

.. toctree::

   drop
   iswig
   mt

Area & Volume Integrators
-------------------------

These constructions produce no surface grid and cannot carry a solvation model.
They integrate the area (and volume) of a cavity, as reference values or for area-based terms.

.. toctree::

   mc
   numsa
