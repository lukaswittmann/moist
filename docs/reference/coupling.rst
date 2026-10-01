.. _coupling:

Host coupling protocol
======================

A quantum-mechanical host and MOIST share the evaluation of the environment's energy and its derivatives.
The host owns the electronic-structure calculation, including its density matrix, basis functions and integrals.
MOIST owns the cavity and the model components, such as PCM, PV and GOSTSHYP.

The exchange has three parts:

#. **Host to MOIST:** MOIST requests quantities on its cavity grid, such as the solute's electrostatic potential or    Gaussian moments.
   The host evaluates the requested quantities for its current electronic state and supplies them without multiplying by MOIST's charges or amplitudes.
#. **Inside MOIST:** the model evaluates its components using those quantities and the cavity.
   It computes the energy contribution and the coefficients needed to differentiate it, including the dependence on cavity geometry.
#. **MOIST to host:** MOIST returns energy contributions and *energy-derivative coefficients*.
   The host contracts these coefficients with its own integrals or derivatives to assemble the component's contribution to the Fock matrix and nuclear gradient.
   MOIST also adds the part of the nuclear gradient it computes internally.

The API calls the returned energy-derivative coefficients a **response**.
A **coupling** holds the requests and the host's answers for a model evaluation; a **response** holds the coefficients returned to the host.
The same exchange is used in Fortran, C and Python; the names and usage is kept very similar as much as possible.
Language-specific calls and examples are in :ref:`Fortran <coupling-fortran>`, :ref:`C <coupling-c>` and :ref:`Python <coupling-python>`.

Energy derivatives and response
-------------------------------

For a host quantity :math:`q_i`, MOIST supplies its energy derivative :math:`w_i = \partial E_\mathrm{MOIST}/\partial q_i`.
The host supplies the other factor in the chain rule.
For density-matrix element :math:`P_{\mu\nu}` and nuclear coordinate :math:`R_A`, the corresponding contributions are

.. math::

   \Delta F_{\mu\nu}^{(q)}
      = \sum_i w_i \frac{\partial q_i}{\partial P_{\mu\nu}},
   \qquad
   \Delta g_A^{(q)}
      = \sum_i w_i \frac{\partial q_i}{\partial R_A}.

The index :math:`i` includes the grid points and any vector or tensor components.
These expressions describe the host contractions: the host takes its derivatives at fixed cavity geometry.
MOIST handles the cavity dependence and supplies any additional coefficients needed when the cavity follows the density.

For example, PCM requests the total nuclear plus electronic potential :math:`\phi_i` on the cavity grid, evaluated with the Gaussian Coulomb operator specified below.
MOIST solves the polarization problem and returns ``potential_adjoint``, containing :math:`w_{\phi,i} = \partial E_\mathrm{PCM}/\partial\phi_i`.
For stationary CPCM and COSMO, these weights are the induced surface charges.
The host contracts them with :math:`\partial\phi_i/\partial P_{\mu\nu}` to obtain the PCM Fock contribution, and with the explicit nuclear-coordinate derivatives of the potential for the host part of the gradient.
Those derivatives include the electronic sign of the potential; the weights require no extra charge factor.

Typical SCF and gradient workflow
---------------------------------

Construct the cavity and model components, update the model for the molecular structure, and create a coupling. That coupling can be reused throughout the SCF calculation and the subsequent gradient evaluation.

At each SCF density:

#. **Update a density-dependent cavity.** Supply the current density to its    density source and update the model before evaluating it.
   A cavity determined only by nuclear geometry can be reused while the geometry is unchanged.
#. **Evaluate the energy.** Prepare the energy phase, answer its requests and obtain MOIST's energy contribution.
   Preparing this phase discards answers from the previous density.
#. **Build the Fock contribution.** Prepare the response phase, answer any additional requests and obtain the response coefficients.
   Contract every response item with the corresponding density-matrix derivatives and add the results to the host Fock matrix.

For the nuclear gradient, use the model and host quantities at the converged density and current geometry.
Prepare the gradient phase, answer its additional requests and evaluate it.
MOIST adds its internally computed contribution to the gradient accumulator and returns response items for the host's remaining contractions.

On a density-dependent cavity, also retain the ``density`` item from the response phase at that same density and geometry.
Contract its weights with the basis-center derivatives of the density and its spatial derivatives, and add this contribution to the nuclear gradient.
The gradient phase does not return that item again. Save those weights before replacing the response with the gradient-phase result.
If no response has been evaluated for the current state, evaluate it first.

Protocol outline
----------------

Each phase follows the same request-and-answer loop. Here, ``prepare_*`` and ``get_*`` denote the matching operations for the selected phase:

.. code-block:: text

   update model
   create coupling

   for each evaluation phase:
       prepare phase
       while a request with missing outputs is available:
           read the request and its grid inputs
           compute and answer each missing output
       evaluate phase
       if the phase returns a response:
           for each response item:
               contract its coefficients with the host derivatives
               stop if the item is unsupported

Energy evaluation adds to an energy accumulator.
Response evaluation returns coefficients for the Fock contribution.
Gradient evaluation both adds to a gradient accumulator and returns coefficients for host gradient contractions.
The density weights retained from the response phase supply the additional contraction described above.

Grid inputs
-----------

Grid inputs are :doc:`fields` of the model: ``xyz``, ``xi0`` and the other per-point arrays of its evaluation domain.
The coupling keeps no copy of the grid and has no grid accessor of its own. A model update invalidates every coupling, so a staged coupling always refers to the current grid.

.. _coupling-requests:

Requests and outputs
--------------------

Requests, outputs and response items are addressed by their name in every language.
The table uses Fortran shapes; C and Python reverse the dimensions, so ``dphi_dr`` is ``(ngrid, 3)`` and ``mt`` is ``(ngrid, 3, 3)``.

.. list-table:: Requests, their grid inputs and their outputs
   :header-rows: 1

   * - Request
     - Cavity inputs
     - Host outputs
   * - ``point_potential``
     - ``xyz``
     - ``phi(ngrid)``, ``dphi_dr(3, ngrid)``
   * - ``gaussian_potential``
     - ``xyz``, ``xi0``
     - ``phi``, ``dphi_dr``, ``dphi_dxi(ngrid)``
   * - ``gaussian_moments``
     - ``xyz`` and the request's own ``width(ngrid)``
     - ``gt(ngrid)``, ``pt(3, ngrid)``, ``mt(3, 3, ngrid)``, ``rt(3, ngrid)``

Both potential requests carry the total nuclear plus electronic potential.
``dphi_dr`` is its derivative with respect to the evaluation point, without a charge multiplication or electric-field sign change.

PCM uses the Gaussian operator ``erf(xi*r)/r`` consistently in the potential, Fock matrix and nuclear gradient. ``dphi_dxi`` differentiates the inverse length ``xi`` in bohr**-1, not the exponent ``xi**2``.
The point operator is a separate request for components requiring bare Coulomb potentials.

Gaussian moment widths are exponents in bohr**-2, chosen by GOSTSHYP as ``pi*ln(2)/a`` from the per-point area.
They are a :doc:`field <fields>` of the request, not of the cavity: read them from the request, never recompute them.

PCM requires ``phi`` in every phase.
It requires both derivatives for the gradient, and for the response only when the cavity geometry follows the density.
GOSTSHYP requires ``gt`` and ``pt`` in every phase; ``mt`` and ``rt`` for the gradient, and for the response only on a density-dependent cavity.
Zero pressure or zero scale registers no moment request.
PV registers no request (as it requires none).

.. _coupling-response:

Response items
--------------

A response is the list of items one phase hands back, accumulated over the components of the model.
The host contracts each with its own derivative of the quantity it is conjugate to.

.. list-table:: Response items
   :header-rows: 1

   * - Item
     - Arrays
     - Present when
   * - ``potential_adjoint``
     - ``w_phi(ngrid)``
     - Any electrostatic component is present
   * - ``density``
     - ``w_rho(ngrid)``, ``w_grad_rho(3, ngrid)``, ``w_hess_rho(3, 3, ngrid)``
     - The cavity surface follows the host density (the isodensity level sets),
       response phase only
   * - ``gaussian_amplitude``
     - ``w_overlap(ngrid)``, ``w_normal_deriv(ngrid)``
     - A GOSTSHYP component is present

The potential adjoint ``w_phi = dE/dphi`` is the induced surface charge of a stationary PCM (CPCM, COSMO); a non-symmetric response matrix (IEF-PCM, SS(V)PE) emits the symmetrized adjoint instead, which the host contracts identically.
Density adjoints already include the level-set scale; do not apply it again.
A geometric cavity produces no ``density`` item at all: that absence is the physical answer and is distinguishable from an item of zeros.

``density`` is a response-phase item; the gradient phase does not emit it.
Its weights are ``dE/drho`` at fixed nuclei, the same object in both phases, so MOIST forms them once and leaves the contraction out of the gradient phase.
On a density-backed cavity the host therefore runs the response phase before the gradient and carries the weights over: contracted with its own ``d rho/dP`` they complete the Fock matrix, with the basis-center derivative ``d rho/dR`` the nuclear gradient.
Going straight to the gradient phase loses that term, and nothing reports the loss.

The host reads a response the way it walks a coupling. The response's ``next()`` makes each item current once per pass, in the order the model produced them, and returns false after the last one, rewinding for the next pass; on an empty response it returns false at once.
The item's name tells which item is current, and its arrays keep the names of the table and are its :doc:`fields`.
``get_response`` and ``get_gradient`` refill the response and restart the walk.
Contract every item the host knows and stop on any other: ``get_*`` reports an unanswered request, but nothing reports a response item the host skipped, and its term would silently be missing from the Fock matrix or the gradient.

.. _coupling-phases:

Lifecycle and answer reuse
--------------------------

After its first successful update, a model can create couplings.
A coupling belongs to that model; several may coexist, and each can serve every phase.
The language APIs describe ownership and lifetime requirements.

Any phase can be staged directly. ``get_*`` is the completeness check: it requires the matching staged phase and all its required outputs, fails by name on a wrong or absent staging and on missing outputs (those of the first incomplete request), and on failure leaves the accumulator unchanged.
Evaluation does not consume answers.

.. list-table:: State transitions
   :header-rows: 1

   * - Operation
     - Effect
   * - Create coupling after model update
     - Requests declared, nothing staged: ``next()`` returns false.
   * - ``prepare_energy``
     - Clear all answers, select energy, reset the cursor.
   * - ``prepare_response`` / ``prepare_gradient``
     - Retain all answers, select the phase, reset the cursor. Only a model
       update or ``prepare_energy`` discards answers.
   * - ``next()`` returns true
     - Move the cursor to the next request, in declaration order, with an
       output the staged phase requires and no valid answer for it.
   * - ``next()`` returns false
     - End the pass: the cursor rewinds and no request is current. The next
       call starts a new pass.
   * - Leave a loop early
     - Keep the cursor; the next ``next()`` resumes the pass.
   * - Answer an output
     - Update that output's missing state; the cursor stays.
   * - Update model
     - Invalidate every coupling's answers and staging and reset its cursor,
       even at unchanged size.

``next()`` visits each request at most once per pass and evaluates the missing state when it is called, so a second loop retries whatever is still missing, such as a rejected answer, and ends at once when nothing is.
On a coupling that is not staged, ``next()`` returns false without an error; ``get_*`` then reports the staging problem.

A request is current only after ``next()`` returned true.
It identifies the requested quantity, its named outputs and any request-specific inputs, such as the Gaussian moment ``width``; the inputs are its :doc:`fields`.
The language APIs describe whether the requirement state is read live or from a snapshot.

A missing output is a required output without a valid answer.
Asking whether an output the request does not have is missing gives false, not an error: ``answer`` rejects the unknown name and ``get_*`` names the missing outputs, so a misspelling cannot pass unnoticed.
A rejected answer (wrong name, rank, shape or a non-finite value) is reported on submission, leaves that output missing until a valid retry, and does not touch the other outputs.
An answer without a current request (before the first ``next()``, after ``next()`` returned false, after a new ``prepare_*`` or a model update) is an error that points to ``next()`` and changes no answer.

Prepare again after a model update.
For a new host density, start with ``prepare_energy`` to discard old answers; also update density-dependent cavities before preparing.
``next()`` never stops at a complete request, so replacing a valid answer, for example the potential of a perturbed density, also starts with ``prepare_energy``.

Energy and gradient accumulation
--------------------------------

``get_energy`` and ``get_gradient`` add their contribution to caller-owned accumulators in every language. Initialize them to zero for a fresh result, or to an existing contribution to add MOIST's term.
Repeated calls add repeated contributions; on failure the accumulator is unchanged.
Responses are replaced by the current phase's items; they are not accumulators.

The host contracts the potential adjoint, the density weights and the Gaussian amplitudes with its own integral and basis derivatives at fixed surface geometry.
MOIST accounts for cavity motion.

