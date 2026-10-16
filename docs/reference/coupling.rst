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

Construct the cavity and model components, update the model for the structure, and create a coupling. Reuse it throughout the SCF calculation and gradient evaluation.

At each SCF density:

#. **Build the cavity.** Follow the :ref:`coupling-build` with the current density.
   A cavity determined only by nuclear geometry can be reused while the geometry is unchanged.
#. **Evaluate the energy.** Prepare the energy phase, answer its requests and obtain MOIST's energy contribution.
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

   update model (structure)
   create coupling
   prepare build
   while a request with missing outputs is available:
       answer it (density matrix, nuclear charges)
   build model

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
Read them from the model; the coupling keeps no grid copy. Grid outputs require a built surface.

.. _coupling-build:

Build phase
-----------

``update(structure)`` prepares the geometry. Cavities that need host inputs defer their surface construction to the build phase.
Run it after each structure update and each density change on a density-dependent cavity:

#. Call ``prepare_build(coupling)``.
#. Answer the requests with the usual ``next()`` / ``answer`` loop.
#. Call ``build(coupling)`` to complete the cavity and components.

Internal isodensity cavities request ``density_matrix``; both isodensity implementations request ``nuclear_charges`` when core fill is enabled.
Callback density is evaluated directly, outside the coupling.
Geometric cavities and callbacks without core fill are built by ``update``: their request loop is empty and ``build`` does not rebuild the surface.
For those callbacks, repeat ``update`` after each density change.

Standalone cavities use the same protocol in every language.
See the :ref:`Fortran <isodensity-build-fortran>`, :ref:`C <isodensity-build-c>` and :ref:`Python <coupling-python>` examples, and :ref:`coupling-phases` for answer reuse and invalidation.

.. _coupling-requests:

Requests and outputs
--------------------

Requests, outputs and response items are addressed by their name in every language.
The table uses Fortran shapes; C and Python reverse the dimensions, so ``dphi_dr`` and ``d2phi_dr_dxi`` are ``(ngrid, 3)``, while ``d2phi_dr2`` and ``mt`` are ``(ngrid, 3, 3)``.

.. list-table:: Requests, their grid inputs and their outputs
   :header-rows: 1

   * - Request
     - Cavity inputs
     - Host outputs
   * - ``point_potential``
     - ``xyz``
     - ``phi(ngrid)``, ``dphi_dr(3, ngrid)``, ``d2phi_dr2(3, 3, ngrid)``
   * - ``gaussian_potential``
     - ``xyz``, ``xi0``
     - ``phi``, ``dphi_dr``, ``dphi_dxi(ngrid)``, ``d2phi_dr2``, ``d2phi_dr_dxi(3, ngrid)``
   * - ``gaussian_moments``
     - ``xyz`` and the request's own ``width(ngrid)``
     - ``gt(ngrid)``, ``pt(3, ngrid)``, ``mt(3, 3, ngrid)``, ``rt(3, ngrid)``
   * - ``nuclear_charges``
     - none
     - ``z(natom)``
   * - ``density_matrix``
     - the request's own ``nbasis``; the basis order of ``moist_get_isodensity_cart_layout``
     - ``dcart(nbasis, nbasis)``, symmetric

Both potential requests carry the total nuclear plus electronic potential.
``dphi_dr`` is its derivative with respect to the evaluation point, without a charge multiplication or electric-field sign change. ``d2phi_dr2`` is its spatial Hessian; ``d2phi_dr_dxi`` differentiates that same spatial derivative with respect to the Gaussian inverse width ``xi``.

PCM uses the Gaussian operator ``erf(xi*r)/r`` consistently in the potential, Fock matrix and nuclear gradient, including the normal derivative :math:`v_n=\mathbf n\cdot\nabla\phi` needed by some variants. ``dphi_dxi`` differentiates the inverse length ``xi`` in bohr**-1, not the exponent ``xi**2``.
The point operator is a separate request for components requiring bare Coulomb potentials.

Gaussian moment widths are exponents in bohr**-2, chosen by GOSTSHYP as ``pi*ln(2)/a`` from the per-point area.
Every moment uses the unit-integral Gaussian ``G=(width/pi)**1.5 * exp(-width*|r-C|**2)``.
The same prefactor applies to ``gt``, ``pt``, ``mt`` and ``rt``, and to the host's AO response and nuclear integral derivatives. Ordinary orbital normalization is different.
Zero-width requests represent inert zero-area points and have zero moments.
They are a :doc:`field <fields>` of the request, not of the cavity: read them from the request, never recompute them.

``nuclear_charges`` are atomic numbers minus ECP core electrons; all-electron hosts supply the atomic numbers.
``density_matrix`` uses the internal isodensity level set's Cartesian-monomial layout; its size ``nbasis`` is an input of the request, like the moment widths, not a count of the snapshot.
Neither has a response item.

PCM requires ``phi`` in every phase.
It requires both derivatives for the gradient, and for the response only when the cavity geometry follows the density.
The gradient also requires ``z``, reused from the build phase when available.
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
     - ``w_phi(ngrid)``, ``w_dphi_dr(3, ngrid)``
     - Any electrostatic component is present
   * - ``density``
     - ``w_rho(ngrid)``, ``w_grad_rho(3, ngrid)``, ``w_hess_rho(3, 3, ngrid)``
     - The cavity surface follows the host density (the isodensity level sets),
       response phase only
   * - ``gaussian_amplitude``
     - ``w_overlap(ngrid)``, ``w_normal_deriv(ngrid)``
     - A GOSTSHYP component is present

The potential adjoint ``w_phi = dE/dphi`` is the induced surface charge of a stationary PCM (CPCM, COSMO); a non-symmetric response matrix (IEF-PCM, SS(V)PE) emits the symmetrized adjoint instead, which the host contracts identically. Variants that use :math:`v_n` also return ``w_dphi_dr``. Otherwise this field is absent from the native field list and is ``None`` in Python. The host contracts this vector with its derivative of the raw spatial potential gradient, including the electronic contribution to the Fock matrix and the nuclear contribution to the gradient.
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

After building, any evaluation phase can be staged directly. ``get_*`` requires the matching phase and all required outputs; errors name the missing phase or outputs and leave the accumulator unchanged.
Evaluation does not consume answers.

.. list-table:: State transitions
   :header-rows: 1

   * - Operation
     - Effect
   * - Create coupling
     - Requests declared, nothing staged: ``next()`` returns false.
   * - ``prepare_build``
     - Clear answers, select build, reset the cursor. A deferred surface is
       discarded; evaluation requires a successful build.
   * - ``build``
     - Complete the model, retain build answers, end staging and invalidate
       other couplings. A failed deferred build leaves no surface.
   * - ``prepare_energy``
     - Clear answers except build inputs, select energy, reset the cursor.
   * - ``prepare_response`` / ``prepare_gradient``
     - Retain answers, select the phase, reset the cursor.
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
Native rejection of an output's rank, shape or values clears that output without touching others; an unknown name changes nothing.
Python validates shapes and density-matrix symmetry before submission; those errors leave any previous answer intact (:ref:`coupling-python`).
An answer without a current request (before the first ``next()``, after ``next()`` returned false, after a new ``prepare_*`` or a model update) is an error that points to ``next()`` and changes no answer.

Prepare again after a model update.
For a new host density, follow the SCF workflow above.
``next()`` never stops at a complete request, so replacing a valid answer, for example the potential of a perturbed density, also starts with ``prepare_energy``.

Energy and gradient accumulation
--------------------------------

``get_energy`` and ``get_gradient`` add their contribution to caller-owned accumulators in every language. Initialize them to zero for a fresh result, or to an existing contribution to add MOIST's term.
Repeated calls add repeated contributions; on failure the accumulator is unchanged.
Responses are replaced by the current phase's items; they are not accumulators.

The host contracts the potential adjoint, the density weights and the Gaussian amplitudes with its own integral and basis derivatives at fixed surface geometry.
MOIST accounts for cavity motion.
