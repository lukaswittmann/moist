# PCM migration: agent handoff

## Goal and source state

Port the PCM work from the **current live `../moist_dev/` working tree** into the current `moist_pcm` branch for a single PR targeting moist main. Do **not** commit `moist_dev` first; the uncommitted working tree is explicitly the source of truth. This is a port onto the post-PR-#18 channels/requests API, not a merge or cherry-pick.

All PCM models are in scope for this one PR: CPCM/COSMO, IEF-PCM, SS(V)PE, DPCM including original form, LPB, SLIC, and nonlocal.

The generated kernels originate from `../impl-code-gen/pcmgen/`, commit `b4467da`. There are Gaussian and point variants, with complex variants where generated. We are importing the generated implementation but are free to move files, rename files/modules, and update `use` statements. The hash/provenance checking machinery is being removed rather than preserved through renames. Keep a simple generator/provenance header unless there is a reason not to.

## Settled architecture

CPCM and COSMO move to the generated conductor kernel. The old `pcm/amat/` implementation is removed completely, including the public C symbols, Python helpers, and tests that directly depend on it. Cavity chain-rule tests formerly using `amat` should be rewritten through the conductor component; PySCF A-matrix comparisons become energy/charge-level comparisons. `moist_contract_pcm_nuclear_gradient` stays because it is a separate cavity-level API.

For models needing the normal derivative \(v_n\), use the existing raw-potential protocol rather than introducing a new request kind. Internally,
\[
v_n = \mathbf n\cdot\nabla\phi.
\]
Extend the potential request with `d2phi_dr2(3,3,ngrid)` and, for Gaussian requests, `d2phi_dr_dxi(3,ngrid)`. Extend the adjoint response with `w_dphi_dr(3,ngrid) = w_vn*n`. The host remains responsible only for raw electrostatic primitives.

The chosen interpretation is that \(v_n\) is the normal derivative of the **same smeared Gaussian potential** used for \(\phi\). This differs deliberately from the old dev behavior for Gaussian nuclear \(v_n\). The user selected protocol A; the smeared interpretation was then carried forward without objection. Treat that as the intended design, but it is the one decision that was not separately answered with an explicit “yes.”

Components may configure cavity requirements before `cavity%update`. Add a component `configure_cavity(cavity,error)` hook and a cavity interface along the lines of `cavity%configure(normal=.true., curvature=.true., error=error)`. Unsupported requested properties should fail early rather than silently no-op. DROP enables/computes requested normals/curvatures; iSwiG supplies stored normals and analytic curvature. Prefer accessors returning actual `k1,k2` rather than only mean curvature.

Use **one concrete type per model**, not dev's single layer type with integer variant dispatch. Naming follows the current repository convention: `solvation_model_component_*` types and `new_component_*` constructors. `phi_source`, `input_potential`, and the complex `epsilon_response` constructors are dropped. Complex generated states remain available internally for later response work.

Solver API gets a clean break. Public choice becomes essentially `{dense, iterative}`, where iterative means restarted GMRES with block-Jacobi preconditioning. Remove CG entirely. Retain the storage concept `{auto, dense, matrix_free}` from the dev solver infrastructure; the clean-cut proposal also adopted the dev-style defaults of approximately `tol=1e-8`, `maxiter=1000`, with automatic storage.

Expose one constructor per model in Fortran/C/Python. Physical model parameters are constructor arguments; numerical solver/Newton controls live in the PCM options/parameter type. LPB takes `kappa` directly; `debye_kappa(...)` remains a Fortran helper rather than entering the C ABI.

Retain `surface_weight_tangent` as a **Fortran-only, component-level hook**, FD-tested but not wired into the public API/general model flow yet. Generated tangent/response tests remain the kernel-level Hessian coverage. Full Hessian host derivatives are out of scope.

PySCF should support the new protocol end to end in this PR, including the second potential derivatives, `w_dphi_dr` Fock contraction, and nuclear gradients for DPCM/SLIC/nonlocal. The Fortran point-charge test host must likewise answer the new raw outputs.

## Tests and regression policy

Before deleting native `amat`, freeze CPCM/COSMO goldens from the current `moist_pcm` implementation on DROP and iSwiG. Also reuse the applicable dev parity/layer goldens. The dev parity data contain a cavity fingerprint, which helps distinguish cavity drift from PCM-kernel drift.

For models without suitable goldens, use FD tests for gradients/surface weights/tangents plus physics/property checks: Born-sphere limits broadly, Debye–Hückel for LPB, SLIC sphere formulas for both charge signs, LPB at \(\kappa=0\) versus IEF on a sphere, DPCM versus DPCM-original where they should coincide, and the nonlocal-to-local \(\lambda\to0\) limit. Add PySCF end-to-end checks and then freeze target-repo goldens as regression tripwires.

Accept the generated-kernel compile cost for now; do not add model-selection build switches or lower per-file optimization in this PR.

## Target layout

Use the established path/module-name correspondence. The intended structure is `src/moist/model/component/pcm/` containing `type.f90`, `solvers.f90`, `electrostatics.f90`, `local.f90`, `conductor.f90`, `layer.f90`, model-specific files (`cpcm.f90`, `cosmo.f90`, `iefpcm.f90`, `ssvpe.f90`, `dpcm.f90`, `lpb.f90`, `slic.f90`, `nonlocal.f90`), plus generated kernels split into:

`pcm/kernel/gauss/<model>.f90` → modules `moist_model_component_pcm_kernel_gauss_<model>`

`pcm/kernel/point/<model>.f90` → modules `moist_model_component_pcm_kernel_point_<model>`

Internal generated procedure/type names stay unchanged. Unit tests are flattened to the repository's `test/unit/test_<suite>.f90` convention. `pcm/amat/` disappears.

Docs should gain one PCM-model page per model, update the C reference for constructors/options/new potential outputs/new adjoint response/removal of `amat`, update PySCF docs, and add a changelog entry.

## Implementation order

1. Port math infrastructure: symmetric/complex LAPACK wrappers and Newton `custom_step`; adapt to the target repository's ILP64 wrappers rather than blindly copying dev.
2. Freeze native CPCM/COSMO goldens and import dev parity/layer golden data.
3. Add cavity/component configuration and curvature access.
4. Extend the host protocol with the second derivatives and `w_dphi_dr`, including C/Python interfaces and the Fortran test host.
5. Vendor/rename the generated Gaussian/point kernels and port `solvers.f90`, removing CG.
6. Rewrite PCM base/conductor, switch CPCM/COSMO to kernels, pass parity tests, then remove `amat` everywhere.
7. Add IEF-PCM, SS(V)PE, DPCM, LPB, SLIC, and nonlocal with FD/property/golden coverage.
8. Finish public C/Python constructors/options and full PySCF \(v_n\) support.
9. Update documentation and changelog.

The requested execution strategy is to parallelize independent early work with subagents/git worktrees. The planned first wave was: math/solver/kernel infrastructure, golden freezing, cavity hooks, host protocol, and PySCF integral work; base rewrite and concrete models come in later waves after their prerequisites land. The transcript shows the plan and agent prompts being prepared, but does **not** show evidence that those agents were actually launched.

## Current repository/build state

A baseline build was eventually obtained successfully. The working recipe required forcing the `jonquil`/`toml-f` subproject fallback because the conda `jonquil.pc` linkage was missing `toml-f`; this was identified as a pre-existing environment issue rather than a PCM regression. With that fallback, the baseline is **green: 69/69 tests pass**, including the PySCF suites. PySCF 2.13.1 is available in the environment.

The plan file has already been written/updated and the baseline recipe added to it. Locate that plan/scratchpad before starting duplicate planning work.

## Immediate takeover action

First inspect the current git/worktree state to see whether any subagents or branches were created after the transcript ended. If none exist, launch the initial independent worktrees from the settled plan. Preserve the green baseline and integrate in the dependency order above. The most important sequencing constraint is to capture native CPCM/COSMO goldens **before** deleting `amat/`.
