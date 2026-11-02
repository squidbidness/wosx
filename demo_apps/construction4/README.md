# Construction 4: experimental CPU 2D implementation

This independent solver implements the practical tangent-disk recurrence in
`arh_phd/construction4_codex_mathematical_handoff.md`. It does not change WoSt.
It drops the cap error C* and normal-displacement error D0. Finite recursive
surface weight is one; unbiasedness, global convergence, almost-sure termination,
and finite variance are not established. Monte Carlo standard error describes
sampling error, not discretization bias.

The library uses WoSX's convention **-Delta u=s**. In the handoff, Delta u=f,
so s=-f. With constant source the implemented updates are:

- Interior disk radius rho: add s rho^2/4, then sample its circle uniformly.
- Near Neumann boundary: add (gN+s r/2)/M_R*, then sample the complementary arc.
- Near Dirichlet boundary: project and return gD plus accumulated contributions.

The tangent disk is tangent at the current interior point y, not its physical
boundary projection y0. Its center is y-r n(y0). The solver requires one combined whole-boundary nearest-point callback returning
distance, closest point, outward normal, boundary type, normal uniqueness, and a
reliable inside/outside sign. Each visited walk position uses exactly one such
query, including the final termination check. It reuses the returned point and
normal without re-projecting. Boundary type is chosen by the nearest point; it
does not search a farther Dirichlet boundary before a nearer Neumann boundary.
Each proposed disk center uses one additional whole-boundary query to certify
containment/clearance, including rejected radius trials. The analytic rectangle
adapter needs no ray, silhouette, or star-radius operations for Construction 4.
A separate engine-mesh benchmark is described below; mesh support is not enabled
in the original image-based `basic_2d_demo`. No speculative Robin or generalized operator kernel is used.

## Build and Run

From the repository root (the graphics options avoid unnecessary desktop dependencies):

```bash
cmake -S . -B build -DWOSX_ENABLE_GPU_SUPPORT=OFF -DWOSX_BUILD_DEMO_APPS=ON \
  -DWOSX_BUILD_BINDINGS=OFF -DPOLYSCOPE_BACKEND_OPENGL3_GLFW=OFF \
  -DPOLYSCOPE_BACKEND_OPENGL3_EGL=OFF -DTCM_BUILD=OFF
cmake --build build --target construction4_benchmark -j4
./build/demo_apps/construction4_benchmark --self-test
./build/demo_apps/construction4_benchmark --solver both --case polynomial --walks 6000 --seed 20261030
```

Domain: (-1,1)x(0,1), bottom Neumann, other sides Dirichlet. Polynomial:
u=1+0.3z+0.2x^2+0.1z^2, s=-0.6, gN=-0.3. At (0,0.1), u=1.031.
Also available: `--case harmonic` with u=1+0.3z+0.2(x^2-z^2), s=0;
`--case constant` with u=1, s=gN=0.

The benchmark runs the unchanged WoSt solver on the same analytic rectangle
geometry/PDE. WoSt retains its original separate boundary queries. It does not
benchmark FCPW BVH traversal performance. CSV counts distinguish:

- `point_nearest_queries`: one whole-boundary query per Construction 4 visited
  position (steps + 1 per completed walk, including final termination).
- `center_nearest_queries`: one whole-boundary query per proposed center,
  including radius rejections.
- `total_nearest_queries`: both of the above for Construction 4; separate
  distance and projection calls for unchanged WoSt.
- WoSt ray, boundary-sampling, and star-radius calls: separate operation types.

The Construction 4 adapter searches all four rectangle edges in one operation;
it does not call the separate Dirichlet/Neumann query callbacks. The initial
point's combined query already supplies membership information, so no additional
initial certificate is needed. Query counters are atomic for parallel evaluation.
Initial sample metadata is computed analytically outside the measured walk.
Neumann encounters count tangent steps for Construction 4 and reflecting-boundary
states for WoSt. These are algorithm-specific events.

Baseline defaults: d=0.01, r0=0.4, radius reduction 1/2, epsilonD=epsilonN=0.01.
Use `--cap`, `--radius`, `--tol-dirichlet`, `--tol-neumann`, `--max-steps`,
`--walks`, and `--seed` to change them. Cap fallback uses d_eff=min(d,1.9r)
only when d>=2r. Radius rejections, cap adjustments, min/max effective cap,
termination failures, signed errors and sampling SE are reported. Keep r0 fixed
when sweeping cap size. Repeat across seeds and compare matched error/cost;
a single low signed error can reflect bias cancellation and sampling noise.

Step-limit and geometry failures abort a batch without committing its sample
statistics, rather than reporting zero values or survivor-conditioned means.
Completed walks before failure are counted in diagnostics. RNG state and attempt
counters are not rolled back. The benchmark refuses comparison means if either
solver fails. Fixed seeds are reproducible for a fixed executable/toolchain;
WoSt and Construction 4 consume different random trajectories.

## Validation and Scope

`--self-test` checks positive-eta Green/cap identities including E0 and EC,
independent finite-volume quadrature, complementary mass and scalar volume limits,
inverse-CDF distribution and chord cutoff, a Stage B analytic continuation mean,
constant-solution recursion, mode guards, failure handling, near-junction cap
fallback, polynomial/harmonic manufactured recursion against unchanged WoSt, and
seeded serial/parallel reproducibility including query counts. Manufactured
regressions allow six sampling standard errors plus an explicit absolute bias
allowance of 0.01; passing is not evidence that this bias vanishes under refinement. It does not prove the general
singular limit interchanges. Finite-source tests retain the handoff's Delta u=f
notation explicitly; production solver callbacks always use s=-f.

Only constant-source, zero-absorption, single-sided, watertight mixed interior
problems and solution estimates are currently accepted. Set PDE::isSourceConstant
truthfully, including for zero source. Gradient, Robin, screened Poisson, variable
source, pure Neumann, double-sided/exterior, roulette, splitting, and contribution
suppression modes are unsupported. A geometry callback must provide unique outward
Neumann normals and a conservative signed clearance in the combined nearest-point result. Runtime guards cannot prove
that a user-provided callback is correct.


## Engine Mesh Comparison

The separate `construction4_mesh_benchmark` target uses the supplied engine mesh
and centroid-based Dirichlet/Neumann labels without editing the supplied loader,
data, or WoSt solver. Six zero-length closing segments are omitted from the
Construction 4 BVH and coincident endpoints welded. The nonzero boundary is
unchanged. One full-boundary FCPW BVH supplies nearest point, segment label and
signed distance. Segment normals are used for tangent steps; vertex pseudonormals
are used only for signing distances.

Before any comparisons, the adapter checks closed oriented topology,
non-adjacent crossings/contacts, every edge's inward/outward orientation, nearest
query signs against independent winding numbers (including near vertices), and
BVH distances against brute-force segment distances. Numerical validation is not
a proof of robustness for arbitrary geometry.

Build and run from the repository root:

```bash
cmake -S . -B build
cmake --build build --target construction4_mesh_benchmark -j4
./build/demo_apps/construction4_mesh_benchmark --validate-only
./build/demo_apps/construction4_mesh_benchmark --case manufactured
./build/demo_apps/construction4_mesh_benchmark --case engine
```

Use `--root` if launched elsewhere. `--config` is relative to root. Other options:
`--points`, `--walks`, `--seed`, `--cap`, `--radius`, `--tol-dirichlet`,
`--tol-neumann`, `--max-steps`, `--reference-walks`, `--reference-tolerance`.
Defaults: 16 fixed interior points, 6000 walks per point, d=epsilonD=epsilonN=0.01,
r0=0.4. Identical points are selected from a fixed 64x64 lattice outside a
boundary strip of width 0.02 times the mesh diagonal. Results concern those
selected interior points, not a full image or uniform error across the domain.

`manufactured` uses the engine geometry and partition with
u=1+0.3z+0.2x^2+0.1z^2, s=-0.6, exact Dirichlet values and Neumann data grad(u).n.
Errors are measured against the exact polynomial. `engine` uses the supplied
textured Dirichlet data with zero source and homogeneous Neumann data: the
supplied configuration disables both contributions despite containing nonzero
textures. WoSt keeps those suppression flags in engine mode; Construction 4 uses
equivalent zero callbacks. Variable forcing is not enabled.

In engine mode, the reference is independently seeded WoSt with 65536 walks per
point and tolerance 1e-4. Differences are reference discrepancies, not exact
errors. Reference SE and combined sampling SE are reported. The reference has
unknown residual bias. RMSE across points and RMS sampling SE are distinct.

Mesh mode explicitly enables ordinary-WoS fallback when a nearest Neumann point
has a nonunique corner normal, using the already obtained distance. This adds no
query and is separately counted. It does not switch to WoSt. The mesh center
acceptance check also requires an extra clearance margin of
8*float_epsilon*mesh_diagonal. Degenerate steps or failed radius searches abort
rather than reporting survivor-conditioned means. The corner policy and
approximate tangent rule have no global convergence guarantee.

The lumped WoSt query total includes distance, projection, star-radius, ray,
boundary-sampling AND visibility calls (the older rectangle tables did not
instrument visibility). Construction 4 queries include visited-point and center
nearest searches. Center checks, cap fallbacks, corner fallbacks and radius
rejections are also reported per walk. API calls have different costs; wall times
must accompany counts. Single-threaded solve timings exclude initial sample
metadata, mesh loading, BVH builds, geometry validation, evaluation-point
selection and reference computation for both methods. Reference times are
reported separately in per-point rows. Repeat seeds and refine tolerances before
claiming matched-accuracy performance.
