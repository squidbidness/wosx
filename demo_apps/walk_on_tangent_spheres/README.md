# Walk on Tangent Spheres: experimental CPU 2D implementation

This independent solver implements the practical tangent-disk recurrence in
`arh_phd/walk_on_tangent_spheres_codex_mathematical_handoff.md`. It does not change WoSt.
It drops the cap error C* and normal-displacement error D0. Finite recursive
surface weight is one; unbiasedness, global convergence, almost-sure termination,
and finite variance are not established. Monte Carlo standard error describes
sampling error, not discretization bias.

The library uses WoSX's convention **-Delta u=s**. In the handoff, Delta u=f,
so s=-f. With constant source the implemented updates are:

- Interior disk radius rho: add s rho^2/4, then sample its circle uniformly.
- Near Neumann boundary: add (gN+s r/2)/M_R*, then sample the complementary arc.
- Near Dirichlet boundary: project and return gD plus accumulated contributions.

The tangent disk is tangent at the current interior point y. WoTS uses the
existing `GeometricQueries` interface and the same loaded representation as WoSt.
`projectToBoundaryWithType` queries the existing absorbing and reflecting
structures and retains the closer result's partition label. Equal distances
prefer absorbing/Dirichlet, as in the original projection callback. The legacy
`projectToBoundary` delegates to the typed callback and discards the label.

The approach direction n=(y0-y)/|y0-y| supplies the normal used for both the
center c=y-r*n and the Neumann-data callback. The artificial source y'=y-eta*n
lies on the same line, so n also points from y' toward y0. The primitive normal
returned by projection is ignored, including at corners. This is a documented
corner convention; its directional-derivative discrepancy remains part of the
omitted error, not a theorem about face-specific Neumann traces. An exactly
coincident Neumann point y=y0 has no defined approach direction and fails
explicitly. There is no face-normal selection or corner-WoS fallback.

For every proposed center, existing `computeDistToBoundary` gives unsigned
clearance and `insideDomain` determines membership. Their defaults use the two
partitioned structures and ray-intersection parity, respectively. No extra BVH,
mesh welding, segment removal, adjacency, pseudonormals, or mesh-validation
pipeline is required or performed. The optional extra `containmentMargin`
defaults to zero; the old adapter's fixed float-BVH margin is no longer imposed.
Callback arithmetic still has finite precision, and this refactor does not
establish a general numerical containment or convergence guarantee.

The solver checks initial membership once per point's batch. Each subsequent
visited position uses typed projection; the final check is included. The
analytic rectangle implements membership analytically, while the supplied mesh
uses the original WoSX membership callback. A separate engine benchmark is
available; the original `basic_2d_demo` solver selector remains unchanged.

## Build and Run

From the repository root (the graphics options avoid unnecessary desktop dependencies):

```bash
cmake -S . -B build -DWOSX_ENABLE_GPU_SUPPORT=OFF -DWOSX_BUILD_DEMO_APPS=ON \
  -DWOSX_BUILD_BINDINGS=OFF -DPOLYSCOPE_BACKEND_OPENGL3_GLFW=OFF \
  -DPOLYSCOPE_BACKEND_OPENGL3_EGL=OFF -DTCM_BUILD=OFF
cmake --build build --target walk_on_tangent_spheres_benchmark -j4
./build/demo_apps/walk_on_tangent_spheres_benchmark --self-test
./build/demo_apps/walk_on_tangent_spheres_benchmark --solver both --case polynomial --walks 6000 --seed 20261030
```

Domain: (-1,1)x(0,1), bottom Neumann, other sides Dirichlet. Polynomial:
u=1+0.3z+0.2x^2+0.1z^2, s=-0.6, gN=-0.3. At (0,0.1), u=1.031.
Also available: `--case harmonic` with u=1+0.3z+0.2(x^2-z^2), s=0;
`--case constant` with u=1, s=gN=0.

The benchmark runs unchanged WoSt on the same analytic rectangle geometry/PDE.
It does not benchmark FCPW BVH traversal. CSV counts distinguish logical checks
from actual partition-query calls:

- `point_nearest_queries`: typed projection invocations, steps+1 per completed
  walk. By default each invokes one absorbing and one reflecting projection.
- `center_nearest_queries`: proposed-center clearance checks, including rejected
  radius trials. By default each invokes both boundary-distance callbacks.
- `total_nearest_queries`: actual absorbing/reflecting distance and projection
  calls; compound callbacks are not counted again.
- `membership_checks`: initial and proposed-center membership invocations. The
  rectangle computes these analytically; mesh membership performs ray searches.
- WoSt ray, boundary-sampling, and star-radius calls: separate operation types.

Counters are atomic for parallel point evaluation. Initial sample metadata is
computed outside measured walks. Neumann encounters count tangent steps for WoTS
and reflecting-boundary states for WoSt, which are algorithm-specific events.

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
WoSt and Walk on Tangent Spheres consume different random trajectories.

## Validation and Scope

`--self-test` checks positive-eta Green/cap identities including E0 and EC,
independent finite-volume quadrature, complementary mass and scalar volume limits,
inverse-CDF distribution and chord cutoff, a Stage B analytic continuation mean,
constant-solution recursion, mode guards, failure handling, near-junction cap
fallback, polynomial/harmonic manufactured recursion against unchanged WoSt, and
seeded serial/parallel reproducibility including query counts. Manufactured
regressions allow six sampling standard errors plus an explicit absolute bias
allowance of 0.01; passing is not evidence that this bias vanishes under
refinement. It does not prove the general singular limit interchanges. Finite-source tests retain the handoff's Delta u=f
notation explicitly; production solver callbacks always use s=-f.

Only constant-source, zero-absorption, single-sided, watertight mixed interior
problems and solution estimates are currently accepted. Set PDE::isSourceConstant
truthfully, including for zero source. Gradient, Robin, screened Poisson, variable
source, pure Neumann, double-sided/exterior, roulette, splitting, and contribution
suppression modes are unsupported. Geometry providers must supply the existing
projection, distance and membership callbacks. Runtime guards cannot prove that
user-provided callbacks are correct.

## Engine Mesh Comparison

`walk_on_tangent_spheres_mesh_benchmark` loads the supplied `ModelProblem<float>`
and uses its existing Dirichlet/Neumann meshes and acceleration structures for
both solvers. It adds no geometry representation or preprocessing. The old mesh
adapter and its topology, winding, normal and brute-force distance checks have
been removed entirely; `--validate-only` is no longer an option.

```bash
cmake -S . -B build
cmake --build build --target walk_on_tangent_spheres_mesh_benchmark -j4
./build/demo_apps/walk_on_tangent_spheres_mesh_benchmark --case manufactured
./build/demo_apps/walk_on_tangent_spheres_mesh_benchmark --case engine
```

Run from the repository root or set `--root`; `--config` is relative to root.
Options include `--points`, `--walks`, `--seed`, `--cap`, `--radius`,
`--tol-dirichlet`, `--tol-neumann`, `--max-steps`, `--reference-walks`, and
`--reference-tolerance`. Defaults: 16 fixed interior points, 6000 walks per point,
d=epsilonD=epsilonN=0.01 and r0=0.4. Both methods use identical evaluation points
from a 64x64 lattice, selected with the existing membership/distance queries
outside a boundary strip of width 0.02 times the model bounding-box diagonal.
Results concern selected interior points, not full-domain uniform error.

`manufactured` uses u=1+0.3z+0.2x^2+0.1z^2, source s=-0.6, exact Dirichlet values,
and the callback grad(u).n. WoSt supplies its sampled face normal; WoTS supplies
its approach direction, giving an explicit radial extension of the manufactured
normal derivative at corners. Exact errors and Monte Carlo SE are reported.

`engine` uses supplied textured Dirichlet values with zero source and homogeneous
Neumann data, matching the supplied configuration's disabled source/reflecting
contributions. WoSt retains those suppression flags; WoTS uses equivalent zero
callbacks. Variable source is not enabled. An independently seeded WoSt reference
uses 65536 walks per point at tolerance 1e-4. Differences are reference
discrepancies, not exact errors. Reference and combined sampling SE are reported;
reference discretization bias remains unknown. RMSE and RMS sampling SE differ.

`query_statistics.h` only instruments the existing callbacks, in place. It is
not a geometry adapter. Both lumped query totals include underlying distance,
projection, star-radius, ray-intersection, boundary-sampling, visibility and
all-hit ray searches. Membership invocations are reported separately; their
underlying all-hit searches are already included in totals and are not counted
twice. Default 2D membership performs two ray directions and searches each
boundary partition. Center checks, cap fallbacks and radius rejections are also
reported per walk. The former corner-fallback statistic is removed.

Solve timings are single-threaded and exclude model loading, shared BVH builds,
point selection, initial point metadata and reference computation for both
methods. Reference times are reported separately. All input geometry is passed
through the same original loader. Unsupported modes, zero steps and failed
radius searches abort without reporting failed-batch survivor means. Repeated
seeds and tolerance refinement remain necessary for accuracy/performance claims.
