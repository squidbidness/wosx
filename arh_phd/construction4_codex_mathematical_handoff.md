# Construction 4: Mathematical and computational handoff for Codex

**Status:** research specification, not a proof of a convergent PDE algorithm.  
**Prepared:** 2026-10-30.  
**Authoritative source:** `construction_4_checkpoint_revised_2026-10-26.md` (saved checkpoint).  
**Related prototypes:** `stage_a_construction4.py`, `stage_a_limit_interchanges.py`, `stage_b_single_step.py`, `stage_c_recursive_wos.py`, `stage_c_radius_heuristic.py` and later sweep files.  
**Scope:** two-dimensional mixed Dirichlet–Neumann Poisson problems first; extension to three dimensions, Robin data, and general forcing is future work.

> **Instructions to Codex:** Preserve the mathematical distinctions between exact equations, conditional limiting equations, and approximations. Do **not** implement an unproved limiting interchange as though it were a theorem. Do **not** conflate the tangent-ball method below with the subsequently explored *generalized boundary-operator kernel*; the latter is a separate, speculative line of work. Never silently drop error terms in descriptions of an “exact” estimator. Keep the existing WoSt implementation unchanged during initial integration.

## 1. Motivation and intended algorithm

Ordinary Walk on Spheres (WoS) solves a Dirichlet problem by sampling uniform points on admissible spheres and terminating near Dirichlet boundary. Mixed Dirichlet–Neumann conditions require a rule for continuing rather than terminating at Neumann boundary encounters. Construction 4 seeks to use **ordinary spheres tangent at the current near-boundary walk point**, with a Green identity whose source point approaches that tangent point. A limiting flux is highly concentrated near the tangent point; instead of sampling extremely short steps, partition the sphere into a short-chord cap `C` and its complement `R`, approximate the unknown solution as locally constant on the cap, then solve algebraically for the solution at the original walk point. The complementary flux, properly normalized, defines **one recursive jump of weight exactly 1**, whose chord distance is at least `d`.

Potential advantages to investigate: nonbranching recursion; no special ray-intersection or star-radius geometry queries for the new solver; explicit cap error; a nearest-boundary/clearance-query-only geometry interface. These are motivations and hypotheses, **not** proven comparative advantages over Walk on Stars (WoSt).

## 2. Notation, geometry, and conventions

* `Ω ⊂ ℝ^n` is the physical domain, initially `n=2`; `Γ_D` and `Γ_N` are its Dirichlet and Neumann boundary parts.
* Solve `Δu=f` in `Ω`, `u=g_D` on `Γ_D`, and `∂_n u=g_N` on `Γ_N`, where `n` is the **outward** normal of `Ω`.
* Classical WoS reaches `y∈Ω` close to `Γ_N`; `y_0` is a nearest point on `Γ_N`, and `δ=|y−y_0|≤ε_N`. In a smooth nearest-point tubular neighborhood, `y = y_0−δ n(y_0)`.
* Choose an admissible **closed** ball `B=B_r(c) ⊂ \overline Ω` whose sphere `S=∂B` is tangent at **the current point** `y` (not generally at `y_0`), with `c=y−r n(y_0)` (center points inward). Assume its interior is inside `Ω` and its surface does not cross any physical boundary. Smoothness/normal uniqueness is required for this construction.
* Introduce the artificial source `y'=y−η n(y_0) ∈ B`, `0<η<2r`; take `η→0+` **with `r` and the positive cap chord radius `d` fixed**. Normal derivatives `∂_{n_x}` on `S` use the outward normal of the *ball* (distinct from the physical normal except at `y`).
* `C={x∈S:|x−y|≤d}`, `R=S\C`. Require `0<d<2r` so `R` has positive measure. Here `R` names the complementary part of the surface, **not** the ball radius, which is `r`.
* `G_η(x)=G_B(x,y')` is the zero-Dirichlet Green function for the ball in the **sign convention** `Δ_x G_η=δ_{y'}`, `G_η|_S=0`. Consequently `G_η≤0` in the interior and its outward-normal flux `P_η(x)=∂_{n_x}G_η(x)` is nonnegative on `S` (the Poisson density). `Δu=f`.
* `M_C(η)=∫_C P_η dS`, `M_R(η)=∫_R P_η dS`; `M_C+M_R=1` because `∫_S ∂_{n_x}G_η=1`.
* `E_C(η)=∫_C [u(x)−u(y)]P_η(x)dS` and `E_0(η)=u(y)−u(y')−η g_N(y_0)`. **This definition of `E_C` uses `u(y)`, not `u(y')`.**

## 3. Exact finite-η theorem: algebraic identity (ESTABLISHED)

**Theorem 1 (finite-source cap identity).** Suppose `B` lies in `Ω`, `u∈C²(B)∩C(\overline B)` (or sufficient Green-representation regularity), `y'∈B`, and `0<d<2r`. Then

```
u(y') = ∫_S u(x) P_η(x) dS_x + ∫_B G_η(x) f(x) dx.
```

Splitting `S=C∪R`, adding and subtracting `u(y)` on `C`, and applying the definitions gives

```
u(y') = M_C u(y) + ∫_R u(x)P_η(x)dS_x
        + ∫_B G_η(x) f(x)dx + E_C(η).
```

Since `u(y)=u(y')+η g_N(y_0)+E_0(η)` and `M_R=1−M_C`, exact rearrangement yields

\[
\boxed{
 u(y') = \frac{\eta}{M_R(\eta)}\left[
 M_C(\eta)\left(g_N(y_0)+\frac{E_0(\eta)}{\eta}\right)
 +\int_R u(x)\frac{P_\eta(x)}{\eta}\,dS_x
 +\int_B\frac{G_\eta(x)}{\eta}f(x)\,dx
 +\frac{E_C(\eta)}{\eta}
 \right].}
\]

**Proof:** Green's representation plus the exact cap split and one algebraic substitution above. This is an exact identity for all positive `η` where `M_R>0`. No approximation and no limit interchange have occurred. The known boundary datum `g_N(y_0)` can be evaluated at the nearest physical boundary point; `E_0` captures its discrepancy from the derivative at `y`.

## 4. Limiting results already justified (ESTABLISHED under stated hypotheses)

Define the pointwise normalized limits away from `y`:

\[
G^*(x):=\lim_{\eta\downarrow0}G_\eta(x)/\eta,\quad
P^*(x):=\lim_{\eta\downarrow0}P_\eta(x)/\eta.
\]

**Theorem 2 (complementary flux convergence).** For a fixed smooth ball, fixed cap `d>0`, and bounded `u` on `R`, the normalized Poisson flux converges uniformly on `R` to `P*`. Hence

\[
M_R^*:=\int_R P^*\,dS>0,\quad
\lim_{\eta\downarrow0}\frac{M_R(\eta)}\eta=M_R^*,\quad
\lim_{\eta\downarrow0}\frac\eta{M_R(\eta)}=\frac1{M_R^*},
\]

\[
\lim_{\eta\downarrow0}\int_R u(x)\frac{P_\eta(x)}\eta\,dS_x
=\int_R u(x)P^*(x)\,dS_x.
\]

**Proof:** The ball Poisson kernel is smooth in source and field coordinates on the compact set separated from `y` by chord distance `d`; Taylor expansion in source depth `η` gives uniform convergence of `P_η/η`. Uniform convergence permits termwise integration; positivity gives `M_R*>0`. This excludes the singular cap and does not justify *any* cap or volume interchange.

**Theorem 3 (normal-displacement error).** If `u` is differentiable along `n(y_0)` at `y`, then

\[
\lim_{\eta\downarrow0}\frac{E_0(\eta)}{\eta}
= \nabla u(y)\cdot n(y_0)-g_N(y_0).
\]

**Proof:** `y'=y−η n`, so `u(y')=u(y)−η∇u(y)·n+o(η)`.

**Theorem 4 (unit-weight one-step sampling).** On `R` set `p_R(x)=P*(x)/M_R*`. Then `p_R` is a probability density and

\[
\frac1{M_R^*}\int_R u(x)P^*(x)dS_x=\mathbb E[u(X)],\quad X\sim p_R.
\]

**Proof:** By Theorem 2, `P*≥0` and its mass on `R` is `M_R*`; the equality is the definition of expectation. This proves unit weight for the **recursive surface term alone**, not finite variance or almost-sure termination of the full algorithm.

## 5. The checkpoint representation (CONDITIONAL, not fully proved)

Assuming that the outstanding scalar limits exist and that the finite-η equality passes to the limit, the saved checkpoint gave the representation

\[
\boxed{
\begin{aligned}
u(y)={}&\mathbb E[u(X)]+\frac{g_N(y_0)}{M_R^*}
+\frac{\nabla u(y)\cdot n(y_0)-g_N(y_0)}{M_R^*}\\
&+\frac1{M_R^*}\lim_{\eta\downarrow0}\int_B\frac{G_\eta(x)}\eta f(x)\,dx\\
&+\frac1{M_R^*}\lim_{\eta\downarrow0}\int_C[u(x)-u(y)]\frac{P_\eta(x)}\eta\,dS_x.
\end{aligned}}
\]

**Important:** For this equation the volume and cap-error limits deliberately remain **outside** their integrals. It is not yet a generally proved formula for an arbitrary smooth mixed-boundary domain. The two first numerator terms involving `g_N` and derivative discrepancy sum to `∇u(y)·n/M_R*`, but the split is kept because it distinguishes known Neumann data from approximation error.

For programming convenience define the two *scalar limits*, if they exist,

\[
V^*:=\lim_{\eta\downarrow0}\int_B (G_\eta/\eta)f\,dx,\qquad
C^*:=\lim_{\eta\downarrow0} E_C(\eta)/\eta.
\]

Then write exactly the conditional representation as `u(y)=E[u(X)]+[g_N(y_0)+V*+D_0+C*]/M_R*`, with `D_0=∇u(y)·n−g_N(y_0)`.

## 6. Candidate kernel-interchanged representation (UNPROVED IN GENERAL)

The tempting further identities are

\[
V^*\stackrel{?}=\int_B G^*(x) f(x)\,dx,
\qquad
C^*\stackrel{?}=\int_C [u(x)-u(y)]P^*(x)\,dS_x.
\]

**Do not claim these as established from the saved checkpoint.** The source approaches the boundary of the ball, so `G_η/η` develops a near-boundary singularity inside `B`; a suitable dominated-convergence or uniform-integrability argument is needed. On the cap, `P*` is singular at `y`; tangential odd terms may need symmetric cancellation, a principal-value formulation, or subtraction of a first-order Taylor polynomial. A conditional proof for special regularity classes or an explicit ball calculation is not a general theorem. In particular, a raw absolute integral of `[u(x)−u(y)]P*` need not exist for a generic smooth `u` in two dimensions if the tangential derivative at `y` is nonzero.

**Numerical evidence:** `stage_a_limit_interchanges.py` independently compared finite-η normalized integrals over `B`, `R` and `C` with independently computed proposed kernel-limit integrals for the symmetric polynomial test in §8. Differences decreased approximately linearly with η over η from `0.1` to `0.0005`. This validates one smooth symmetric example numerically, **not** the general interchanges.

## 7. Two-dimensional tangent-disk formulas (EXPLICIT BALL CALCULATIONS)

Use local `(s,z)` coordinates with `y=(0,0)`, inward direction `+z`, ball center `c=(0,r)`, `S` parametrized by

\[
x(\theta)=(r\sin\theta,\;r(1-\cos\theta)),\quad -\pi\le\theta\le\pi,
\quad |x(\theta)-y|=2r|\sin(\theta/2)|.
\]

Take `y'=(0,η)` and set `a=2 arcsin(d/(2r))` (`0<a<π`). `C` corresponds to `|θ|≤a`, `R` to `a<|θ|≤π`.

The **finite-η outward flux per unit arc length** is

\[
P_\eta(\theta)=\frac{2r\eta-\eta^2}
{2\pi r[\eta^2+2r(r-\eta)(1-\cos\theta)]}.
\]

Away from `θ=0`, its normalized limit is

\[
P^*(\theta)=\frac1{\pi |x(\theta)-y|^2}
=\frac1{4\pi r^2\sin^2(\theta/2)}.
\]

Integrating `ds=r dθ` on `R` gives

\[
\boxed{M_R^*=\frac1{\pi r}\cot\left(\frac a2\right)
=\frac{\sqrt{4r^2-d^2}}{\pi r d}.}
\]

Thus the density of the **angle** is

\[
p_\Theta(\theta)=\frac{rP^*(\theta)}{M_R^*}
=\frac1{4\cot(a/2)\sin^2(\theta/2)},\quad a<|\theta|\le\pi.
\]

**Inverse-CDF sampler:** Draw independent `U,V~Uniform(0,1)`. Set `cot(θ_abs/2)=(1−U)cot(a/2)`, equivalently `θ_abs=2 atan2(1,(1−U)cot(a/2))`; set signed `θ=θ_abs` if `V<1/2`, otherwise `−θ`. Jump to `X=x(θ)` translated and rotated to the actual tangent ball. Then `|X−y|≥d` by construction.

For **constant** `f` in a two-dimensional radius-`r` disk, the volume identity is exactly

\[
\int_B G_\eta f\,dx
=\frac f4[(r-\eta)^2-r^2],\quad
V^*= -\frac{fr}{2}.
\]

This computes the scalar limit **without** needing to exchange the limit and the integral, and is the justified analytic specialization used in the Stage C prototype. For `f=0`, the volume term vanishes. Spatially varying `f` requires a separate unbiased/consistent estimator and analysis.

## 8. Proof-of-concept manufactured PDE and numerical results (EMPIRICAL)

The Python test domain was the rectangle `Ω=(-1,1)×(0,1)`, with Neumann on its bottom (`z=0`) and Dirichlet on other three sides. The manufactured solution was

\[
u(s,z)=1+0.3z+0.2s^2+0.1z^2,
\quad f=\Delta u=0.6,
\quad g_N=\partial_{n_\Omega}u|_{z=0}=-0.3.
\]

At evaluation point `(0,0.1)`, the exact solution is `1.031`. At a near-bottom tangent step the practical **approximate** recursion was

\[
\boxed{\widehat u(y)=\widehat u(X)+\frac{g_N(y_0)-fr/2}{M_R^*},\qquad X\sim p_R.}
\]

The approximation **drops** `D_0=∇u(y)·n−g_N(y_0)` and `C*`, and therefore introduces bias. At ordinary interior WoS steps with sphere radius `ρ`, add the exact constant-forcing correction `−f ρ²/4`; terminate near Dirichlet boundary by projection and evaluating `g_D`. This is a single-chain recursion with additive known contributions, not branching.

Stages and evidential limits:

* **Stage A:** Verified the exact positive-η Green identity by deterministic quadrature on a disk for the manufactured polynomial; separately numerically checked three proposed integral limits in the symmetric test. This does *not* prove the general `B` or `C` limit interchanges.
* **Stage B:** Sampled the analytic `p_R` using the inverse CDF, compared its mean against deterministic quadrature; observed ordinary Monte Carlo `N^{-1/2}` scaling. `u(X)` was still evaluated from the known analytic function.
* **Stage C:** Replaced known `u(X)` by the recursive mixed-boundary WoS solver; observed termination in tested finite runs and solution estimates broadly near the analytical reference. Does not prove almost-sure termination, unbiasedness, or convergence under joint parameter refinement.

### Radius-selection heuristic used later

At each near-Neumann point choose proposed `r0` and inward center `c=y−r n_Ω(y_0)`. Query closest-boundary/signed interior clearance `q(c)`. Accept if center is inside and `q(c)≥r` (with numerical margin); otherwise set `r←r/2` and retry. Proposed radius is **not** necessarily maximal. Whenever halving causes `d≥2r`, the prototype used `d_eff=min(d,1.9r)`; this changes the effective cap and must be recorded in benchmarks. Geometry-query counts include initial point queries and center-clearance trials. In arbitrary geometry, the closest-point API must ensure inside/outside and normal information are well-defined; mere unsigned distance is insufficient to certify containment.

For the rectangle, the chosen **current baseline** is `d=0.01`, `r0=0.4` (equivalent to `k=r0/d=40`), halve on rejection; a nominal single boundary tolerance of `0.01` was used for Neumann activation and Dirichlet termination in the standalone Python tests. In the `d=0.01` radius sweep with `6000` walks per ratio, `k=40` gave roughly **77.1 mean steps and 120.5 geometry queries per walk**, solution estimate `1.0296`, exact `1.031`, Monte Carlo SE about `0.0041`. These are *specific empirical results*, not general performance guarantees. A separate two-parameter sweep explored `d∈{0.1,0.05,0.01}` and `tol∈{0.2,0.1,0.05,0.01}` with `r0=0.4`, displaying signed errors and sample SE. Near `tol≈d` some small observed signed errors arose, but this can reflect cancellation of biases and Monte Carlo noise; no proven optimal `tol/d` relationship exists.

**Crucial benchmark control:** When varying `d`, keep `r0=0.4` fixed if the purpose is cap-size convergence. Do not inadvertently enforce fixed `k=40` and change `r0` with `d`. When varying radius, hold `d` fixed. Where feasible, use distinct Dirichlet and Neumann tolerances and record cap fallback statistics.

## 9. Open mathematical problems and conjectures (NOT ESTABLISHED)

1. **Volume interchange:** rigorous sufficient conditions (regularity of `f`, ball dimension, boundary-layer integrability) for `lim∫_B (G_η/η)f=∫_B G*f`. Prove it, or avoid it using an exact finite-η calculation or a specifically justified source sampler.
2. **Cap-error interchange:** analyze `lim E_C(η)/η`, the singular flux `P*` on `C`, tangential odd-term cancellation, and any symmetric principal-value interpretation. Derive rigorous cap-error bounds in terms of `d`, `r`, and derivatives of `u`.
3. **Normal-displacement error:** bound `D_0=∇u(y)·n(y_0)−g_N(y_0)` in terms of `δ=dist(y,Γ_N)`, including curved boundaries and projection regularity.
4. **Global error:** control accumulation of dropped `C*` and `D_0` over multiple near-Neumann encounters. A small one-step error does not automatically imply a small global bias.
5. **Termination and variance:** almost-sure walk termination on a mixed domain; finite expected steps/geometry queries; control of estimator variance when the unit recursive weight is coupled to accumulated source and Neumann increments.
6. **Parameter scaling:** quantify jointly how `d`, `ε_N`, `ε_D`, `r`, and geometry affect bias, variance and cost. The earlier heuristic local scalings `normalized cap error ~ O(d²)` and `normalized displacement error ~ O(δd)` require hypotheses and are **not proved general/global bounds**; low measured error near `d=ε_N` may be bias cancellation.
7. **Higher dimensions:** derive and validate `p_R` and `M_R*` in 3D before enabling the 3D solver; 2D inverse CDF must not be reused unchanged.
8. **General forcing and geometry:** formulate a compatible volume term for variable `f`, a cap sampler/estimator as needed, and robust geometry for non-smooth features, concavity, and Neumann/Dirichlet junctions.
9. **Comparative novelty:** investigate cap-subtraction/regularization precedents in singular boundary integral equations and WoS literature before claiming that the construction is novel. The novel *combination* may be the candidate contribution, not any single algebraic ingredient.

**Present main research question:** Can the tangent-sphere cap decomposition be turned into a demonstrably convergent, efficient, nonbranching Monte Carlo method for mixed Dirichlet–Neumann Poisson problems, with quantitative bias bounds and a fair comparison against WoSt? The immediate mathematical obstacle is controlling the singular cap and volume limits and the neglected errors; the immediate software milestone is a small, validated CPU 2D solver integrated alongside, not inside, WoSt.

## 10. Implementation handoff / verification checklist

1. Build existing WoSX baseline and run an existing 2D WoSt test first. Do not change WoSt's numerical method.
2. Implement independent `WalkOnTangentSpheres<T,2>` using the same geometry, PDE, RNG and statistics interfaces; restrict initial PDE tests to `f=0` and constant `f`. Keep sign conventions explicit and verified against manufactured solutions.
3. Unit test 2D `P_η`, `P*`, `M_R*` and inverse-CDF sampling against high-precision quadrature and distributional tests. Check sample coordinates, chord cutoff, positive mass, and numerical stability when `d/r` is small.
4. Unit test finite-η exact identity, including `E_C` and `E_0`, before using the approximate recursion. Keep independent tests of volume and cap limit interchanges marked *experimental* unless proofs are supplied.
5. Separately test a Neumann step with analytic `u(X)` (Stage B), then with recursion (Stage C); inspect signed bias and sampling SE. Unit recursive *surface* weight should be 1.
6. Implement nearest-boundary point and signed center-clearance radius tests; count all geometry queries and rejections. Ensure the accepted ball is geometrically admissible; log effective cap radius (`d_eff`). Guard against `d_eff≥2r`, zero radius, infinite loops, and accidental crossing of Dirichlet features.
7. Use identical manufactured PDEs, geometry, start points, and reference solutions for WoSt and new solver. Compare matched accuracy (bias + sampling error) as well as matched number of walks, wall clock, steps, query types, query counts, Neumann encounters, and failure/termination rates. Do not treat one generic query as equivalent to all WoSt ray/silhouette operations.
8. Start with the documented baseline `d=0.01`, `r0=0.4`, reduction `1/2`, `ε_N=ε_D=0.01` in the original rectangle, then sweep one parameter at a time. Port the Python experiment and regress C++ statistics to its means *within expected Monte Carlo error*, not to identical random trajectories.
9. Keep clearly separated in code: exact finite-η representation; special-case proven formulae (`V*=-fr/2` for constant `f` in 2D); assumptions about general limit kernels; and approximate practical correction dropping `E_C` and `D_0`.

## 11. Deferred idea: generalized boundary-operator kernel (SPECULATIVE; NOT Construction 4)

A later discussion proposed a different zero-Dirichlet kernel `K` satisfying formally

\[
\Delta_x K(x,y')=\delta(x-y')-\eta\,\mathcal B_x\delta(x-y'),\quad K|_S=0,
\]

with `\mathcal B` encoding Dirichlet, Neumann, or Robin conditions. Distributional derivatives involve formal-transpose/sign and variable-coefficient issues. A speculative flux-mass split might produce a Robin absorption probability. This has **not** been mathematically established as a unifying Monte Carlo algorithm; it did not recover certain Dirichlet termination. The research scope was intentionally returned to the original zero-Dirichlet Green function `G_η` and Neumann Construction 4 for the first paper. Do not substitute `K` for `G_η` or build Robin behavior unless separately requested.

## 12. Source files supplied to Codex

* `construction_4_checkpoint_revised_2026-10-26.md` — authoritative saved mathematical state.
* `stage_a_construction4.py` — finite-η exact identity and special 2D disk limits.
* `stage_a_limit_interchanges.py` — numerical checks of three limit interchanges.
* `stage_b_single_step.py` — 2D inverse-CDF sampling validation.
* `stage_c_recursive_wos.py` — first independent recursive solver.
* `stage_c_radius_heuristic.py` — signed clearance query and halving algorithm.
* `stage_c_radius_sweep_extended.csv`, `stage_c_d_tol_sweep.csv`, `stage_c_reduction_all_five.csv` — empirical tuning results. They are diagnostic evidence, not proofs.

**Stop and report any inconsistency between a mathematical claim and the code rather than silently editing the derivation to fit the implementation.**
