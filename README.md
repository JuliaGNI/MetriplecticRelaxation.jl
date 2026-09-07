# MetriplecticRelaxation

Reproduction of the **two-dimensional** numerical examples of *Metriplectic relaxation to
equilibria* (CNSNS 110076), built on [PoissonBrackets.jl](https://github.com/JuliaGNI/PoissonBrackets.jl)
and [SimpleSplines.jl](https://github.com/JuliaDEC/SimpleSplines.jl).

| case | § | what it shows |
|:--|:--|:--|
| **A1–A4** | 4.1–4.2 | relaxation on the periodic torus, double bracket and projector bracket |
| **B1–B3** | 5.4 | reduced Euler on `[0,1]²` with homogeneous Dirichlet conditions |
| **C1–C2** | 5.5 | Grad–Shafranov equilibria, rectangle and mapped disk |

The three-dimensional Beltrami example of §6.4 is out of scope: it needs a de Rham complex of
compatible spline spaces that SimpleSplines does not have.

## A deliberate deviation from the paper

The paper discretises §4 with a **Fourier spectral method** and §5 with **P2 Lagrange elements in
FEniCS**. This repository uses **B-splines throughout**. The numbers therefore reproduce the
paper's *claims* — energy conserved to machine precision, monotone entropy decay, convergence to
the stated closed-form references and the rates ≈ 1 and ≈ 1/2, `τ_h = (ℓ_h/2π)²`, and the
eigenvalues `λ = 0.030302` and `λ = 0.002599` — and not its bit patterns. Nothing here should be
compared figure-pixel to figure-pixel with the published plots.

Where the paper leaves a choice open — time step, final time, stopping criterion, solver
tolerance, quadrature order — the choice is recorded in `CHANGELOG.md` next to the number it
produced.

**§4 is run twice**, once with B-splines and once with the paper's own Fourier spectral method,
because two discretisations sharing no code path is the only way to tell a statement about the
*equation* from a statement about the *method*. `scripts/converge.jl` turns their agreement into
a refinement statement: the spline difference from the spectral reference must fall at the
order of the spline space.

## §4 — what is implemented, and what it found

`src/torus.jl` holds the problem itself, free of any discretisation; `src/spectral.jl` and
`src/spline.jl` are the two solvers; `src/diagnostics.jl` the quantities both report.

All four runs reproduce their claims. The numbers below are the spectral runs at the paper's own
`256²` and `Δt`; `CHANGELOG.md` carries the spline column, the choices behind each, and the
full tables.

| run | claim | measured |
|:--|:--|--:|
| A1 | `H` conserved | 9.79e-16 |
| A1 | `τ_h = (ℓ_h/2π)²`, from the decay rate | rel **6e-5 … 1e-3**, `r² = 1.000000` |
| A1 | `u → u_∞(h)`, contour average conserved | 9.2e-8 … 4.7e-6 of peak |
| A1 | relaxation is **incomplete**: `ω(T) ≠ ω_η` | 5.95× its own norm away |
| A2 | `H` conserved | 1.66e-15 |
| A2 | `S` plateaus **above** `S_η` | **+179 %** of `S_η` |
| A3 | `S → S_η` | 2.35e-10 |
| A3 | `S − S_η` rate (**exact: 1**) | **1.00047**, `r² = 1.0000000` |
| A3 | `‖ω(t)−ω(T)‖` rate (**exact: ½**) | **0.50679**, `r² = 0.9999970` |
| A3 | spline vs spectral, final state | **1.42e-07** |
| A4 | both readings of `u₀` give rate ≈ 1 | 1.02186 / 1.02213 |

And the control: the spline difference from the spectral reference falls at **4.12, 4.08, 4.05**
for A1–A3 — the order cubic B-splines have — so the two agreeing is a refinement statement, not
a coincidence at one resolution.

Three things came out that are about the **paper** rather than about the reproduction. Each has a
script that settles it, and none is a change made silently.

- **A factor of 2 is missing from the §4.2 evolution equation as printed.** The equation below
  `eq:projector-brackets` reads `∂_t u = −[u − u_Ω − H(u)‖φ‖⁻²φ]`; the manuscript's own
  `eq:L2-projector` and `eq:Euler_H_periodic` give `2H` there, and only `2H` reproduces its own
  `eq:SS-projector`. With the printed factor the energy is not conserved, at 70 % of scale.
  The analytic test case one page earlier is unaffected, which is what makes it look like a
  transcription slip. — `scripts/verify_projector_factor.jl`

- **A4's initial condition as printed is discontinuous on the torus**, by 8.5 % of its own peak:
  `eq:initial_gaussian` is written unmodified on `T²` and A4's Gaussian sits `π/2` from the
  boundary with `w₂ = 1` and amplitude 1.8. A1 reaches `1.2e-25` there and A2/A3 `5.2e-5`, so
  A4 is the only run affected. Which reading is intended is not stated, so **both are run**, and
  neither changes a conclusion. Two independent controls: periodising the Gaussian improves the
  spline/spectral agreement on `S(0)` by **51 336×** (9.73e-05 → 1.90e-09), and it restores the
  spline convergence order from **1.10** to **4.05** where every other run is at 4.
  — `scripts/verify_spline.jl`, `scripts/converge.jl`

- **A4's early entropy relaxation rate is *faster* than its late one**, where Fig. 7's discussion
  says slower: **1.248** early against **1.022** late, agreeing to four digits between the two
  discretisations. The linearised spectrum predicts it — `cos(2x₂)` is a `λ = 4` mode, decaying
  at `3/2` against the `λ = 2` mode's `1`. Not necessarily an error: the manuscript's reasoning
  is about the PL constant of `eq:shrunk-cone`, a *lower bound* an actual rate may exceed. What
  does not reproduce is reading it as a description of the observed rate. — `scripts/run_a4.jl`

And one that sharpens a claim rather than correcting it: the manuscript reads its relaxation
rates off a semi-log plot as "≈ 1" and "≈ 1/2". Linearising the projector flow about a relaxed
state gives a spectrum of **exactly** `{1 − 1/λ}` over the eigenvalues of `−Δ`, verified against
the Jacobian to ten digits. So `≈ 1/2` is the `λ = 2` mode exactly, `≈ 1` follows because
`S − S_η` is quadratic, and the whole `λ = 1` eigenspace is neutral — which is
`eq:u-eta_Euler_periodic`'s three-parameter family. — `scripts/verify_projector_rates.jl`

`τ_h` likewise has a closed form the manuscript does not give: on the contours of
`h = cos²x₁ sin²x₂` the integral `eq:relaxation-time` is a complete elliptic integral, so
`ℓ_h = π/(√h·agm(1,√h))` and `τ_h = 1/(4h·agm(1,√h)²)`. A1 therefore measures the decay rate
against a closed form rather than against another numerical procedure.

## §5.4 — what is implemented, and what it found

`src/euler.jl` holds reduced Euler on `[0,1]²` with homogeneous Dirichlet conditions, under the
collision-like div–grad bracket. One file rather than §4's problem-plus-two-solvers: this bracket
has no spectral counterpart in this tree, so the cross-discretisation control §4 has is simply not
available. `ImplicitMidpoint` **is** Crank–Nicolson for this field, so each step is a Newton solve
through PoissonBrackets' existing `NewtonSolver`; no solver code is written here.

`scripts/verify_euler.jl` settles the problem, the space, the bracket and both references in
**101 checks** before any run is believed — including that §5 states `w²` and not `w`, with the
`w`-reading control that must fail, and the mismatched-generator and indefinite-mobility controls
that must break degeneracy and positivity respectively.

All three runs use **26² cells of degree 2** and not the paper's `64²` P2, and the reason is
measurable: this bracket is nonlocal, so its operator is dense and one Newton matrix is `N` dense
`N×N` assemblies. Measured, one step takes **1.0 s** at `N = 256`, **7.7 s** at `N = 576` and
**36 s** at `N = 1024`; `64²` would be `N = 4096` and some 40 minutes *per step*. Nothing §5.4
claims needs that: the discrete first Dirichlet eigenvalue is already within `3.4e-5` of `2π²` at
**8** cells. What does need resolution is the narrow direction, `w₁ = 0.1`, and 26 is where the
`L²` projection of the initial Gaussian stops oscillating below zero — measured, not chosen.

| run | claim | measured |
|:--|:--|--:|
| B1 | `H` conserved | 2.45e-12 relative, **1.37e-15 absolute** |
| B1 | `S` monotone | worst step **−3.31e-10** of `S₀` |
| B1 | `S → S_η = λ₁,₁H₀` | excess **+6.56e-09** |
| B1 | `ω → λ₁,₁φ`, fitted `λ` | **19.7392146901** vs `λ_h = 19.7392146644`, `2π² = 19.7392088022` |
| B1 | ...and the residual `‖ω−λφ‖/‖ω‖` | **7.25e-05**, from 6.27e-01 at `t = 0` |
| B2 | `S` flat, **then** monotone | held to **8.10e-03** of the fall for `t < 8.5`, then fell 99.997 % |
| B2 | ...against B1 as control | B1 breaks at the **first** sample, holding nothing |
| B2 | the instability, `(S,S)` | `1.21e-05 → 1.30e-02`, a factor **1070** |
| B2 | relaxes to the **lowest** mode | `λ: 510.84 → 19.74213`, `λ_h = 19.73921` |
| B1, B2 | `|Δλ|/λ_h = 0.2474 rel²` | **0.24739** and **0.24740** — a measured law, five digits |
| B3 | the `Δt` threshold | monotone to **0.02**, inadmissible from **0.04**; tightens with the mesh |
| B3 | `ω → e^{λφ−1}`, `λ = (M+S)/2H₀` | `λ = 15.0344`; residual **4.24e-03** from 4.56e-01 at `t = 0` |
| B3 | ...and the multiplier it assumes away | fitted `μ = +0.0071`, `μ/λ = 4.7e-04` |
| B3 | the boundary layer, at margins `0,h,2h,4h` | `3.56e-02, 9.33e-03, 4.24e-03, 9.12e-04` |
| B3 | the two-space control (16², 25 steps) | `μ = −0.35` in `V_D` against **−1.45** in `V` |
| B1–B3 | `H` conserved | `2.45e-12`, `3.21e-14`, `1.39e-14` relative |
| B1–B3 | `S` monotone, worst step | `−3.31e-10`, `−5.75e-07`, `−8.56e-11` of `S₀` |

### The two findings, and both are about the manuscript

- **§5.4's closed-form references are the `μ = 0`, `c = 0` case of a four-multiplier equilibrium,
  and only a homogeneous-Dirichlet *state* space forces that.** The equilibrium condition is
  `δS/δu = λφ + c·x + μ`, with `μ` and `c` the multipliers of the bracket's other Casimirs — the
  mass `∫u` and the momenta `∫x_k u`, each a Casimir because its variational derivative has
  vanishing gradient. In `V_D` the constants and the linears are not in the space, the discrete
  Casimirs are absent, and both references become exact; in a space that contains them the mass
  is conserved and both acquire an extra term. The reproduction runs both and measures the
  difference — `run_b3.jl`'s control table — and pays for the exact references with a mass that
  drifts 46 % over B1.

- **B3's printed initial condition is not admissible for its own entropy.** `s = y log y` is
  undefined at `y ≤ 0` and `eq:M-condition` gives it the mobility `M = y`, which the same equation
  requires to be positive. The printed Gaussian is strictly positive as a *function* — but with
  `w₁² = 0.01` it decays to `1.4e-11` of its peak at the far corner, which a Galerkin scheme
  cannot tell from zero: the projected state starts at `6e-13` and the **first** step drives the
  iterate to `-6.6e-5`, where the entropy raises a `DomainError` before any step completes. B3
  therefore carries a positive background of 1 % of its peak.
  The continuum does not have this problem — `M = ω` vanishes with `ω`, so the equation is
  degenerate-parabolic and positivity-preserving, and it is the Galerkin projection that loses it.

## §5.5 — what is implemented, and what it found

Two files, because this section has a classical counterpart §5.4 lacks. `src/takeda.jl` holds the
Grad–Shafranov problem itself and the **classical** solver — the finite-volume `−Δ*` and the
iteration of Takeda & Tokuda Eqs. (2.111)–(2.112) — plus C2's geometry: `disk_map`, its
triangulation, and `disk_eigenvalue`. `src/gradshafranov.jl` holds the **relaxation** on the
rectangle: `GradShafranovBox`, the entropy and its references, `gs_rayleigh`, `gs_fit`, and
`SECTION55_RUNS`. Both feed §4's and §5.4's diagnostics, so the numbers are comparable across
sections by construction.

`scripts/verify_takeda.jl` settles the classical side in **38 checks** across seven sections, ten
of them controls; `scripts/verify_gradshafranov.jl` settles the relaxation in **42 checks** across
ten sections, also ten controls. Between them they pin the measure, the profile, the mobility and
both references before `run_c1.jl` is believed.

C1 runs at **18 × 21 cells of degree 2** with `Δt = 0.0625` over 400 steps. `Δt` is not a
preference: measured on the run's own space, `‖Δt − Δt/2‖ = 2.738e-04` against
`‖Δt/2 − Δt/4‖ = 7.390e-05` is a ratio of **3.705**, implicit midpoint's second order. The first
pass of this run used `Δt = 0.25`, passed every other check, and was still wrong about the
trajectory — `H`, monotonicity and the equilibrium are all `Δt`-independent, so nothing else in
the output would have caught it.

| run | claim | measured |
|:--|:--|--:|
| C1 | `H` conserved | **5.324e-12** relative over 400 steps, `Δt`-independent by derivation |
| C1 | `S → λ_h H₀` | `(S_T − λ_h H₀)/λ_h H₀ = **−5.3172e-12**` |
| C1 | `u/(Cr²+D) → λψ` | **6.9980e-05**, from `6.3955e-01` at `t = 0` — the space's own floor |
| C1 | `S/H = λ_h`, the sharp arrival | `0.03023479920660` vs `λ_h`, to **1.95e-15** |
| C1 | ...and the fit's second-order law | `\|Δλ\|/(λ_h·rel²) = **0.3463**` |
| C1 | against Takeda's `64²` iteration | `λ_h = 0.030234799207` vs `0.0302248584`, **+3.29e-04** |
| C1 | against the continuum eigenvalue | **+5.73e-06** |
| C1 | against the printed `0.030302` | **+2.22e-03** — the paper's own grid error, same side |
| C1 | the two-space control (10 × 12, 240 steps) | `V` misses `eq:gs-ref` by **2.5171e-01** against `V_D`'s **4.2761e-04**, a factor **589** |
| C2 | the reference eigenvalue, `P₁` disk | `λ = 0.0025990851` against the printed `0.002599` |

### The findings, and both are again about the manuscript

- **The manuscript prints two profiles, and mistaking them costs a factor 4.6**; and its printed
  `λ` values encode its own mesh rather than physics — `0.030302` sits `+2.23e-03` above the
  continuum eigenvalue, and the reproduction disagrees with the paper by essentially exactly
  that, from the same side. Every measured `λ` here lies **above** its continuum limit, as a
  conforming Galerkin discretisation requires.

- **C1's mass drifts by +83.7 %, and that is the price of the exact reference.** The constant
  function is not in `V_D`, which is the same absence that forces the equilibrium's multipliers
  to zero and so makes `eq:gs-ref` exact. §5.4's B1 drifts +46 % for the same reason. The
  cheaper eigen-based control reports nothing at all — the pencil `(MΛ, W)` carries no mass
  constraint, so its lowest eigenvector is the `μ = 0` member in *both* spaces; what separates
  them is the **flow**, not the fixed-point problem.

### C2's relaxation is deferred, and nothing was regularised

Two obstructions, either sufficient alone. The map of `eq:mapping` collapses the whole circle to
one point at `s = 0` — measured, the spread of `disk_map(0, θ)` over sixteen angles is exactly
`0.000e+00` — so a tensor-product spline space is not even `C⁰` there; imposing single-valuedness
needs a **polar-spline** construction that SimpleSplines does not have. And PoissonBrackets'
only two-dimensional space is `TensorSplineSpace`, while `CollisionBracket` reads `∇ψ` at *its*
quadrature points, so the `P₁` triangulation that gives the reference eigenvalue cannot carry the
flow. There is no `ε` floor on `s`, no puncture at the origin, and no modified basis: C2's
reference eigenvalue is reproduced and its relaxation is simply not run.

## Layout

| directory | holds | committed |
|:--|:--|:--|
| `src/` | the model definitions: domains, brackets, entropies, references, diagnostics | yes |
| `scripts/` | **all** code: run drivers, shared includes, option parsing | yes |
| `runs/` | **all** data a run produced (`.jld2`, `.h5`, `.csv`) | no |
| `results/` | **all** rendered output (`.pdf`, `.png`, reports) | no |

`runs/` and `results/` are regenerable and are not tracked; `/runs` and `/results` are
root-anchored in `.gitignore`.

## Running

```sh
julia --project=scripts scripts/run_all.jl          # every run driver, in order
julia --project=scripts scripts/run_all.jl run_a1.jl  # one of them
```

Every driver takes `--runs-dir` and `--results-dir`, so the same script can write into a scratch
tree or a talk's figure directory without being edited.

The test suite runs coarse, fast versions of the same claims:

```sh
julia --project=. --check-bounds=auto test/runtests.jl
```

`Pkg.test()` also works, but it defaults to `--check-bounds=yes`, which inflates the timings
roughly fourfold and invalidates any precompile image shared with a run driver. Prefer the form
above when the numbers or the wall-clock matter.

## Dependencies

`PoissonBrackets` and `SimpleSplines` are unregistered, so `Project.toml` carries a `[sources]`
table for both — a Julia 1.11 feature, which is why the `julia` floor here is 1.11 rather than the
1.10 used elsewhere in this tree.

Both are resolved from their **git remotes at `rev = "main"`**, not from the sibling checkouts in
`../../Packages/`. A relative path escapes the repository and so does not exist on a CI runner:
with the path form every job failed in `buildpkg` before a single test ran. The cost is that
resolution follows `main` on the remote, so a local change to either package is invisible here
until it is pushed — develop against a sibling with `Pkg.develop(path=...)` in a *local*
environment, and do not put the path back into `Project.toml`. Its comments say the same thing
at length, and they are load-bearing: `Pkg.instantiate` and `Pkg.develop` both strip them.
