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
julia --project=. -e 'using Pkg; Pkg.test()'
```

## Dependencies

`PoissonBrackets` and `SimpleSplines` are unregistered, so `Project.toml` carries a `[sources]`
table for both — a Julia 1.11 feature, which is why the `julia` floor here is 1.11 rather than the
1.10 used elsewhere in this tree. `PoissonBrackets` is resolved from the local sibling checkout at
`../../Packages/PoissonBrackets`.
