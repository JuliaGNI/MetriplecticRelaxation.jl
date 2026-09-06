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

Two things came out that are about the **paper** rather than about the reproduction. Both have a
script that settles them, and neither is a change made silently.

- **A factor of 2 is missing from the §4.2 evolution equation as printed.** The equation below
  `eq:projector-brackets` reads `∂_t u = −[u − u_Ω − H(u)‖φ‖⁻²φ]`; the manuscript's own
  `eq:L2-projector` and `eq:Euler_H_periodic` give `2H` there, and only `2H` reproduces its own
  `eq:SS-projector`. With the printed factor the energy is not conserved, at 70 % of scale.
  The analytic test case one page earlier is unaffected, which is what makes it look like a
  transcription slip. — `scripts/verify_projector_factor.jl`

- **A4's initial condition as printed is discontinuous on the torus**, by 8.5 % of its own peak:
  `eq:initial_gaussian` is written unmodified on `T²` and A4's Gaussian sits `π/2` from the
  boundary with `w₂ = 1` and amplitude 1.8. A1 reaches `1.2e-25` there and A2/A3 `5.2e-5`, so
  A4 is the only run affected. Which reading is intended is not stated, so **both are run**.
  — `scripts/verify_spline.jl`, `scripts/converge.jl`

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
