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
