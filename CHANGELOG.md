# Changelog

All notable changes to MetriplecticRelaxation.jl are recorded here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and the project adheres to
[Semantic Versioning](https://semver.org/spec/v2.0.0.html).

This repository begins with this file, so there is no history predating it and nothing has been
reconstructed from `git log`.

Two rules specific to an experiment repository, and they are the reason this file matters more
here than in a library:

- **When a measured number changes, the entry says why** — a fix, a dependency bump, a different
  tolerance, a different machine. A results table that silently differs from last month's is the
  failure this prevents.
- **Every choice the paper leaves open is recorded next to the number it produced.** The published
  §5 runs state no time step, no final time, no stopping criterion and no solver tolerance; each
  one picked here is a choice, not a datum, and is written down as such.

## [Unreleased]

### Added

- **`src/torus.jl` — the §4 problem definition.** The periodic domain, the anisotropic Gaussian
  of `eq:initial_gaussian`, the prescribed Hamiltonian `h = cos²x₁ sin²x₂` of `eq:islands-h`,
  the closed-form entropy minimisers of both test cases, and the contour geometry behind the
  relaxation time. Deliberately free of any discretisation, so that the spectral and the spline
  solver are two independent solvers of the *same* stated problem and their agreement is
  evidence about the equation rather than about a shared code path.

  **`τ_h` has a closed form, which the manuscript does not give.** `eq:relaxation-time` defines
  `ℓ_h = ∮ ds/|∇h|` and says it can be estimated "from a sample of points on the considered
  contour", obtained by integrating the flow of `X_h`. For this particular `h` it can be had
  exactly. On the island centred at `(π, π/2)`, writing `s = x₁−π` and `t = x₂−π/2` gives
  `h = cos²s cos²t`, and the substitution `sin s = √(1−h) sin θ` maps the contour integral onto
  the complete elliptic integral of the first kind:

      ℓ_h = (2/√h) K(1−h) = π / (√h · agm(1,√h))     τ_h = 1 / (4h · agm(1,√h)²)

  Both limits the manuscript describes in words fall out: `τ_h → 1/4` at the island centre —
  which is also the period `π` of the harmonic approximation `h ≈ 1−s²−t²`, whose flow rotates
  at angular frequency 2 — and `τ_h → ∞` logarithmically at the separatrix. This is an exact
  reference where the paper has a sampled estimate, so A1's `τ_h` check is against a closed form
  rather than against another numerical procedure.

  **Recorded choice — the final time `T`, which the manuscript never states.** `Δt` is the
  paper's own (`1e-4` for A1, `1e-3` for A2–A4). For the reduced Euler runs Fig. 6 places the
  vertex of the cone at `t ≈ 5`, so `T = 10` is taken, leaving as much trajectory again past the
  vertex to measure the asymptotic rates on. A1 has no such landmark and is given `T = 20`,
  which is 80 `τ_h` at an island centre and resolves the island interiors while leaving the
  separatrix stalled — the behaviour Fig. 1 shows.

  **Recorded choice — the Gaussian is not periodised.** The manuscript writes
  `eq:initial_gaussian` unmodified on `T²`, with no summation over images, so this does too. The
  resulting jump across the boundary is `1e-69` for §4.1's `w₂ = 0.4` but `5.2e-5` for the
  `w₂ = 1` of the reduced Euler runs. Both discretisations absorb it the same way, so it does
  not bias the spectral-versus-spline comparison, but it does put a floor under how well either
  represents the stated initial condition.

- **`scripts/verify_torus_geometry.jl`** — 24 checks on that geometry, run before any solver was
  built on it. The closed form agrees with quadrature of the arclength integral to **9e-16 …
  8.7e-15** over `h ∈ [1e-5, 0.9]`; `ℓ_h` is confirmed to be the *period* of the `X_h` flow, the
  manuscript's own characterisation, with the orbit closing to **1.6e-14**; and the contour
  average reproduces a field that is already a function of `h` to **1e-12** and agrees with the
  flow-time average to **8.6e-14**. Two controls that must fail do fail: `agm(1,h)` in place of
  `agm(1,√h)` is off by **16 %**, and after half a period the orbit has moved **1.57**.

- **`src/spectral.jl` — the Fourier solver, which is the manuscript's own discretisation.**
  `poisson_periodic`, `canonical_bracket`, `hamiltonian_field`, the trapezoidal `integrate`,
  the two §4 vector fields, and classical RK4. This is the reference half of the pair; the
  spline Galerkin half shares no code path with it, so their agreement is evidence about the
  equation rather than about either method.

  **The paper's own `scripts/torus.jl` could not be depended on as a module**, contrary to the
  plan's expectation. It lives in a manuscript repository with no package boundary and no UUID,
  so `using` it would mean a hard-coded absolute path in committed code or a copy of the file.
  (The plan also names it `scripts/src/torus.jl`; the actual path has no `src/`.) The same
  operators are therefore built here and **checked against PoissonBrackets' `spectral_grid`**,
  which computes derivatives from a dense differentiation matrix `D = Re(V diag(ik) W)` rather
  than by FFT — a genuinely independent implementation. They agree to **8.2e-14 … 1.0e-13**.

  **Recorded choice — the state variable is `ω = u − u_Ω`, not `u`.** This is the manuscript's
  own alternative, stated at `eq:Poisson-eq-periodic` ("equivalently, we could have chosen the
  phase space to be the subspace of functions satisfying `u_Ω = 0` and `u = ω`"), and it is what
  makes the projector bracket well posed: its generating field `φ` is mean-free by construction.
  Both vector fields preserve the zero mean exactly — measured at **1.2e-18 … 1.7e-17**
  normalised across all four runs — so `u_Ω` is a constant carried alongside and added back only
  for plotting.

- **`scripts/verify_spectral.jl`** — 27 checks: the operator agreement above, the Poisson solve
  (`−Δφ = ω` on eigenmodes to **5.6e-16**, `φ` mean-free to **1.5e-18**), the identity
  `[h,[h,u]] = ∇·(X_h⊗X_h ∇u)` of `eq:parallel-diffusion` to **4.7e-15** with `∇·X_h = 0` at
  **4.8e-15** and `X_h·∇h = 0` **exactly**, and the semi-discrete conservation laws: energy
  conserved at **1.7e-17** (double) and **1.3e-16** (projector) normalised, entropy strictly
  decreasing, and `dS/dt = −∫|X_h·∇ω|²` matching to 11 digits.

### Found

- **A factor of 2 is missing from the §4.2 evolution equation as printed.** The equation just
  below `eq:projector-brackets` reads `∂_t u = −[u − u_Ω − H(u)‖φ‖⁻²φ]`. The manuscript's own
  definitions give `2H` there: `eq:L2-projector` sets `c(u,v) = ‖φ‖⁻²(φ,v)`, and with
  `v = δS/δu = ω` and `eq:Euler_H_periodic`'s `H = ½(φ,u)` this is `c = 2H/‖φ‖²`.

  This is **internal to §4.2** — the two printed equations cannot both be right. Substituting
  the same `c` into `(S,S) = (ω, Π_H ω)` reproduces `eq:SS-projector`,
  `(S,S) = 2S − 4H₀²/‖φ‖²`, only with the factor 2.

  `scripts/verify_projector_factor.jl` settles it numerically. With the printed factor the
  energy is **not** conserved — `dH/dt = −9.51`, i.e. **70 %** of `‖φ‖‖f‖` — while with the
  derived factor it is **9.4e-17**; and `−dS/dt` matches `eq:SS-projector` to **8.1e-16** with
  the factor 2 against a **192 %** relative error without it. The analytic test case one page
  earlier is the mirror image and is unaffected: `eq:analytical_H` has no `½`, so the printed
  factor is correct there (**7.7e-17**, against **38 %** for the factor 2). That asymmetry is
  what makes this look like a transcription slip rather than a different convention — a
  convention would have moved both.

  **Why it would survive review unnoticed:** the factor-1 field is still symmetric and still
  dissipates entropy (`dS/dt = −16.1` against the correct `−5.50`), so it looks healthy in both
  plots Fig. 4 shows. What it does is relax to a state of the *wrong energy*, hence to the wrong
  member of the family `eq:u-eta_Euler_periodic`. This reproduction uses the derived factor, and
  the discrepancy is a candidate finding for the manuscript rather than a change made silently.

### Fixed

- **A catastrophic cancellation in the contour quadrature, found by its own refinement check.**
  The direct parameterisation forms `sin t = √(1 − h/cos²s)`, whose bracket is a difference of
  two nearly equal numbers at the turning points, where `cos²s → h`. Near the separatrix that
  cost most of the mantissa: at `h = 1e-4` the residual against the closed form was **1.3e-10**
  and *grew* with the node count — `1.6e-11` at `n = 200` rising to `4.6e-9` at `n = 3200` —
  which is the signature of round-off rather than of truncation, and is what identified the
  cause instead of prompting a wider tolerance.

  The cure is algebraic and exact: `cos²s = h + (1−h)cos²θ` and `sin t = √(1−h)|cos θ|/cos s`
  are both sums of positive terms with no cancellation anywhere. The residual fell to
  **2.1e-15** at `h = 1e-4` and now stays at round-off down to `h = 1e-5`, so the check runs at
  a `1e-13` tolerance rather than the `1e-10` that the defective version could not meet.

- **The repository itself.** Scaffolding only: `Project.toml` with the dependency set the
  reproduction needs and hand-written `[compat]` bounds, the `src/` – `scripts/` – `test/` layout
  with root-anchored `/runs` and `/results` in `.gitignore`, the `Checks` PASS/FAIL harness copied
  from PoissonBrackets' `scripts/check.jl`, an empty `scripts/run_all.jl` driver list, and the
  shared CI workflows and git hooks installed from `~/Research/Knowledge/AI/githooks/`. No physics,
  no brackets and no run drivers yet.

  The `julia` floor is **1.11**, not the 1.10 used across the rest of the tree: `PoissonBrackets`
  and `SimpleSplines` are both unregistered and are reached through a `[sources]` table, which
  Julia 1.10 ignores rather than honours, so a 1.10 environment fails to resolve rather than
  resolving to something older.

  Both sources are **paths to the local siblings**, not `{rev = "main", url = …}`. That follows
  the rule `PoissonBrackets/scripts/Project.toml` already states — the url form resolves against
  a cached clone and pins the environment to whatever was last fetched — and it is also the only
  form that resolves here at all: Julia 1.13.0-rc3's `Pkg` segfaults (signal 11) in
  `Updating git-repo` for a git-url source, offline as well as online and with the cached clone
  already at the current commit.
