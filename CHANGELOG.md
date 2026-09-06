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

  **Recorded choice — the final time `T = 10`, which the manuscript never states.** `Δt` is the
  paper's own (`1e-4` for A1, `1e-3` for A2–A4). For the reduced Euler runs Fig. 6 places the
  vertex of the cone at `t ≈ 5`, so `T = 10` leaves as much trajectory again past the vertex to
  measure the asymptotic rates on. A1 has no such landmark; `T = 10` is 40 `τ_h` at an island
  centre and resolves the island interiors while leaving the separatrix stalled — the behaviour
  Fig. 1 shows. A1 still costs the most, because its `Δt = 1e-4` is set by the explicit
  stability limit rather than by accuracy and buys it ten times the step count.

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

- **`src/spline.jl` — the B-spline Galerkin solver, which is the deliberate deviation.** The
  same equations in a periodic tensor-product spline space, through PoissonBrackets'
  `TensorSplineSpace`, `DoubleBracket`, `ProjectorBracket` and `MetriplecticFlow`. The numbers
  do not match the spectral run bit for bit and must not be reported as if they did; what is
  reproduced is the claims.

  `PoissonMap` applies `Λ = K⁻¹M` as a cached sparse factorisation rather than as an assembled
  matrix. `Λ` is dense — at `64²` cells that is 4096×4096, 134 MB, streamed once per
  Runge-Kutta stage, and a `128²` run would need 2.1 GB and be impossible. `K` is singular on
  the torus, and the obvious regularisation `K + bbᵀ/|Ω|` destroys the sparsity because
  `bᵢ = ∫Φᵢ` is dense; the **bordered** system adds one row and column instead. Its Lagrange
  multiplier comes out at **3.3e-16**, which is the check that the state really is mean-free.

  `LinearHamiltonian` and `EllipticEnergy` are two `DiscreteHamiltonian`s PoissonBrackets does
  not carry. A1's energy `H = (h−h_Ω, u)` is *linear*, and the package has no linear
  Hamiltonian: `MassCasimir` has exactly the right structure but its docstring defines it as
  the total mass, so using it would put a false statement in the code. `EllipticEnergy` exists
  because `QuadraticHamiltonian(M*Λ)` would assemble the dense matrix `PoissonMap` avoids, its
  constructor checking symmetry with `isapprox(A, A')`. Both are candidates for the package.

  **A1 does not go through `vectorfield`.** Its `h` is prescribed, so `X_h` never changes and
  the whole operator is a constant sparse matrix — **8.5 %** filled at `24²` — reducing the
  right-hand side to `−M⁻¹Aω̂`. That matters because A1 runs at `Δt = 1e-4` and so takes ten
  times the steps of the other three: at ~10 ms per generic evaluation its `4·10⁵` stages would
  be over an hour. The assembled path is checked against `vectorfield(MetriplecticFlow)` and
  agrees to **9.2e-15**, which is what makes it a substitution rather than a second
  implementation.

- **`scripts/verify_spline.jl`** — 44 checks. All four brackets are symmetric, positive
  semi-definite, and degenerate on the flow's own energy at **8.1e-16 … 1.3e-15**, giving
  `dH/dt` at **2.8e-17 … 3.7e-16** normalised with `dS/dt < 0` throughout. `MΛ` is symmetric to
  **3.7e-15**. The spline and spectral vector fields agree to **3.5e-4 … 4.8e-3** for A1–A3.

- **`src/diagnostics.jl` — the quantities every run reports, written once against both
  discretisations.** `Diagnostics` is the only place that knows which energy a run has, so the
  drivers and figures never branch on it. Energy, entropy, vorticity mass, `‖φ‖²`, the `Trace`
  time series, the cone residuals of `eq:theoretical-limits`, exponential rate fitting and the
  scatter data.

  **`best_fit_euler` has a closed form; the manuscript's three-phase search is not needed.**
  Expanding `eq:u-eta_Euler_periodic` gives `ω_η = A(a₁cos x₁ + b₁sin x₁ + a₂cos x₂ + b₂sin x₂)`
  with `A = √H₀/π` fixed by the energy and `a₁²+b₁²+a₂²+b₂² = 1`. The three phases parameterise
  exactly the unit sphere in those four coefficients, and the modes are `L²`-orthogonal on `T²`,
  so `‖ω − ω_η‖² = ‖ω‖² − 2A v·p + 2π²A²` is minimised over `‖v‖ = 1` at `v = p/‖p‖` — the
  normalised projection.

  Checked against the search it replaces: a `60³` scan over `(θ₀,θ₁,θ₂)` cannot beat it, coming
  within **2.5e-4** relative, and doubling the scan resolution to `120³` cuts the gap **30×**
  (**1.53e-4 → 5.04e-6`), which is what identifies the gap as the scan's own resolution rather
  than a defect in the closed form. A state already on `𝔠_η` is fitted to **3.4e-16**.

  **Recorded choice — the rate-fit window is `t ∈ [0.5T, 0.9T]`.** The manuscript's own reading
  of Fig. 6 is that the rate is asymptotic: for the second initial condition "the initial
  entropy relaxation rate is slower, but approaches ≈ 1 as the trajectory approaches the vertex
  of the cone". Fitting from `t = 0` would average the transient with the asymptote and report
  neither. `fit_rate` also returns `r²`, because a rate quoted without it says nothing about
  whether the decay was exponential at all, and drops samples below `1e-12` of the peak so that
  a round-off tail does not drag the slope toward zero.

- **`scripts/verify_diagnostics.jl`** — 30 checks, including the fit-vs-scan comparison above,
  exact recovery of known exponential rates (**1e-12**), an algebraic decay correctly showing
  `r² = 0.9957`, and the control that a *rising* entropy is caught. The spline and spectral
  energies agree to **1.2e-7** for A1–A3.

- **`scripts/run_a1.jl` … `run_a4.jl`, `runner.jl`, `projector_run.jl`, `converge.jl`,
  `figures.jl`.** Each driver takes `--runs-dir` and `--results-dir` and derives no output path
  from `@__DIR__`. Runs are serialised to `runs/` with `Serialization` — a stdlib, needing no
  `[compat]` bound, for data that is regenerable and read by nothing outside this repository —
  and each driver writes its own numbers to `results/<name>.md`, which is where prose quotes
  them from.

  **Recorded choice — A1 gets 128 spline cells, the others 64.** A1's Gaussian is the narrowest
  of the four (`w₁ = 0.25`), and the `L²` projection error of `u₀` onto the spline space is
  `1.05e-2` at 32 cells, `6.39e-4` at 64, `9.94e-5` at 96 and `3.11e-5` at 128. The contour
  deviations A1 tracks fall to ~`1e-4` of the peak before `T`, so 64 cells would put the
  measurement inside its own discretisation error.

  **Recorded choice — A1's `τ_h` contours are `h ∈ {0.5, 0.4, 0.3, 0.2}` on the upper island,
  and the band is forced.** The Gaussian sits *on* the separatrix (0.1 above it), so the low-`h`
  contours carry the amplitude and the high-`h` ones near the island centres are empty. Measured
  on `u₀`: the along-contour deviation is `1.2e-4` at `h = 0.9` but `T/τ_h = 34`, against
  `2.3e-1` at `h = 0.05` where `T/τ_h = 0.6`. Above `h ≈ 0.6` the "decay" would be
  discretisation noise; below `h ≈ 0.15` nothing relaxes within `T`, which is the manuscript's
  own point about the separatrices. The fit window is `[1.5τ_h, 4.5τ_h]` per contour, since
  `τ_h` varies by a factor 3.7 across them.

### Changed

- **A1's spectral right-hand side hoists `X_h` out of the time loop.** `h` is prescribed and
  fixed there, but the nested-bracket form recomputes `∂₁h` and `∂₂h` on every evaluation —
  four of its eight derivative applications. `parallel_diffusion` takes the precomputed field
  instead. The two agree **exactly** (0.00e+00), and A1's `256²` spectral run drops from
  **69 minutes to 41**. That matters only for A1, whose `Δt = 1e-4` is set by the explicit
  stability limit and buys it ten times the step count of the other three.

- **The projector bracket's relaxation rates are exact, not approximate.** The manuscript reads
  them off a semi-log plot: "exponential relaxation of entropy with exponential rate ≈ 1", and
  "the fact that ω has relaxation rate ≈ 1/2 is a consequence of the simple choice of the
  entropy function". Both can be derived.

  Linearising the projector flow about a relaxed state `ω* = φ*` — and keeping in mind that `H`
  is a function of `ω` and varies under an arbitrary perturbation, even though it is conserved
  *along* the flow — gives `∂_t ε = −[ε − Λε − (δc)φ*]` with
  `δc = [(φ*,ε) − (φ*,Λε)]/H₀`. For `ε` in the eigenspace of `−Δ` with eigenvalue `λ ≠ 1`, `δc`
  vanishes and

      ∂_t ε = −(1 − 1/λ) ε .

  `scripts/verify_projector_rates.jl` checks this against the Jacobian and gets `1 − 1/λ` to
  **10 digits** for λ = 2, 4, 5, 8, 9, with the off-mode residual at `1.5e-11 … 5.2e-11`. So
  `≈ 1/2` is the λ = 2 mode **exactly**, and `≈ 1` follows because `S − S_η` is quadratic in the
  perturbation. The whole λ = 1 eigenspace is **neutral** — that is `eq:u-eta_Euler_periodic`'s
  three-parameter family, and the manuscript's "minimally degenerate" made quantitative.

  Two consequences for how the rates are measured, both recorded:

  - **`T = 20` rather than 10 for A2–A4.** The λ = 4 mode decays at 3/4, only 1/4 faster than
    λ = 2, so it contaminates a vorticity fit as `e^{−t/4}`: **8.2 %** of the signal at `t = 10`
    and **0.67 %** at `t = 20`. A `T = 10` run measured **0.597** for what is exactly 1/2.
  - **The vorticity rate is fitted *earlier* than the entropy rate**, over `[0.35T, 0.55T]`
    against `[0.65T, 0.85T]`. `‖ω(t) − ω(T)‖` uses `ω(T)` as a stand-in for the limit — the
    manuscript's own construction — but `ω(T)` is not the limit, and the resulting log-slope is
    `−1/2 − (1/2)e^{−(T−t)/2}/(1 − e^{−(T−t)/2})`. That predicts **0.545** at `t = 15` and
    **0.502** at `t = 9`; the same run measures **0.548** and **0.507**. The bias, not mode
    contamination, is what dominates a late vorticity fit.

### Fixed

- **A rate fit could measure its own resolution floor.** `S − S_η` does not decay forever at a
  finite resolution: it settles on a floor set by how well the mesh represents the relaxed
  state. Measured on A4 at 20 spline cells, where it settles at `8.5e-8`, a late fit returned
  **0.457** for a rate that is exactly **1** — while the spectral run at the same time, with a
  lower floor, returned **1.022**. A floor expressed relative to the *maximum* cannot catch
  this, because the plateau's height depends on the resolution rather than on the initial
  excess. `settled_floor` scales it to the series' own terminal value instead, discarding only
  the last stretch of a well-resolved run and the whole plateau of a floor-limited one.

### Found

- **A4's initial condition as printed is discontinuous on the torus, by 8.5 % of its peak.**
  `eq:initial_gaussian` is written unmodified on `T²`, with no summation over periodic images.
  How much that matters depends on how far the centre sits from the wrap, and A4 is the one run
  where it is not negligible: its Gaussian is centred at `x₂ = 3π/2`, only `π/2` from the
  boundary, with `w₂ = 1` and amplitude 1.8, so it still has the value
  `1.8·exp(−(π/2)²) = 0.153` there — **8.5 %** of its own peak. A1 reaches only `1.2e-25` at
  the wrap and A2/A3 `5.2e-5`.

  This was found by the spline/spectral cross-check failing for A4 alone, at **6.9e-2** against
  the **2.5e-4** of A2, whose Gaussian is otherwise identical. The control settles the cause:
  periodising A4's Gaussian — the only change — brings the disagreement to **3.85e-4**, a
  **180×** improvement landing exactly in the band A1–A3 occupy. So the disagreement is the
  stated initial condition's own discontinuity and not a defect in either discretisation.
  Periodising changes nothing for A1 (**1.2e-25**) or A2/A3 (**5.2e-5**).

  Which reading the manuscript intends is not stated, so **both are run for A4 and both are
  reported**: the literal form as the primary, being what is printed, and the periodised form
  alongside it. `Gaussian` takes an `images` argument and `periodise` switches between them.

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
