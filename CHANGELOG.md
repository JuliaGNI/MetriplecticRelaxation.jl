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

  **Superseded for A2–A4 further down this section: `T = 20`.** `T = 10` was measured to read
  0.597 for a rate that is exactly 1/2, and the entry that changed it records why. A1 keeps
  `T = 10`.

- **`scripts/verify_torus_geometry.jl`** — 26 checks on that geometry, run before any solver was
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

- **`scripts/verify_spectral.jl`** — 29 checks: the operator agreement above, the Poisson solve
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

- **`scripts/verify_spline.jl`** — 52 checks. All four brackets are symmetric, positive
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

  **Recorded choice — there are two rate-fit windows, `t ∈ [0.65T, 0.85T]` for the entropy and
  `t ∈ [0.35T, 0.55T]` for the vorticity.** The manuscript's own reading of Fig. 6 is that the
  rate is asymptotic: for the second initial condition "the initial entropy relaxation rate is
  slower, but approaches ≈ 1 as the trajectory approaches the vertex of the cone". Fitting from
  `t = 0` would average the transient with the asymptote and report neither.

  They differ because two errors pull in opposite directions as the window slides. Contaminating
  modes decay faster than the one being measured and inflate an *early* fit; using `ω(T)` as a
  stand-in for the limit inflates a *late* one, by the amount `reference_bias` predicts. The
  entropy fit is quadratic in the perturbation, so its contamination clears four times faster in
  the exponent and it can afford the later window; the vorticity fit cannot, and takes the
  earlier one, where the reference bias is still 0.002. `projector_run.jl` then asserts the late
  vorticity excess *is* the predicted bias, which is what identifies it as the `ω(T)` stand-in
  rather than a failure of the rate to be 1/2.

  `fit_rate` also returns `r²`, because a rate quoted without it says nothing about whether the
  decay was exponential at all, and drops samples below `1e-12` of the peak so that a round-off
  tail does not drag the slope toward zero.

- **`scripts/verify_diagnostics.jl`** — 31 checks, including the fit-vs-scan comparison above,
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

- **`src/euler.jl` — the §5.4 problem and its discretisation.** Reduced Euler on `[0,1]²` under
  the collision-like div–grad bracket of `eq:collision-bracket`, through PoissonBrackets'
  `CollisionBracket`, `MetriplecticFlow` and `ImplicitMidpoint` — which for this field **is**
  Crank–Nicolson, the manuscript's own method. `EulerSquare` is the space and its assembled
  operators, `EulerSpec` and `SECTION5_RUNS` the three runs, `GibbsEntropy` the `s = ω log ω`
  functional PoissonBrackets does not carry, and `eigenmode_fit`, `euler_entropy_floor`,
  `gibbs_lambda`, `gibbs_fit` and `gibbs_residual` the closed-form references.

  Unlike §4 there is one file rather than a problem definition plus two solvers: the
  collision bracket has no spectral counterpart in this tree, so a second file would separate
  nothing and the cross-discretisation control §4 has is simply not available here.

  **Recorded choice — `Λ` is an assembled dense matrix, where §4's `PoissonMap` is a cached
  factorisation.** Opposite choice, opposite reason. §4 applies `Λ` once per Runge–Kutta stage
  and never differentiates it, so a factorisation wins by the whole memory bandwidth. Here
  `CollisionBracket`'s analytic `metric_derivative` reads `∂φ̂/∂ω̂ₘ` as **column `m` of `Λ`**, so
  the columns have to exist. That also settles the caching question the plan raises: `K` is
  fixed in time and is factorised exactly **once**, in the constructor, before the time loop.

  **Recorded choice — `26²` cells of degree 2, not the manuscript's `64²` P2, and the reason is
  measurable.** The collision bracket is nonlocal, so its operator is dense and one Newton
  matrix is `N` dense assemblies of size `N²` — `O(N³)`. Measured on this machine at degree 2,
  one implicit-midpoint step takes **1.0 s** at `N = 256`, **7.7 s** at `N = 576` and **36 s**
  at `N = 1024`; the manuscript's `64²` is `N = 4096` and some 40 minutes *per step*. No §5.4
  claim needs that resolution — the relaxed states are the lowest Dirichlet eigenmode and a
  smooth exponential of it, and the discrete first eigenvalue is already within `3.4e-5` of
  `2π²` at **8** cells.

  What *does* need resolution is the **narrow direction of the initial condition**, `w₁ = 0.1`,
  and **26 is that threshold, measured**: it is the coarsest mesh on which the `L²` projection of
  the printed Gaussian stops oscillating below zero, `min ω₀` running `-1.4e-2`, `-5.8e-4`,
  `-4.2e-7`, `+6.0e-13` at 12, 16, 20 and 26 cells. A statement about the peak, not about the
  boundary: refining the *other* axis does not move it at all.

- **`scripts/verify_euler.jl`** — **101 checks** on all of that, run before any B run was
  believed, in six sections. `w²` and not `w` (below); the space (`|Ω| = 1`, `K` positive
  definite with `λ_min = 19.73925`, the discrete first eigenvalue converging **16×** on doubling
  from 8 to 16 cells, `φ = Λω` vanishing on `∂Ω` at **exactly** 0, `MΛ` symmetric at exactly 0);
  the bracket on that space for all three runs (degeneracy **4.2e-16 … 5.0e-16**, `metric_apply`
  against `metric_matrix` at **2.0e-13 … 4.2e-13**, the vector field orthogonal to `∂H/∂û` at
  **1.3e-16**, `(S,S) = −Ṡ` to 11 digits); both entropies' analytic gradients and Hessians
  against central differences (**1.4e-8** and **3.6e-10** worst); and the two references.

  **Three controls that must fail, and do.** A bracket generated by an unrelated prescribed
  field is still a clean metric bracket (**1.7e-16**) and breaks the *flow's* degeneracy
  (**1.7e-3**, a ratio of **9.8e12**) — the trap `MetriplecticFlow`'s own docstring warns about.
  An indefinite mobility loses positive semi-definiteness while keeping symmetry and degeneracy,
  which is why positivity is a separate assertion. And reading §5's printed `w²` as `w` puts the
  e-folding distance at 0.01 instead of 0.1, an error of order one.

  Two calibrations worth recording, because both first appeared as failures of a badly designed
  check rather than of the code. `degeneracy_residual` normalises by `max|𝔾| max|g|`, which
  bounds a *single* product rather than the `N`-term sum with its cancellation, so a complete
  misalignment still normalises to `1e-3`; the control therefore asserts the **ratio** to the
  clean value, not an absolute floor. And a mobility `x₁ − 0.5` integrates to zero over the
  square, so the bracket's moment `m₀ = ∫M dμ` vanishes and its recentring divides by it —
  which destroys the degeneracy for a reason unrelated to the sign. `x₁ − 0.3` is the control
  that isolates the sign.

- **`src/diagnostics.jl` now serves three discretisations rather than two.** `Diagnostics` takes
  its run specification as a type parameter so that an `EulerSpec` fits alongside a `RunSpec`,
  and `entropy` dispatches on the run's own entropy — the generic method is the quadratic one
  and would report a plausible number for B3 rather than an error, which is why B3 gets a method
  instead of a flag. Added `entropy_plateau`, which is what B2's claim needs and monotonicity
  cannot supply: the time at which `S` has fallen by 1 % of its total fall, and how flat it was
  before then.

### Changed

- **`scripts/check.jl` drops `failures` and `reset_failures!`.** The follow-up the export-list
  comment in that file deferred, now made in all three repositories that carry a copy of the
  harness — `PoissonBrackets/scripts/`, `Papers/Metriplectic Relaxation to Equilibria/scripts/`
  and here — so that the copies stay diffable rather than differing by two names.

  **Zero callers, re-measured, and by a second method.** The call-site table in the file
  reproduces row for row, `fmt`'s 94 in PoissonBrackets included. Beyond that count, none of the
  60 `using .Checks: …` lines that import the harness names either name — 19 here, 10 in the
  paper, 31 in PoissonBrackets, and none anywhere else. That is the position a `name(` count
  cannot see, and the one place an unused export could still have been in real use.

  The pair only ever made sense for an in-process runner that reads the tally between suites and
  then clears it, and `run_all.jl` is not that: it gives every script its own subprocess.
  `summary` stays the only read path, and `_failures` stays private with no accessor and a comment
  saying why, so it is not reinstated on the assumption that the state was left unreachable by
  oversight. **No number and no count moves**, and the failure path is checked directly rather
  than assumed: with a deliberate failing check, all three copies still record the label and still
  exit 1.

- **`SpectralTorus` no longer stores the bare wavenumber grids.** `k₁` and `k₂` were written by
  the constructor and read by nothing: every operator applies one of `ik₁`, `ik₂` or `negk²`,
  and `Δ⁻¹` is built from the locals. They are now locals. Measured,
  `Base.summarysize(SpectralTorus(256))` falls from **5 243 472 to 4 194 896 bytes** — the
  1.0 MiB two `256²` `Float64` grids cost — on every grid the manuscript's own resolution
  builds. **No number moves:** `verify_spectral.jl` is 29/29 against PoissonBrackets'
  independent differentiation matrix and the suite is 401/401, unchanged.

- **`README.md` gained the `§5.5` section it was missing.** The file documented `§4` and `§5.4`
  and stopped, though §5.5 is roughly half of #2 — `takeda.jl`, `gradshafranov.jl`, both their
  verification scripts and `run_c1.jl`. It now carries C1's results table, the two findings, and
  why C2's relaxation is deferred, in the shape the other two sections use. Its *Running* section
  also stops recommending `Pkg.test()` without qualification: that defaults to
  `--check-bounds=yes`, which inflates timings roughly fourfold and invalidates a precompile
  image shared with a run driver.

- **`run_b3.jl`'s sections are numbered consecutively.** They ran 0, 0b, 1, 2, 3, 5, 6 — there
  was no section 4, and the preamble's cross-reference pointed at the gap. Now 0, 0b, 1, 2, 3, 4,
  5.

- **Two comments say what to do rather than what once happened.** `verify_gradshafranov.jl`'s
  warning against perping `d` twice recorded that "it was made once and cost an 0.85 relative
  disagreement" — project history, which belongs in this file; it now states the mechanism and
  the size of the disagreement as a consequence instead. And the module file now records that
  `field` is imported from PoissonBrackets' internals: it is neither exported nor `public` there,
  nor in SimpleSplines, so a rename upstream breaks this package silently and the fix belongs
  upstream rather than here.

- **A1's spectral right-hand side hoists `X_h` out of the time loop.** `h` is prescribed and
  fixed there, but the nested-bracket form recomputes `∂₁h` and `∂₂h` on every evaluation —
  four of its eight derivative applications. `spectral_rhs` is now built once outside the loop,
  mirroring `spline_rhs`, rather than being rebuilt on every step. The two agree **exactly**
  (0.00e+00). An earlier incomplete hoist kept the computation out of the Runge-Kutta stages
  only, leaving it inside the time loop — that version did 18 derivative applications per step
  and ran in **41 minutes**. The final correct version hoists completely out of the time loop,
  achieving 16 applications per step, and A1's `256²` spectral run drops from **69 minutes to
  35 minutes**. That matters only for A1, whose `Δt = 1e-4` is set by the explicit stability
  limit and buys it ten times the step count of the other three.

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

- **`[sources]` names the GitHub remotes instead of relative sibling paths.** `PoissonBrackets`
  and `SimpleSplines` were `{path = "../../Packages/…"}`, which points outside the checkout and
  therefore does not exist on a runner: every CI job failed in `julia-actions/julia-buildpkg`
  with "expected package PoissonBrackets [2a9ffa85] to exist at path
  `/home/runner/work/Packages/PoissonBrackets`" before a single test ran. The whole matrix was
  gated on a directory layout only the author's machine has, and had been since the repository
  was created — `main` failed identically. Both are now `{rev = "main", url = …}`, in the root
  and in `scripts/`; `MetriplecticRelaxation = {path = ".."}` stays a path because it is inside
  the checkout.

  **This trades away what the path form was for.** Resolution now follows `main` on the remote,
  so a local edit to either sibling is invisible here until it is pushed. Develop against a
  working tree with `Pkg.develop(path=…)` in a scratch environment rather than by restoring the
  path here.

  **No measured number moves.** At the time of the switch both siblings' remote `main` matched
  the local working tree except for one commit touching `.githooks/pre-commit` in each, so the
  source resolved against is byte-identical to what produced the tables below. Verified by
  `git ls-remote` against the cached refs rather than assumed.

- **Dependencies updated to the latest available.** `GeometricBase` 0.14.10 → **0.14.11**, `JSON`
  1.7.1 → **1.8.0**, `Parsers` 2.8.8 → **3.0.0**; every other registered dependency was already
  at its registry maximum, and no `[compat]` bound is binding. The suite is unchanged at **133
  passed** and all **177** verification checks still pass, so none of these moved a result.

- **`spectral_step!` and `spline_step!` lost the `!`.** Neither mutates its argument — both
  docstrings said so outright while the name promised the opposite. Renamed to `spectral_step`
  and `spline_step` at all fourteen sites.

- **`SpectralTorus` precomputes its three spectral multipliers.** `∂₁`, `∂₂` and `laplacian`
  rebuilt `im .* k₁`, `im .* k₂` and `−(k₁² + k₂²)` on every call. These are per-grid constants,
  and at the manuscript's `256²` each is a megabyte; the operators are applied sixteen times per
  RK4 step, so this was **~17 MB of a step's allocations** — larger than the 12.06 MB the `X_h`
  hoist above recovers. They are now fields, built once in the constructor.

### Fixed

- **The review follow-ups of #1 and #2 were re-measured against the same dependency trees the
  §4 and §5 results were, and nothing moved.** `Manifest.toml` is gitignored and `[sources]`
  follows `main` on both remotes, so the resolved trees are recorded rather than committed;
  verified from the manifest's own `git-tree-sha1` before any of the work below, in both the
  root and the `scripts` project:

  - **PoissonBrackets**, tree `23dba270c9b1c10f9576a9a09043a0f4790d0cd8`;
  - **SimpleSplines**, tree `6c199ca136d88712bf7505129153538a324e5c03`.

  Both upstream tips have since moved — PoissonBrackets to `2fc0221`, SimpleSplines to
  `c8b6432`/v0.1.0 — and the manifests did not, so **the loaded code is still the measured
  code.** The three cheap verification scripts reproduce unchanged: `verify_euler.jl` 101/101,
  `verify_takeda.jl` 38/38, `verify_gradshafranov.jl` 42/42, all exit 0, with
  `λ_takeda = 0.030224858357` after 65 sweeps, `λ_h = 0.030234799207`, C2's
  `λ = 0.0025990851` and the continuum `0.0302346260` all to every digit already recorded.
  A1–A4 were re-run in full at the manuscript's own `256²` and reproduce every number in the
  *Results* sections below, again to every digit.

- **`run_a1.jl`'s contour-average row goes from `2e-3` to `5e-5`.** Re-measured over both
  central islands and all four contours: **9.21e-08, 1.26e-07, 2.46e-07, 3.91e-07, 1.51e-06,
  1.53e-06, 2.06e-06, 4.71e-06** of the initial peak. The old constant was **425×** the worst
  of them. `5e-5` is eleven times it. There is no law to set this against — the average is
  conserved exactly in the continuum and *not* by the semi-discrete flow, so the residual is
  truncation error — which is why the comment now says what the number is a multiple of.

  **Must-fail control:** the same run on 48 cells instead of 128, same degree 3, gives
  **9.92e-05, 1.32e-04, 1.79e-04, 2.29e-04, 2.74e-04, 4.05e-04, 5.90e-04** on seven of the
  eight rows — every one **passes** the old `2e-3` and **fails** the new `5e-5`. (The eighth,
  `h = 0.20` on the upper island, lands at 3.30e-05 and still passes; the row was tightened to
  the worst measurement, not to the best.)

- **`run_a2.jl`'s cross-discretisation row goes from `5e-2` to `5e-4`.** Re-measured,
  **6.978e-05** — the old constant was **716×** it. `5e-4` is seven times the measurement.
  It stays deliberately looser than A3's new `1e-6`: A2's double bracket relaxes
  *incompletely*, so its final state is not the one member of `eq:u-eta_Euler_periodic` that
  A3's conserved energy fixes, and the two discretisations are left disagreeing at their own
  truncation error rather than at their agreement on `H₀` — measured, A2's final state is
  **0.9132** of its own norm away from that family where A3's is 3.33e-05.

  **Must-fail control, isolating exactly the row that changed:** halving only the spline
  resolution — 32 cells instead of 64, same degree 3, same `64²` spectral reference — gives
  **4.300e-03**, which **passes** the old `5e-2` and **fails** the new `5e-4`, while the
  `H₀`/`S(0)`/`S(T)` rows beside it still pass at 6.80e-06, 3.39e-05 and 9.24e-05.

- **`projector_run.jl`'s three numeric bounds now bound something.** All three were constants
  read off one run, and two of them by three to five orders.

  The cross-discretisation rows switch on the run, for the reason `verify_spline.jl` already
  switches: A4 as printed is the only run whose two discretisations are handed *different*
  initial data, because `eq:initial_gaussian` is discontinuous on `T²` by 8.5 % of A4's own
  peak. Re-measured on this manifest, and identical to what #2 recorded:

  | run | final state | `H₀` | `S(0)` | `S(T)` |
  |:--|--:|--:|--:|--:|
  | `a3` | 1.419e-07 | 5.31e-08 | 2.78e-08 | 5.31e-08 |
  | `a4` as printed | 5.367e-04 | 1.29e-04 | 9.73e-05 | 1.29e-04 |
  | `a4-periodic` | 1.264e-07 | 5.94e-09 | 1.90e-09 | 5.91e-09 |

  So `1e-1` and `5e-2` become `5e-3`/`1e-3` for A4 as printed and `1e-6`/`5e-7` for the other
  two — 7.0× to 9.4× the worst measurement each has to admit, on quantities that are
  deterministic.

  **Measured, the old constants admitted a badly degraded discretisation and the new ones
  reject it.** `run_a3.jl --spectral 64 --cells 16 --degree 1` — two spline degrees and a
  factor 4 in cells below the run's own space — gives final state **1.284e-02**, `H₀`
  **2.57e-02**, `S(0)` **2.78e-02**, `S(T)` **1.31e-02**. Every one of those **passes** the old
  `1e-1`/`5e-2` and **fails** the new bounds.

  The vorticity-rate row's lower bound goes from `0.47` to **0.495**. Its own comment gives the
  reason as a finite-window curvature term "worth a few times 1e-3", which supports about
  0.495, not 0.47 — an order of magnitude more slack than the stated reason, so the two are
  reconciled in favour of the reason. Measured, the smallest rate in this window is A3's
  **0.50679** (both discretisations), A4's is 0.59607/0.59605 and A4-periodic's 0.59678, so
  0.495 clears the smallest by 2.3 %. Not fragile: even the degraded control run above still
  measures **0.50635**.

- **`verify_spline.jl`'s vector-field agreement tolerance switches on the run, as its
  initial-condition row above it already did.** One `1e-1`, set by A4, left A1–A3 asserted
  **20× to 280×** loose: measured, A1 gives `4.80e-03`, A2 `2.69e-03`, A3 `3.51e-04` and A4
  `8.54e-02`. A1–A3 now carry `1e-2`, twice A1's measured worst, and A4 keeps `1e-1` because
  its printed initial condition is discontinuous on `T²` by 8.5 % of its own peak — which
  section 7 of the same script isolates as the cause rather than assuming it.

  **The tightening is a bound and not a record, and that is measured.** Degrading only the
  spline space — degree 1 instead of 3 at the same 64 cells, or 16 cells instead of 64 at the
  same degree 3 — takes A3's row to `2.86e-02` and `4.33e-02`. Both **pass** the old `1e-1`
  and **fail** the new `1e-2`, so the row now rejects a discretisation the old one waved
  through. 52/52, unchanged, and no reproduced number moves.

- **`reset_failures!`'s docstring described a runner that does not exist.** It said "Only
  `run_all.jl` needs this, between scripts", but `run_all.jl` gives every script its own
  subprocess, so each starts with an empty tally and nothing is ever reset. **Zero callers in
  any of the three repositories that carry a copy of this harness** — here,
  `PoissonBrackets/scripts/` and `Papers/Metriplectic Relaxation to Equilibria/scripts/` — by
  grep in all three and by `trace_path(direction="inbound")` on the code graph in the two that
  are indexed. `failures` is dead in all three as well.

  Both are **left in place**, and the export list now says why: it is the shared harness's API
  rather than this repository's usage, so that a converted script moves between the copies
  unedited. The measured call-site counts are recorded there — `check_exact` (0 here, 53 in the
  paper, 33 in PoissonBrackets), `fmt` (0, 0, 94) and `normerr` (0, 0, 4) have no local caller
  and are **not** dead; dropping the two that are is one change across three repositories, not
  a change to the downstream copy.

- **Three of `verify_projector_factor.jl`'s rows could not fail, and are now reports plus one
  row that can.** The two "still dissipates entropy" rows asserted `dS/dt < 0` for `κ = 1` and
  `κ = 2`. With `H₀ = (φ,ω)/2` and `S₀ = (ω,ω)/2` the field gives
  `-dS/dt = 2S₀ - 2κH₀²/‖φ‖²`, and Cauchy-Schwarz gives `H₀²/‖φ‖² ≤ S₀/2`, so
  `-dS/dt ≥ (2-κ)S₀` — non-negative for every `κ ≤ 2`, and strict once section 1 has shown
  `φ ∦ ω`. **The bound is `(2-κ)S₀` and not `2S₀ - κS₀/2`, so it covers `κ ≤ 2` and not
  `κ ≤ 4`:** measured on this field, `κ = 3` and `κ = 4` give `-dS/dt = -5.067` and `-15.638`,
  i.e. entropy *production*. Both κ the manuscript's discrepancy is between are inside the
  bound, which is why the rows could not fail. Measured `H₀²/(‖φ‖²S₀) = 0.3967` against `1/2`.

  The third was `Π_H φ = 0`, whose comment claimed that *applying* the projector escaped the
  vacuity of writing the formula out at `v = φ`. It does not: `Π`'s body **is** that formula,
  so the coefficient is `x/x`, exactly `1.0`, and `φ .- 1.0 .* φ` is exactly zero elementwise.
  Measured, `all(iszero, Π(ψ))` holds for `ψ = φ`, for a random field and for `φ` scaled by
  `1e±9`.

  All three now print as report lines in the six-space form `verify_takeda.jl` and
  `verify_gradshafranov.jl` already use, and the closed form itself is asserted per `κ` in
  their place — **`rel 2.21e-16` and `8.07e-16`**. That row is the one with a must-fail
  control: on a field whose projector term carries a sign slip it fails at **rel 5.68e-01 and
  8.85e-01**, while the two rows it replaces still *pass* on that same broken field. Net
  14 checks to 13.

- **The `SimpleSplines` compat bound follows its 0.1.0 release. Compat only — no code changed
  and no number below moves.** SimpleSplines' `main` carried `version = "1.0.0-DEV"`, which the
  bound `"1"` satisfied; on 2026-09-07 it released `0.1.0`, a *downgrade* in numbering rather
  than the `1.0.0` the bound anticipated. Because `[sources]` resolves it from `rev = "main"`,
  the release propagated on the next resolve and every Julia job in #2's CI matrix then died in
  `buildpkg`, before a single test ran, with

  ```
  ERROR: LoadError: empty intersection between SimpleSplines@0.1.0 and project compatibility 1
  ```

  The bound is now `"0.1"`. **None of the recorded results is affected**, and that is measured
  rather than assumed: SimpleSplines' release commit `c8b6432` touches only its `Project.toml`
  version string and its own changelog, so the *source* tree that resolves now is identical to
  `d77a472` — and therefore to what every number in the *Results* sections below was produced
  against. The dependency-state section further down still records the tree-SHA that was
  actually resolved for those runs, and is deliberately left as it was.

  That release is a version on `main`, **not** a registration: SimpleSplines is still absent
  from General, with no tags. So this changes nothing about why `[sources]` is needed here and
  in `scripts/Project.toml`, nor about the `julia = "1.11"` floor those tables force — the
  `Project.toml` comment explaining both remains correct as written.

  This also needed the same one-line bound in **PoissonBrackets**
  ([`JuliaGNI/PoissonBrackets.jl#7`](https://github.com/JuliaGNI/PoissonBrackets.jl/pull/7)).
  A dependency's `[compat]` takes part in resolution, and
  PoissonBrackets is itself resolved here from `rev = "main"`, so its stale `"1"` kept
  resolution unsatisfiable no matter what this file said. That is the standing cost of tracking
  `main` on two unregistered remotes: an upstream release can turn this repository's CI red with
  no change here at all, as it did in a pull request that does not touch `Project.toml`.

- **Three numbers in the prose did not match what the code does, found in review of #2 and
  re-measured rather than adjusted.** All three were quoted in docstrings or comments and none
  was asserted anywhere, which is how they drifted:

  - **B3's printed initial condition fails on the *first* step, not after five.** The prose said
    "five implicit-midpoint steps at `Δt = 0.02` reach `-6.6e-5`". Re-measured: `min ω₀` starts
    at `+5.959889e-13` and step 1 raises `DomainError` with `-6.554556141995815e-5` — so
    `-6.6e-5` was right, but it is the argument the entropy raised on *inside* the first step and
    no step ever completes. The finding is strengthened, not weakened: the margin is too thin to
    survive one step.
  - **The 20-cell undershoot is `-4.2e-7`, not `-3.9e-7`** — a 9 % drift. Re-measured at degree
    2, `min ω₀` runs `-1.374746e-02`, `-5.830994e-04`, `-4.249944e-07`, `+5.959889e-13` at 12,
    16, 20 and 26 cells. `verify_euler.jl` asserted only the *sign* of each row while printing
    the value, so nothing caught it; it now asserts the three coarser magnitudes to
    `rtol = 1e-3` as well, which takes it from **98 checks to 101**. Measured deviations are
    `1.2e-07`, `1.7e-08` and `4.2e-08`, so the tolerance leaves four orders of platform margin.
    The 26-cell entry stays sign-only, its value being a cancellation residual.
  - **`run_b3.jl` quoted `~5e-12` and `-3.6e-5`** for the same two quantities, and wrote the
    first into the generated `results/b3.md`. Both corrected to the measured `6e-13` and
    `-6.6e-5`, which `README.md` and `B3_FLOOR` already had right.

- **B1's control on the entropy plateau accepted the diagnostic's own failure sentinel.**
  `entropy_plateau` returns `ib = 0` both when the total fall is non-positive and when no sample
  crosses the threshold — "the entropy never fell" — and `run_b1.jl` asserted `ib <= 2`, which
  `0` satisfies. A run that dissipated nothing would have passed the control that exists to show
  B1 dissipates immediately, i.e. the exact opposite of the claim. Now `1 <= ib <= 2`. B2's two
  uses of the same diagnostic were already safe: both bound quantities that are `NaN` in the
  sentinel cases.

- **`entropy_plateau`'s docstring prescribed a threshold the runs cannot meet.** It asked for
  `t_break` to be "a substantial fraction of `T`"; B2 measures `t_break ≈ 8.5` against `T = 150`,
  5.7 % of it, because `T` is set by how long the *decay* takes rather than the plateau. The
  docstring overpromised and `run_b2.jl`'s `t_break >= 4Δt` was the honest statement, so the
  docstring now says that and names the run that implements it.

- **C1's section 4 was labelled as a statement about the relaxation when it is one about the
  discretisations.** Its label read "the **relaxed** λ agrees with Takeda's iteration" while the
  asserted quantity is `gs_eigenvalue(box)`, the spline space's eigenvalue, which does not depend
  on the run — the whole block would pass if the relaxation had never been stepped. It is sound
  because section 3 separately pins `gs_rayleigh(box, tr.final)` to `λ_h` at `1e-10`; the labels
  and the CHANGELOG now say which check carries which claim.

- **A test asserted a floating-point tautology.** `all(isapprox.(u .* (1 ./ u), 1.0; rtol=1e-15))`
  is true of any nonzero float and so asserted nothing. Replaced with `all(>(0), u)`, the premise
  that actually makes the Hessian's weight `1/u` finite. The testset's real claim was already
  asserted on the line above, so no coverage was lost and the suite count is unchanged.

- **`README.md`'s dependency section contradicted `Project.toml`.** It said `PoissonBrackets` "is
  resolved from the local sibling checkout at `../../Packages/PoissonBrackets`"; both sources have
  been git urls at `rev = "main"` since `066e6cc`, because a relative path escapes the repository
  and fails in `buildpkg` on every CI runner. Pre-existing, but `README.md` is edited here.

- **A rate fit could measure its own resolution floor.** `S − S_η` does not decay forever at a
  finite resolution: it settles on a floor set by how well the mesh represents the relaxed
  state. Measured on A4 at 20 spline cells, where it settles at `8.5e-8`, a late fit returned
  **0.457** for a rate that is exactly **1** — while the spectral run at the same time, with a
  lower floor, returned **1.022**. A floor expressed relative to the *maximum* cannot catch
  this, because the plateau's height depends on the resolution rather than on the initial
  excess. `settled_floor` scales it to the series' own terminal value instead, discarding only
  the last stretch of a well-resolved run and the whole plateau of a floor-limited one.

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

- **A1's `X_h` was never hoisted out of the time loop, only out of the Runge-Kutta stages.**
  `spectral_step!` rebuilt the right-hand side on every step, recomputing `hamiltonian_field`
  10⁵ times over A1 — **12.06 MB** of the **110.6 MB** a step allocates, and 18 derivative
  applications per step against the 16 the docstring promised. `spectral_rhs` now mirrors
  `spline_rhs` and is built once outside the loop. Verified **bit-identical** (max|Δ| = 0 after
  20 steps at 256²), so the runs already completed remain valid; **10.9 %** fewer bytes per step
  and **14 %** less time per step.

  **Five checks asserted nothing.** The worst applied the projector formula with `v = φ`, i.e.
  `φ − (φ,φ)/(φ,φ)·φ`, which is identically zero for any field whatever; it now applies the
  projector and adds orthogonality (`Π_H ω ⊥ φ`), idempotency, and a control that `Π_H` is not
  the identity. Four rows in `scripts/verify_projector_rates.jl` had the literal condition `true`;
  they now assert that the effective fitted rate exceeds 1/2, falls as `T` grows, and tends to 1/2.
  This matters because the file quotes check counts as evidence.

  **The degree study in `converge.jl` was asserting the wrong claim.** It required the observed
  order to equal `p+1`, but the measured orders exceed it in every degree — **4.41, 4.67 and 6.89**
  for degrees 2, 3 and 4 over meshes of 32–96 cells — which is superconvergence on a uniform
  periodic mesh. `p+1` is a lower bound on the rate, not a prediction, so the study now asserts
  that the order is at least `p+1` and increases with the degree. (A4 as printed shows what a real
  failure looks like under the same measurement: order **1.10**.)

- **Two more checks that asserted a theorem rather than a claim.** The pattern the entry above
  describes had survived in two further places, both of them on a *headline* result:

  - `run_a2.jl` tested A2's central claim — that `S` plateaus above `S_η` — with `excess > 0`.
    But `S_η` is the *constrained minimum* of `S`, so `S ≥ S_η` holds for every admissible state
    by definition: the check could not fail, and passes for A3's **complete** relaxation too,
    whose excess is `2.35e-10`. It now asserts the size of the excess, `> 50 %` of `S_η`, which
    the measured **179 %** clears and a relaxed run misses by nine orders of magnitude.
  - `verify_projector_rates.jl` had a check labelled "`S − S_η` is quadratic in `δ`" whose
    condition was again `excess > 0` — true for a linear, cubic or constant dependence equally.
    It computed the ratio to `δ²` at two amplitudes and printed both without ever comparing
    them. It now asserts that the ratio is the same at `δ = 1e-3` and `δ = 1e-4` to under a per
    cent, which is what distinguishes quadratic from anything else.

- **A comment in `projector_run.jl` was refuted by the check three lines below it.** It claimed
  the fitted vorticity rate "must fall as the window moves later" and that "one rising with time
  would contradict the derivation" — while the very next check asserts a later window reads
  **higher** (0.507 → 0.548), which is the reference bias the entry above derives. Two effects
  push the rate above 1/2 and they move in opposite directions as the window slides; the comment
  now says so, and explains why the lower bound is 0.47 rather than 0.5.

- **`PoissonMap` factorised the bordered matrix twice.** `PoissonMap{T, typeof(lu(B)), …}(lu(B),
  …)` calls `lu` in the type parameter and again in the value — `typeof` is an ordinary call and
  evaluates its argument. Every spline torus paid for two LU decompositions of an `(N+1)`-square
  sparse matrix and kept one.

- **`parse_options` accepted an option as an option's value.** `--runs-dir --results-dir out`
  bound the flag *name* as a path and `mkpath` then created a directory called `--results-dir`;
  one had accumulated in the repository root. A missing or option-shaped value now raises.

- **Two inert lines, and one misleading pointer.** `run_a1.jl`'s contour-average summary had the
  literal condition `true` without the `[REPORTED]` marker the file's seven other such lines
  carry, and `converge.jl` multiplied a log-ratio by `log2(2)`, which is 1.

  `check.jl` pointed at a `torustools.jl` with no hint that it is in a *different* repository.
  It is real — `PoissonBrackets/scripts/torustools.jl` — and `check.jl` is a harness ported
  between three repositories here, so the reference was right and only the location was
  missing. It now says where the file is and that this repository has none. Worth recording
  because the first reading of it was that the file did not exist at all: a grep confined to
  this repository says exactly that, and would have justified deleting a correct comment.

- **Two of this file's own records had gone stale.** The "recorded choice" for the final time
  still read `T = 10` for every run, superseded further down by `T = 20` for A2–A4; and the
  rate-fit window was recorded as `[0.5T, 0.9T]` where the code uses `[0.65T, 0.85T]` for the
  entropy and `[0.35T, 0.55T]` for the vorticity. Both now state what the code does, and the
  window entry says why there are two.

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

- **A4's early entropy relaxation rate is *faster* than its late one, where the manuscript says
  slower.** Fig. 7's discussion reads: "The cosine terms shift the initial condition closer to
  the boundary. As a consequence the initial entropy relaxation rate is slower, but approaches
  ≈ 1 as the trajectory approaches the vertex of the cone." Measured over `t ∈ [0.2, 2]` against
  `t ∈ [13, 17]`:

  | | early | late |
  |:--|--:|--:|
  | spline | **1.24825** | 1.02011 |
  | spectral | **1.24820** | 1.02186 |

  The two discretisations agree to four digits, so this is a property of the equation. It is
  also what the linearised spectrum predicts: `cos(2x₂)` is a `λ = 4` mode, whose contribution
  to `S − S_η` decays at `2(1 − 1/4) = 3/2` against the `λ = 2` mode's `2(1 − 1/2) = 1`, so an
  initial condition loaded with `λ = 4` must shed entropy *faster* early and slow toward 1.

  **This is not necessarily an error in the manuscript.** Its reasoning is about the constant
  `κ_η = 2a` of `eq:shrunk-cone`, which vanishes as a state approaches the cone boundary — a
  statement about a *lower bound* on the rate, which an actual rate may freely exceed. What does
  not reproduce is the reading of it as a description of Fig. 7's observed rate. The
  reproduction therefore reports the early rate rather than asserting a direction, and
  `run_a4.jl` says so at the check.

- **§5.4's closed-form references are the `μ = 0`, `c = 0` case of a four-multiplier
  equilibrium, and it is the *state* space — not the geometry — that forces it.** The sharpest
  thing the §5 work turned up, and not a spline artefact: it is about which functions the space
  contains.

  The equilibrium condition of `eq:collision-bracket` is `∇(δS/δu) = λ ∇(δH/δu)`, hence
  `δS/δu = λφ + c·x + μ`, with `μ` and `c` the multipliers of the bracket's **other** Casimirs:
  the mass `∫u` and the momenta `∫x_k u`, each a Casimir because its variational derivative has
  vanishing gradient and `Q₂` annihilates it identically. **Both** of §5.4's references are the
  `μ = 0`, `c = 0` case — `ω = λ₁,₁φ` for `s = y²/2`, and `ω = e^{λφ−1}` with
  `λ = (M+S)/2H₀` for `s = y log y`, the second following from the first by one substitution
  into `S = ∫ω log ω`.

  In a homogeneous-Dirichlet space that case is **forced**, because `1`, `x₁` and `x₂` are not
  in the space: the discrete Casimirs are absent and nothing constrains the mass. In a space
  that contains the constants the mass *is* conserved, `μ` is fixed by it, and both references
  acquire an extra term. So all three runs use `V_D`, and the price is a mass that drifts —
  B1's by **+46.0 %**.

  **The `:free` alternative is implemented and B3 measures both**, because the claim deserves a
  control rather than a paragraph. `EulerSquare(n, p; state = :free)` puts `ω` in the plain
  clamped space and still solves for `φ` in `V_D ⊆ V` through the recombination matrix,
  `Λ = E (EᵀKE)⁻¹ EᵀM`. Measured at 16 cells over the same 25 steps, the two spaces separate
  exactly as the table predicts: `μ` shrinks through `-1.24, -0.85, -0.63, -0.47, -0.35` toward
  zero in `V_D` while `V` holds it at **−1.46**, and the residual against the manuscript's
  `μ = 0` reference falls to **0.147** in `V_D` and sticks at **0.41** in `V`. That is the mass
  Casimir being present or absent, seen directly.

  `verify_euler.jl` settles the algebra in both directions: `S = 2λH + (μ−1)M` closes to
  **1e-12** for manufactured states at `(λ,μ) = (3,0)`, `(7.5,0)` and `(3,0.4)`, and
  `gibbs_lambda` returns `λ` exactly at `μ = 0` and **15.374 instead of 3** at `μ = 0.4`.

- **B3's printed initial condition is not admissible for its own entropy, in any space.**
  `s = y log y` is undefined at `y ≤ 0` and `eq:M-condition` gives it the mobility `M = y`,
  which the same equation requires to be positive. The printed Gaussian is strictly positive as
  a *function* — but with `w₁² = 0.01` it decays to `1.4e-11` of its peak at the far corner of
  the square, and a Galerkin scheme cannot tell that from zero. Measured on the `26²` space:
  `min ω₀ = +6.0e-13` at the outermost quadrature node — the Gaussian's own value there, not a
  projection artefact — and **the first implicit-midpoint step at `Δt = 0.02` drives the iterate
  to `-6.6e-5` and raises a `DomainError`**, before any step completes. The admissible set has no
  interior around such a state, and no tolerance fixes that.

  **The continuum does not have this problem, and the reason says what the discretisation is
  losing:** `M = ω` vanishes where `ω` does, so the flux vanishes with it and the equation is
  degenerate-parabolic and positivity-preserving. A Galerkin projection of it is not.

  B3 therefore carries a positive background of **1 % of its peak** (`B3_FLOOR`), which is four
  orders above the undershoot it has to absorb and leaves `min ω₀ = 9.3e-3` — a real margin
  where the printed condition leaves `6e-13`. It costs 12 % of the mass, which §5.4 measures
  nothing against: the reference and its multiplier are both read off the *relaxed* state.
  `run_b3.jl` measures the printed condition failing before it uses the floored one, so the
  departure is justified by the run rather than asserted in a comment. The alternative — a
  positivity-preserving limiter, or evolving `log ω` — is a different scheme, and this
  reproduction does not write one.

### Results — A1, parallel diffusion under the metric double bracket (§4.1, Fig. 1)

`Δt = 1e-4` (the manuscript's), `T = 10`, spectral `256²` (the manuscript's), spline `128²`
cells of degree 3, 100 000 steps.

- **`τ_h = (ℓ_h/2π)²` is confirmed to three or four digits, by measurement rather than by
  restatement.** The manuscript computes `τ_h` from the formula and plots it; this fits the
  decay rate of the along-contour deviation of the run and compares it with `1/τ_h` from the
  closed form. On the upper central island:

  | h | ℓ_h | τ_h | measured rate | 1/τ_h | rel | r² |
  |--:|--:|--:|--:|--:|--:|--:|
  | 0.50 | 5.2441 | 0.6966 | **1.43563** | 1.43554 | 6e-5 | 1.000000 |
  | 0.40 | 6.1651 | 0.9628 | **1.03871** | 1.03868 | 3e-5 | 1.000000 |
  | 0.30 | 7.5782 | 1.4547 | **0.68762** | 0.68744 | 3e-4 | 1.000000 |
  | 0.20 | 10.0945 | 2.5811 | **0.38798** | 0.38742 | 1e-3 | 1.000000 |

  `τ_h` varies by a factor 3.7 across these contours, so the agreement is not a constant being
  matched by a constant.

- **`u` relaxes to the contour average `u_∞(h)`.** The average is conserved to
  **9.2e-8 … 4.7e-6** of the initial peak on both central islands, and the along-contour
  deviation collapses by factors of **120 to 298 000**, tracking `e^{−T/τ_h}` contour by
  contour.

- **The relaxation is incomplete, as §4.1 says.** The relaxed state is **5.95×** its own norm
  away from the fully relaxed `u_η` of `eq:u-eta_analytic` — it "retains some information of
  the initial condition", which is what motivates §4.2.

- Energy conserved at **9.79e-16** (spectral) and **2.12e-13** (spline) relative; entropy
  monotone throughout, `0.19578 → 0.04322` (spectral) and `→ 0.04265` (spline); vorticity mass
  at **5.7e-16**.

### Results — A2, reduced Euler under the metric double bracket (§4.1, Figs. 2–3)

`Δt = 1e-3` (the manuscript's), `T = 20`, spectral `256²` (the manuscript's), spline `64²`
cells of degree 3, 20 000 steps.

| claim | spline | spectral |
|:--|--:|--:|
| `H` conserved, max \|ΔH\|/\|H₀\| | 1.73e-13 | **1.66e-15** |
| `S(0)` | 0.2243696428 | 0.2243696491 |
| `S(T)` | 0.1863056311 | 0.1863057468 |
| `S_η = H₀` | 0.0667843242 | 0.0667843277 |
| `S(T) − S_η` | **+0.1195213** | **+0.1195214** |

- **The entropy plateaus strictly above `S_η`, by 179 % of `S_η` itself.** This is the
  manuscript's incomplete relaxation — "the entropy appears to converge to a value that is
  higher than its constrained minimum" — and it is the reason §4.2 introduces the projector
  bracket. It is confirmed here as the *expected* outcome: a run that drove `S` down to `S_η`
  would contradict the paper, not improve on it. A3 is the same initial condition, grid and
  time step under the projector bracket, and reaches `S_η` to **2.35e-10**; the only difference
  between the two runs is the bracket.

- **The double bracket completes 24.15 % of the entropy reduction available to it**, identically
  in both discretisations. The manuscript quantifies this nowhere — it says only that the limit
  is higher than `S_η` — so the fraction is a result of this reproduction rather than a
  reproduction of a published number. The plateau is genuine: `S` moves by 2.9e-3 of the
  remaining excess over the last 10 % of the run.

- Energy is conserved to machine precision, and entropy is monotone throughout with a worst
  increment of **−3.44e-05**, i.e. always decreasing. The two discretisations agree on the
  final state to **6.98e-05** and on `S(T)` to **6.2e-07**.

### Results — A3, reduced Euler under the projector bracket (§4.2, Figs. 4–6)

`Δt = 1e-3` (the manuscript's), `T = 20`, spectral `256²` (the manuscript's), spline `64²`
cells of degree 3, 20 000 steps. Uses the coefficient `2H/‖φ‖²`, not the `H/‖φ‖²` printed below
`eq:projector-brackets` — see the finding below.

| claim | spline | spectral |
|:--|--:|--:|
| `H` conserved, max \|ΔH\|/\|H₀\| | 9.08e-14 | **4.57e-15** |
| `S → S_η`, as `(S_T−S_η)/(S₀−S_η)` | 2.48e-10 | **2.35e-10** |
| `S − S_η` decay rate (**exact: 1**) | 0.99994 | **1.00047** |
| `‖ω(t)−ω(T)‖` decay rate (**exact: ½**) | 0.50679 | 0.50679 |
| `rS/rω` (**exact: 2**) | 1.973 | 1.974 |
| `ω(T)` on `𝔠_η`, `‖ω(T)−fit‖/‖ω(T)‖` | 3.33e-05 | 3.33e-05 |

Both rate fits have `r² = 0.9999970` or better. **The two discretisations agree on the final
state to 1.42e-07** and on `H₀`, `S(0)` and `S(T)` to `5.3e-08` — five orders of magnitude
closer than A1, because nothing here develops a separatrix layer.

The complete relaxation of §4.2 is confirmed against A2's incomplete relaxation under the same
initial condition, grid and time step: the only difference between the two runs is the bracket.
All three inequalities of `eq:theoretical-limits` hold along the whole trajectory, with worst
residuals of `−2.8e-10` to `−5.8e-10`, i.e. inside the cone throughout.

**Deviation from the manuscript — the vertex is reached at `t = 3.10`, not `t ≈ 5`.** Fig. 6
says "the system … reaches the vertex of the cone already at `t ≈ 5`". The number depends on
what "reaches" means, which the manuscript does not state; the criterion here is a **recorded
choice**: the first time the excess entropy falls below 1 % of its initial value, which is
where the trajectory becomes indistinguishable from the vertex at plot resolution. Both
discretisations give 3.10 independently, so it is not a discretisation artefact. A stricter
threshold would move it later, so this is consistent with the figure rather than in conflict
with it, but the number is not the manuscript's and is recorded as such.

### Results — A4, reduced Euler under the projector bracket, second IC (§4.2, Fig. 7)

`eq:ic-projector2`: `u₀ = cos(2x₂) + u_G`, `1/N = 1.8`, `x₀ = (π, 3π/2)`, `w = (0.3, 1)`.
Everything else is A3's. **Both readings of the ambiguous initial condition were run.**

| claim | literal, spectral | periodised, spectral |
|:--|--:|--:|
| `H₀` | 2.521919349 | 2.527758701 |
| `H` conserved, max \|ΔH\|/\|H₀\| | **7.92e-15** | **7.73e-15** |
| `S → S_η`, as `(S_T−S_η)/(S₀−S_η)` | **5.62e-11** | **5.53e-11** |
| `S − S_η` decay rate (**exact: 1**) | 1.02186 | 1.02213 |
| `‖ω(t)−ω(T)‖` decay rate (**exact: ½**) | 0.59605 | 0.59678 |

- **The two readings differ, but no claim does.** `H₀` differs by **2.18e-3** relative and
  `S(0)` by **1.79e-3** — small but not negligible, and a genuine difference in the physical
  setup rather than in its solution. Both give the entropy rate to within 0.03 of 1 and both
  reach `S_η` at `5.6e-11`, so the ambiguity changes the numbers without changing any
  conclusion. That is the useful outcome: the reproduction does not need the question settled.

- **The periodisation control is decisive at production resolution.** The spline/spectral
  agreement on `S(0)` goes from **9.73e-05** with the printed formula to **1.90e-09** when the
  Gaussian is periodised — **51 336× better**. Nothing else changes between the two runs, so
  the disagreement is the stated initial condition's own discontinuity.

- **`‖ω(t)−ω(T)‖` reads 0.596 rather than A3's 0.507**, and that is A4's initial condition
  rather than a defect: `eq:ic-projector2` contains `cos(2x₂)`, which **is** a `λ = 4`
  eigenmode of `−Δ`. The mode that contaminates a vorticity fit is therefore present at order
  one rather than as a tail, and the fit reads high by correspondingly more. The entropy fit is
  unaffected — 1.022, against A3's 1.000 — because in `‖ε‖²` the modes are separated by 0.5
  instead of 0.25.

- **`src/takeda.jl` — the §5.5 Grad–Shafranov problem and its classical solver.** The
  Herrnegger–Maschke profile constants, the measure `dμ = dr dz / r`, the mobility
  `M = Cr²+D` of `eq:M-condition`, the profile `f(r,y) = (Cr²+D)y`, the Galerkin matrix of
  `−Δ*` and of the profile, and the iteration of Takeda & Tokuda Eqs. (2.111)–(2.112) on a
  uniform finite-difference grid. **It carries no relaxation and depends on no bracket**, which
  is the point: the manuscript's `λ = 0.030302` is quoted as coming from "a standard
  Grad–Shafranov solver", so reproducing it from the metriplectic flow would settle nothing —
  the flow is what is being checked.

  **`f` is linear here, so the "eigenvalue problem" really is one.** `f(r,·)` is defined as the
  inverse of `∂_y s(r,·)`, and for `s = y²/2(Cr²+D)` that inverse is `(Cr²+D)y`. So
  `eq:Grad-Shafranov-equation` reads `−Δ*ψ = λ(Cr²+D)ψ`, a linear symmetric generalised
  eigenvalue problem, and Takeda's iteration — which is written for a general profile — becomes
  inverse power iteration whose fixed point is the smallest eigenvalue. That is what makes an
  algorithmically independent eigensolve available as a check on the iteration.

  **Recorded choice — the normalisation, which is what `λ` comes out of.** The classical scheme
  fixes the flux on the magnetic axis and reads `λ` off the amplitude the elliptic solve
  returns: `−Δ*ψ̃ = f(r,ψᵏ)`, then `λᵏ⁺¹ = ψ_axis / max ψ̃` and `ψᵏ⁺¹ = λᵏ⁺¹ψ̃`. `ψ_axis = 1`,
  stopping tolerance `1e-13` on the relative increment of `λ`, cap 2000 sweeps; measured, 65
  sweeps on the `64²` grid. The reference itself (pp. 22–23 of Takeda & Tokuda 1991) was not
  available to this session, so what is implemented is the scheme of that *class* as the
  manuscript describes it, and the docstring says so.

  **Recorded choice — the discretisation, and why it is not the runs'.** A second-order
  conservative five-point stencil for `div_μ ∇` with `1/r` at the **cell faces**, on the
  manuscript's own "uniform grid of `64 × 64` nodes", with the measure lumped. Face values
  rather than node values because they are shared between neighbours, which makes the matrix
  symmetric exactly and the Rayleigh quotient meaningful. It deliberately shares no assembly
  code with the spline space the relaxation runs on, so agreement between the two is evidence
  about the problem.

- **`scripts/verify_takeda.jl`** — 38 checks in seven sections, **ten of them controls that must
  fail and do**: the measure dropped from `Δ*`, from the profile matrix, from both, on each of
  the two geometries; the current profile substituted for the state profile in the eigensolver
  and again inside the iteration itself; and an unweighted stiffness matrix put through the
  `−Δ*` identity.

- **`src/gradshafranov.jl` — the §5.5 relaxation on the rectangle.** `GradShafranovBox` (the
  space, the `Δ*` solution operator, the energy and entropy forms), `gs_state`, `gs_flow`,
  `gs_fit`, `gs_rayleigh`, `gs_ordinate`, `gs_current`, and `SECTION55_RUNS` with C1's recorded
  choices. Plus `Diagnostics`, `potential`, `entropy` and `scatter_data` methods in
  `src/diagnostics.jl`, so the §5.5 run produces numbers comparable with §4's and §5.4's by
  construction.

  **The state variable is `j = u/r`, not `u`, and that is forced rather than presentational.**
  The manuscript's state is `u = (4π/c) r J_φ` with functional derivatives taken in the
  `L²(μ)` product. `CollisionBracket` forms `𝔾 = M⁻¹𝔸M⁻¹` with the space's own mass matrix, so
  energy conservation — `𝔾 ∂H/∂û = 0` — holds exactly when `M⁻¹∂H/∂û` is the bracket's
  generating field `ψ̂`. In the `u` formulation the elliptic solve is `Λ = (K^μ)⁻¹M^μ`, so
  `∂H/∂û = M^μψ̂` and the degeneracy would need the **μ-weighted** mass matrix in the sandwich,
  which no tensor-product mass operator provides. Substituting `u = rj` cancels the measure
  against the Jacobian: `∫ f δu dμ = ∫ f δj dx`, every functional derivative is numerically the
  same function in both formulations, and the weak `Δ*` solve becomes `Λ = (K^μ)⁻¹M` with the
  **plain** mass matrix on the right. Then `Λᵀ K^μ = M` and `∂H/∂ĵ = Mψ̂` exactly.

  Measured, on the same space and the same initial state: `𝔾 ∂H/∂ĵ` is **8.99e-16** in the `j`
  formulation and **3.74e-03** in the `u` one — a factor of 4·10¹². And `Λᵀ K^μ = M` holds to
  **1.93e-15**. The energy is still the poloidal magnetic energy: `½ĵᵀMΛĵ` equals
  `½ψ̂ᵀK^μψ̂ = ½∫|∇ψ|²dμ` to the last bit, and the entropy `½ĵᵀWĵ` equals
  `∫u²/2(Cr²+D) dμ` to the last bit, both checked rather than asserted.

  The price is recorded: `u_h = r j_h` is a spline times `r` and not itself a spline, so the
  trial space is `rV` where the manuscript's is `V`. Same order, different discretisation.

- **`scripts/verify_gradshafranov.jl`** — 42 checks in ten sections, **ten of them
  controls**, run before the relaxation. **`scripts/run_c1.jl`** and `run_gs` in `scripts/runner.jl`; a
  §5.5 block in `scripts/figures.jl`, whose scatter ordinate is `u/(Cr²+D)` and not the state,
  because `eq:gs-ref` is a statement about that field.

- **`test/runtests.jl` — a `Grad-Shafranov Problem Tests` set of 135 assertions**, coarse-grid
  versions of every §5.5 claim including four controls: the squared widths read as widths, the
  current profile read as the state profile, the `u` formulation's lost degeneracy, and the `dx`
  reading of `eq:entropy-2D`. The suite is 401 tests and green in 90 s.

- **`src/takeda.jl` also carries C2's geometry and its reference eigenvalue.** `disk_map`
  (`eq:mapping`), `DiskTriangulation`, `disk_matrices`, `disk_eigenvalue`, `disk_area`. C2's
  *relaxation* is deferred; its printed `λ = 0.002599` is not, and reproducing it needs no
  polar spline because a triangulation of the **physical** domain has the pole as an ordinary
  node.

### Results — the Grad–Shafranov reference eigenvalue λ (§5.5, step 11)

Dependency state as for §5.4 above, re-confirmed against `Manifest.toml`: **PoissonBrackets
`a13598b`**, tree `23dba270c9b1c10f9576a9a09043a0f4790d0cd8`; **SimpleSplines `c74e37d`**, tree
`6c199ca136d88712bf7505129153538a324e5c03`.

- **The manuscript's `λ = 0.030302` is a discrete number, and the continuum eigenvalue is
  `0.0302346260`.** Nine digits, three routes that share no assembly code. On the rectangle the
  problem *separates* — `Δ*` has no mixed term and `Cr²+D` depends on `r` alone — so
  `ψ = R(r)sin(π(z+9.5)/19)` reduces it to a one-dimensional Sturm–Liouville problem, which a
  1D spline Galerkin solve settles to a spread of **2.23e-11** over degrees 3–4 and 32–128
  cells. A **two-dimensional** Galerkin solve that uses no separation at all agrees to
  **4.9e-09**, which is what makes the separation legitimate rather than assumed. And the
  finite-volume iteration converges onto it at a **fitted order 2.06** over six grids,
  Richardson-extrapolating from 48²/96² to within **2.7e-06**.

- **The printed number is reproduced to 0.25 %, and that is the honest figure.** Takeda's
  iteration on the manuscript's own `64 × 64` node grid gives **0.0302248584**, i.e.
  `−2.55e-03` relative to `0.030302`. The last three digits of the printed value cannot be
  recovered, because recovering them means recovering the authors' element, mesh and stopping
  tolerance from one sentence of description.

  What *can* be said, and is measured: the printed value sits **0.22 % above** the continuum
  limit, and a conforming Galerkin eigenvalue approaches its limit **from above** — so the sign
  is the expected one. Tensor-product `Q₁` on this geometry brackets it: `0.0304253` at 16
  cells, `0.0303193` at 24, `0.0302823` at 32, `0.0302469` at 63, every one above the limit. So
  `0.030302` is what a low-order conforming element gives on a mesh of the stated order, and
  `0.0302346260` is what any of them converges to.

- **The measure is load-bearing in both places, and each is an O(1) error.** `dμ = dr dz/r` has
  to appear in `−Δ*` **and** in the profile matrix. Measured on the same space: the consistent
  pair gives `0.0302346`, no measure anywhere `0.0264256` (−13 %), the measure in `Δ*` only
  `0.0060680` (a factor 5 low), and in the profile matrix only `0.1199882` (a factor 4 high).
  None of the three is a tolerance-scale slip.

- **The manuscript prints two profiles, and mistaking them costs a factor 4.6.** §5.5 gives the
  equilibrium current as `(4π/c)J_φ = λ(Cr + D/r)ψ` and, equivalently, the state condition as
  `u/(Cr²+D) = λψ` — the same statement, because `u = r(4π/c)J_φ`. Reading the current form as
  `f` gives `λ = 0.1402132` against `0.0301943` on a `32²` grid, and it fails the defining
  identity `∂_y s(r, f(r,y)) = y` by **86 %**. Both are checked.

- **`ψ` is single-signed on the interior**, minimum `6.31e-04` against a peak of 1 at `r = 5`,
  which is what makes the scatter plot of `u/(Cr²+D)` against `ψ` a single line through the
  origin rather than two branches. The second axial mode is `0.0364702`, a ratio of **1.206**
  above the first, so the relaxed state's axial structure is unambiguous.

### Results — C1, Grad–Shafranov on the rectangle (§5.5, figs `gsr_*`)

**The recorded choices, all of them:** `Δt = 0.0625`, `T = 25` (400 steps), `18 × 21` cells of
degree 2 (`N = 378`, state space `V_D`, `j ∈ V_D`), quadrature at SimpleSplines' default
`⌈3p/2⌉ = 3` points per cell per axis, the nonlinear solve `Integrator`'s `NewtonSolver` with a
dense LU, run to `default_f_abstol = 3.35e-13` absolute rather than to an iteration count, and
the stopping criterion the excess entropy `(S − λ_h H₀)/λ_h H₀`. §5.5 states none of these. The
reasoning for each, with the measurements behind it, is in `SECTION55_RUNS`; what each produced
is below.

- **`H` is conserved at the Newton residual tolerance, and the bound has a derivation.** `H` is
  quadratic, so the midpoint increment satisfies `ΔH = ∂H/∂ĵ(j̄)·(ĵⁿ⁺¹−ĵⁿ)` exactly; the
  increment is `−Δt 𝔾(j̄)∂S/∂ĵ(j̄)` plus the Newton residual `ρ`, and the first term is
  annihilated by the degeneracy — so `ΔH = ∂H/∂ĵ(j̄)·ρ` per step, **independent of `Δt`**.
  Measured over 400 steps: `max|ΔH|/H₀ = 5.324e-12`, `max|ΔH| = 1.654e-12` on
  `H₀ = 3.106658399115e-01`, against `f_abstol/H₀ = 1.08e-12` — a factor **4.93** on a
  random-walk bound of `√400 = 20`.

- **`S` falls monotonically to `λ_h H₀` and then stops, and three diagnostics land on signed
  round-off because of it.** `S: 2.6564011390e-02 → 9.3929192900e-03`, and
  `(S(T) − λ_h H₀)/λ_h H₀ = −5.3172e-12` — it reaches the floor to twelve digits, and completes
  `1.000000000003` of the reduction available to it. Once there, the *evaluated* differences
  float: the worst entropy increment over the whole trace is `+7.836e-16` relative to `S₀`, and
  the entropy production reaches `−3.78e-18` against a maximum of `4.5348e-01` — a min/max of
  `−8.3e-18`. **Those are the equilibrium being reached, not a defect**, and `run_c1.jl` says so
  by splitting each claim: strict where the quantity is above round-off, sign-to-scale
  everywhere. The `Δt = 0.25` run that preceded this one did *not* reach that floor (its
  smallest production was `8.4e-12`), which is what makes the distinction worth drawing.

  For this entropy monotonicity is a **theorem at any step size**, unlike §5.4's B3: `S` is
  quadratic, so the same midpoint identity gives `S(ĵⁿ⁺¹) − S(ĵⁿ) = −Δt gᵀ𝔾g ≤ 0` by positive
  semi-definiteness alone.

- **The scatter collapses onto `u/(Cr²+D) = λψ`, and reaches the floor set by a projection.**
  `‖u/(Cr²+D) − λψ‖/‖u/(Cr²+D)‖` in `L²(μ)` falls from **6.3955e-01** at `t = 0` to
  **6.9980e-05**, which is *exactly* the number `verify_gradshafranov.jl` measures for the space
  itself (`6.998e-05` at 18×21). It is not zero and should not be: the discrete equilibrium
  satisfies `Π(σj) = λψ` — the `L²` **projection** of the ordinate is proportional to `ψ` — so
  the pointwise residual is the projection error of `σj`, and it falls with the mesh
  (`2.2e-03`, `4.3e-04`, `1.5e-04`, `7.0e-05` at 6×7, 10×12, 14×16, 18×21). The run is asked to
  reach that floor, and it reaches it to four digits.

- **Three estimates of `λ`, and they agree in the order their errors predict.** The Rayleigh
  quotient `S/H = 0.03023479920660` equals the space's own eigenvalue
  `λ_h = 0.03023479920660` to **1.95e-15** — that is the sharp statement that the run arrived,
  both being the smallest eigenvalue of the pencil `(MΛ, W)`. The pointwise fit gives
  `0.030234799155`, `−1.70e-09` relative, and the ratio `|Δλ|/(λ_h · rel²) = **0.3463**` — a
  clean second-order law, the §5.5 counterpart of §5.4's measured `0.2474`.

- **The relaxation and the classical solver agree — by way of the previous point, and not
  independently of it.** `λ_h` is the *spline space's* eigenvalue and does not depend on the run;
  what ties the run to it is `S/H = λ_h` to `1.95e-15` above. The three comparisons below are
  therefore statements about the two **discretisations**, and `run_c1.jl`'s section 4 would pass
  identically if the relaxation had never been stepped. `λ_h = 0.030234799207` against Takeda's
  iteration on the manuscript's own 64×64 node grid, `0.0302248584`: **+3.29e-04** relative,
  which is the finite-volume grid's own error and the larger of the two. Against the continuum
  eigenvalue, **+5.73e-06**. Against the manuscript's printed `0.030302`, **+2.22e-03** — and
  the printed value is `+2.23e-03` from the continuum one, so the reproduction and the paper
  disagree by essentially exactly the paper's own discretisation error, from the same side.

- **`Δt = 0.0625` is where the trajectory becomes resolved, and the coarse-mesh answer does not
  transfer.** Measured on the run's own space at `t = 2.5`: `‖Δt − Δt/2‖ = 2.738e-04`,
  `‖Δt/2 − Δt/4‖ = 7.390e-05`, **ratio 3.705** — implicit midpoint's second order. At
  `Δt = 0.25` the same difference is `1.35e-02` with a ratio of 9.8, i.e. above the asymptotic
  regime: Crank–Nicolson damps a stiff mode by only `4/(λΔt)` per step, so too large a step does
  not blow up, it slows the *apparent* relaxation as `Δt²`. **Nothing else in the run's output
  would have caught it** — `H`, monotonicity and the equilibrium are all `Δt`-independent — which
  is why the first pass of this run, at `Δt = 0.25`, passed every other check and was still
  wrong about the trajectory. The stiffest mode's rate grows like `h⁻²`, so the resolved step
  depends on the mesh and the `10 × 12` study that first suggested `Δt = 0.25` did not carry
  over to `18 × 21`.

- **The mass drifts by +83.7 %**, `∫u dμ = ∫j dx` from `1.0097733919` to `1.8545274685`. Not a
  defect and not asserted on: the constant function is not in `V_D`, which is the same absence
  that forces the equilibrium's multipliers `μ` and `c` to zero and so makes `eq:gs-ref` exact.
  §5.4's B1 drifts +46 % for the same reason.

- **The state-space control separates cleanly, and it has to be a relaxation.** Both spaces
  relaxed for 240 steps on `10 × 12` cells: the free space `V` misses `eq:gs-ref` by
  **2.5171e-01** against `V_D`'s **4.2761e-04** — a factor **589** — and carries the three
  multipliers `(μ, c₁, c₂)` at **2456×** the size, while its mass is held to **4.4e-16** where
  `V_D`'s drifts **+83.5 %**. That is S8's §5.4 finding reproduced on §5.5's problem.

  **The cheaper eigen-based control reports nothing, and that is worth recording.** The pencil
  `(MΛ, W)` carries no mass constraint, so its lowest eigenvector is the `μ = 0, c = 0` member in
  *both* spaces: measured, the two eigenvalues agree to **5.36e-09** and both eigenvectors fit
  `eq:gs-ref` to `4e-04`. What separates the spaces is that the **flow** in `V` conserves the
  mass and the momenta and so cannot reach that member from an initial state whose mass is
  nonzero. A control built on the fixed-point problem instead of on the flow would have passed
  while measuring the wrong thing — the first version of this check did exactly that.

- **Cost:** 400 steps in **1 090 s** at `N = 378` (2.7 s per step), plus 280 steps for the
  step-size study. The manuscript's `64 × 64` would be `N = 4096` and, at the measured `O(N³)`
  scaling of the nonlocal bracket's Jacobian, days.

### Results — C2's reference eigenvalue λ = 0.002599, on the mapped disk (§5.5, step 13, half)

**C2's relaxation is deferred; its reference eigenvalue is reproduced.** The two are separable
because the obstruction is about the *state space*, not about the geometry.

- **`eq:mapping` is transcribed correctly, and three consequences of the printed constants say
  so.** The image of the unit disk is `r ∈ [8, 16]` exactly, `z ∈ [−9.749139, 9.749139]`, of
  area **114.77699** (the Jacobian integrated over the parameter disk) — against C1's
  rectangle's 114, so the two experiments are on domains of deliberately comparable size — with
  the pole at `r = 11.412925`, just inboard of the Gaussian centre `r₀ = 12` the manuscript
  specifies. None of those four numbers is in the paper; all four are what its constants
  produce. The triangulated area approaches the exact one **from below** at second order —
  114.605, 114.734, 114.766 at 32×64, 64×128, 128×256 — because an inscribed polygon is smaller
  than the region it approximates, so a check on it names its mesh.

- **λ = 0.002599 is reproduced to 3.3e-05 relative — all six printed digits.** `P₁` finite
  elements on a triangulation of the *physical* domain, at 64 radial by 128 angular cells
  (8065 degrees of freedom), give **0.0025990851**. The continuum value is **0.0025970**,
  established by a **fitted order 1.9895** over seven refinements
  (0.0027256, 0.0026547, 0.0026296, 0.0026116, 0.0026052, 0.0026007, 0.0025991) and a Richardson
  extrapolation landing within **1.4e-05** of it. So the printed number sits **0.077 %** above
  the continuum limit — the same pattern as C1's `0.030302`, and from the same side, because a
  conforming Galerkin eigenvalue converges from above.

- **A control that behaves differently on the two geometries, reported rather than tuned.**
  Dropping `dμ` from `Δ*` alone moves λ to `0.00020973` (−92 %) and from the profile matrix
  alone to `0.03101267` (a factor 12); but dropping it from **both** moves λ only to
  `0.00255080`, **−2.3 %**. The weight sits on both sides of the Rayleigh quotient and largely
  cancels, and what survives is set by how far `r` varies: a factor of 2 on the disk against a
  factor of 7 on the rectangle, where the same control is −13 %. **C1 is the geometry that
  discriminates between the readings of the measure, and C2 is not** — a check run only on the
  disk would be far weaker than it looks.

### Results — C2's relaxation, deferred, and the two independent reasons

Not a partial result and not a regularised one. Two obstructions, either sufficient on its own:

1. **The pole of `eq:mapping`.** The isogeometric route — solve on the parameter square
   `(s,θ) ∈ [0,1]×[0,2π)`, periodic in `θ`, pulling the metric back through the map — is the
   right one and needs only a Jacobian-weighted assembly on top of `tensor_weighted_matrix`.
   But at `s = 0` the map collapses the **whole circle to a single point**: measured, the spread
   of `disk_map(0, θ)` over sixteen angles is `0.000e+00`. A function on the parameter square is
   therefore single-valued at the pole only if its `θ`-dependence there is constrained, and a
   tensor-product spline space provides no such constraint — so the space is not even `C⁰` at
   the pole, let alone `C¹`. Imposing it needs a **polar-spline** construction, the first two
   rows of the `s`-basis replaced by a three-function pole triangle, which is what
   Zoni & Güçlü build — the very paper the manuscript cites for the map — and which
   SimpleSplines does not have.
2. **PoissonBrackets has no triangular space.** Its only two-dimensional `DiscreteSpace` is the
   tensor-product `TensorSplineSpace`, and `CollisionBracket` reads `∇ψ` at *its* quadrature
   points. The `P₁` triangulation that gives the reference eigenvalue above cannot carry the
   flow, and the parameter square cannot carry the pole.

**Nothing was regularised to get around this.** No `ε` floor on `s`, no puncture at the origin,
no modified basis near `s = 0`. C2 joins the 3D Beltrami case (§6.4) as deferred, and for a
related reason: both need a spline construction the stack does not have — a polar space here, a
de Rham complex there.

### Results — the deliberate deviation converges to the paper's method

`scripts/converge.jl`, against a `192²` Fourier reference after 200 steps. Cubic B-splines are
4th-order accurate in `L²`, and that is what is observed:

| run | 16 | 24 | 32 | 48 | 64 cells | order |
|:--|--:|--:|--:|--:|--:|--:|
| a1 | 1.04e-1 | 1.73e-2 | 1.00e-2 | 1.31e-3 | 3.07e-4 | **4.12** |
| a2 | 3.35e-2 | 1.53e-2 | 5.26e-3 | 5.64e-4 | 1.45e-4 | **4.08** |
| a3 | 3.79e-2 | 1.27e-2 | 5.31e-3 | 6.07e-4 | 1.53e-4 | **4.05** |
| **a4, as printed** | 1.02e-2 | 4.56e-3 | 3.15e-3 | 2.39e-3 | 2.12e-3 | **1.10** |
| a4, periodised | 9.47e-3 | 3.20e-3 | 1.34e-3 | 1.53e-4 | 3.84e-5 | **4.05** |

So the spline runs converge *to* the manuscript's own discretisation, and the agreement between
them is a refinement statement rather than a coincidence at one resolution. **A4 as printed is
the sole exception, at order 1.10**, and periodising its Gaussian — the only change — restores
**4.05**. That is the discontinuity finding confirmed a third independent way, after the
initial-state comparison and the run itself.

**The order moves with the degree, which is what says it is the spline space's own and not an
artefact of the time stepping or of the reference.** Over meshes of 32–96 cells:

| degree | observed order | `p+1` | errors, 32 → 96 cells |
|--:|--:|--:|:--|
| 2 | 4.41 | 3 | 1.18e-2 → 8.83e-5 |
| 3 | 4.67 | 4 | 5.31e-3 → 3.16e-5 |
| 4 | 6.89 | 5 | 4.15e-3 → 2.07e-6 |

Every degree exceeds `p+1`, systematically — superconvergence on a uniform periodic mesh. Since
`p+1` is a lower bound on the rate rather than a prediction of it, converging faster than it is
not a defect, and the study asserts the two statements a broken space would violate while
superconvergence does not: the order is at least `p+1`, and it increases with the degree. A4 as
printed shows what a real failure looks like under the same measurement, at **1.10**.

### Results — where A1's two discretisations differ, and why

A1's spline and spectral final states differ by **7.99e-2** in `L²` where their initial states
differ by **2.21e-5** — a factor **3616**. `scripts/analyse_separatrix.jl` locates it:

| h | share of the squared difference | share of the nodes |
|:--|--:|--:|
| [0, 0.001) | **87.1 %** | 7.6 % |
| [0.001, 0.01) | 5.8 % | 11.6 % |
| [0.01, 1] | 7.1 % | 80.8 % |

**92.9 % of it lies within `h < 0.01`**, on 19.2 % of the nodes — the separatrix. The mechanism
is measured rather than asserted: `max|∇u|` near the separatrix relative to the island interiors
grows from **1.66** at `t = 0` to **90.74** at `t = T`, a 55× steepening.

This is a property of the equation, not of either method, and the manuscript states it in words:
`eq:parallel-diffusion` equalises `u` along the contours of `h` and moves nothing across them,
and at `h = 0` the contours are infinitely long, so `ℓ_h` and `τ_h` both diverge — "the dynamics
at the boundary of the islands is very slow ... the solution remains constant on those boundary
contours". Neither discretisation resolves an unboundedly steepening layer, and they fail to
resolve it differently. Every claim A1 makes — the contour averages, `τ_h`, the incomplete
relaxation — is confirmed independently by both runs.

Split by the global field norm, the total `7.985e-02` is **2.129e-02** from the 80.8 % of the
domain with `h ≥ 0.01` and **7.696e-02** from the 19.2 % nearest the separatrix. Away from the
layer the two discretisations agree at the level their resolutions predict; at it, neither is
converged and they are not converging to each other either.

`run_a1.jl` therefore asserts the comparison **away** from that layer and reports the total
beside it, both normalised by the global field norm. Normalising a masked region by the field
*inside* it is ill-posed here — A1's Gaussian sits on the separatrix, so `u ≈ 0` in the island
interiors, and the first version of the check reported `6.4e-2` for a region carrying 7 % of the
error.

### Results — §5.4, the dependency state every number below was measured against

`Manifest.toml` is gitignored and `[sources]` follows `main` on both remotes, so the resolved
tree is not committed and has to be recorded here instead. Every §5.4 number in the sections
that follow was produced against

- **PoissonBrackets `a13598b`**, tree `23dba270c9b1c10f9576a9a09043a0f4790d0cd8` — the merge of
  `metriplectic-brackets`, which is where `CollisionBracket`, `TensorSplineSpace` and
  `MetriplecticFlow` come from;
- **SimpleSplines `c74e37d`**, tree `6c199ca136d88712bf7505129153538a324e5c03`.

Both were confirmed against the manifest's `git-tree-sha1` before the runs, not assumed.

### Results — B1, the single vortex (§5.4, figs `sv_*`)

`Δt = 1`, `T = 200`, 200 steps, `26²` cells of degree 2 in the homogeneous-Dirichlet space,
`N = 676`, **57 minutes** of wall clock on an otherwise idle machine. Every one of `Δt`, `T`,
the mesh and the degree is a choice of this reproduction and not the manuscript's `64²` P2 —
the reasons are in `SECTION5_RUNS` and `EulerSquare` and are repeated in the *Added* entry
above.

- **`ω` relaxes to `ω = λ₁,₁φ`, and the check is two-sided.** The fitted `λ` reaches
  **19.7392146901** against the space's own eigenvalue `λ_h = 19.7392146644` — a relative
  agreement of **1.30e-09** — and `2π² = 19.7392088022`, from which both differ by **5.9e-06**,
  which is the discretisation error of the 26-cell space and not of the run. The residual
  `‖ω − λφ‖/‖ω‖` falls from **6.267e-01** at `t = 0` to **7.254e-05** at `t = T`.

  Both halves are needed. `λ` alone is a projection coefficient, `(ω,φ)/(φ,φ)`, and its error is
  *second* order in the state's: at `t = 0`, when the state is nothing like the eigenmode, it
  already reads **22.584**, within 14 % of the answer. A check on `λ` alone would pass for a run
  that had barely moved. B2's section quantifies that second-order law.

- **`S` decreases monotonically to `S_η = λ₁,₁H₀`.** Strictly: the worst increment over the 200
  samples is **−3.310e-10** of `S₀`, i.e. every step decreased it. `S` goes
  **2.0775801677e-02 → 1.1027476429e-02** against `λ_h H₀ = 1.1027476356e-02`, an excess of
  **+6.56e-09** relative, and **99.999999 %** of the available reduction. The Poincaré floor
  `S ≥ λ₁,₁H₀` holds at every sample with a minimum margin of **+3.35e-09**.

- **`H` is conserved at the Newton residual tolerance, which is the right way to say it.**
  `max |ΔH|/|H₀| = 2.449e-12` with `H₀ = 5.586583126006e-04` — an **absolute** drift of
  **1.4e-15**. That is machine precision on a quantity of this size, and the relative number is
  larger only because `default_f_abstol` is absolute: `4 max(8,N) eps ‖ω̂₀‖_∞ = 6.0e-13` at
  `N = 676`, and it knows nothing about `H₀`. The first version of this check asserted
  `< 1e-12` on the relative error and failed at `2.449e-12`; the threshold now follows the
  measurement and the check reports both numbers.

- **The mass drifts by +46.0 %**, from `0.0824356619` to `0.1203719455`, and that is not a
  defect. `∫ω` is a Casimir of the continuous bracket, and the discretisation inherits it only
  if the constant function is in the space — which in `V_D` it is not. That same absence is what
  makes `ω = λ₁,₁φ` the exact relaxed state; see the *Found* entry above. The driver reports the
  drift and does not assert on it.

- **The entropy production stays strictly positive and spans nine orders**, from
  `(S,S) = 2.879e-03` at the start to `6.57e-12` at `t = 200`, which is the run running out of
  entropy to dissipate rather than the bracket losing definiteness.

- **`Δt = 1` is a cost choice and this is the measurement that says so.** Five steps at `Δt = 1`
  and ten at `Δt = 0.5` land on states differing by **9.53e-04** relative in `L²`. One step is
  one Newton solve and one Newton matrix is `N` dense assemblies, so at `N = 676` a step costs
  **27 s**; a step size chosen for accuracy rather than for cost would have bought nothing here.

- **B1 is also the control for B2's plateau.** `entropy_plateau` on this trace puts the 1 %
  break at the **first** sample, `t = 1`, with nothing held before it — a run that starts far
  from any equilibrium dissipates from the first step.

### Results — B2, the perturbed equilibrium (§5.4, figs `pe_*`)

`Δt = 0.5`, `T = 150`, 300 steps, same `26²` degree-2 space, `N = 676`, just under two hours on a
machine also running the other two. Half B1's step and three quarters of its horizon, and both
halves of that are the plateau's doing — see `SECTION5_RUNS`.

- **The entropy is flat first and then monotone, which is the two-phase claim, and both halves
  are asserted.** `S` holds within **8.10e-03** of `S₀` — measured against its own total fall —
  for the first **8.5** time units, 18 samples of 301, and then falls. Over the whole run it
  never rises: worst increment **−5.751e-07**, `1.2500144159e-01 → 4.8108181391e-03`, which is
  **99.997 %** of the reduction available to it.

  Monotonicity alone would not say this. **B1, on the same mesh and with the same diagnostic,
  breaks at the first sample with nothing held** — that is the control, and it is what makes
  "flat first" a measurement rather than a description.

- **The mechanism is measured, not asserted: the entropy production grows by a factor 1070**,
  from `(S,S) = 1.214e-05` at `t = 0` to `1.302e-02` at its peak at `t = 18.5`. A plateau
  produced by a bracket that simply dissipates slowly would show no such growth.

- **The initial state really is near an equilibrium, and how near is a property of the mesh.**
  `‖ω₀ − λφ₀‖/‖ω₀‖ = 6.89e-02` with `λ₀ = 510.84` against the exact `52π² = 513.219`. The
  Gaussian is a 1 % perturbation by construction; the rest is the projection error of a mode
  with three wavelengths across 26 cells, and it is the larger of the two. So what seeds the
  instability here is partly the mesh, which is reported rather than hidden — the plateau's
  *existence* does not depend on which, only its length does.

- **It relaxes to the lowest eigenmode, not the one it started on**: `λ` falls from **510.84**
  to **19.74213425** against `λ_h = 19.73921466`, and the state becomes single-signed
  (`ω(T) ∈ [3.26e-05, 0.205]`) where `ω₀` was not (`[−0.999, 1.000]`). `H` is conserved at
  **3.21e-14** relative, **7.81e-18** absolute.

- **The tolerance on the fitted `λ` is `rel²`, and that is a measured law rather than a
  loosened threshold.** `λ` is a projection coefficient, so its error is second order in the
  state's. Measured on the two runs, three orders of magnitude apart in `rel`:

  | run | `‖ω−λφ‖/‖ω‖` | `|Δλ|/λ_h` | ratio |
  |:--|--:|--:|--:|
  | B1 | 7.2545e-05 | 1.3020e-09 | **0.24739** |
  | B2 | 2.4451e-02 | 1.4791e-04 | **0.24740** |

  Five digits of agreement. The first version of B2's check asserted an absolute `1e-4` and
  failed at `1.479e-04`; the check now asserts `|Δλ|/λ_h < rel²`, which both runs satisfy with a
  factor of four to spare and which tightens on its own as a run relaxes.

- **The mass drifts from `8.24e-04` to `7.84e-02`** — the same absence of the constants as in
  B1, and reported the same way.

### Results — B3, the Gibbs entropy (§5.4, figs `ge_*`)

`Δt = 0.02`, `T = 3`, 150 steps, same `26²` degree-2 homogeneous-Dirichlet space, `N = 676`,
plus the step-size sweep, the two-space control and the step-halving check on top. `Δt` is the
one number in §5.4 that is *measured* rather than chosen, and B3's initial condition carries
`B3_FLOOR` — see the two *Found* entries above for why.

- **`Δt` is bounded by admissibility, not by accuracy, and the sweep locates the boundary.** Ten
  steps at each of seven step sizes on the run's own mesh:

  | `Δt` | worst `ΔS/S₀` | `max \|ΔH\|/\|H₀\|` | outcome |
  |--:|--:|--:|:--|
  | 0.005 | −9.179e-03 | 8.149e-15 | monotone |
  | 0.01 | −8.006e-03 | 8.715e-15 | monotone |
  | 0.02 | −7.623e-03 | 8.715e-15 | monotone |
  | 0.04 | — | — | **left the admissible set** |
  | 0.08 | — | — | **left the admissible set** |
  | 0.16 | — | — | **left the admissible set** |
  | 0.5 | — | — | **left the admissible set** |

  The outcome is monotone in `Δt` — three `ok` then four failures, no interleaving — so there is
  a single threshold rather than a scatter, and it lies between **0.02 and 0.04**. That is the
  quantitative content of the manuscript's "sufficiently small time steps must be used", and it
  **tightens with the mesh**: `0.08→0.16` at 12 cells, `0.05→0.1` at 16, `0.02→0.04` at 26. What
  fails is not monotonicity but *admissibility* — `ω` goes negative and `y log y` raises.

  The rows are printed as measurements rather than as checks. The first version asserted each
  one and reported **four failures** for a sweep doing exactly what it exists to do; the claims
  are now the three conclusions drawn from the table.

- **It relaxes to the manuscript's `ω = e^{λφ−1}`.** With `λ = (M+S)/2H₀ = 15.03442967` from the
  run's own final mass, entropy and initial energy, the residual `‖ω − e^{λφ−1}‖/‖ω‖` falls from
  **4.556e-01** at `t = 0` to **4.238e-03** at `t = T`, measured over the interior at a margin of
  two cells.

- **The remaining disagreement is a one-cell boundary layer, and that is shown rather than
  assumed.** `e^{λφ−1}` is `e^{-1} = 0.3679` on `∂Ω` while every state of `V_D` is zero there, so
  no relaxed state in that space can match it in the last cell. The residual at margins
  `0, h, 2h, 4h` is **3.565e-02, 9.333e-03, 4.238e-03, 9.121e-04** — monotone in the margin,
  which is what identifies it as a layer.

- **`μ` relaxes to zero, which is what makes the manuscript's `λ` formula an identity here.** The
  independent two-parameter fit returns `λ_fit = 14.97322000` and **`μ = +0.00710146`**, i.e.
  `μ/λ = 4.7e-04`, and the formula agrees with the fit to **4.09e-03**. The three-quantity
  identity `S = 2λH₀ + (μ−1)M` closes to **2.05e-03** against a fit residual of `1.96e-02` — the
  accuracy the state itself permits, which is what the check is tolerance to.

  The fitted pair does **no better** than `μ = 0` (`3.588e-02, 9.764e-03, 4.545e-03, 2.154e-03`
  at the four margins, against `μ = 0`'s numbers above). That is the sharpest single statement
  that `μ` has genuinely vanished rather than been assumed away.

- **The two-space control, and it separates cleanly.** 25 steps at `Δt = 0.02` on a `16²` mesh in
  each space:

  | space | `∫ω` | fitted `λ` | fitted `μ` | `(M+S)/2H₀` | `‖ω − e^{λφ−1}‖/‖ω‖` |
  |:--|:--|--:|--:|--:|--:|
  | `:dirichlet` | 0.92028530 → 1.15511611 | 17.54191 | **−0.35055** | 14.53294 | **1.473e-01** |
  | `:free` | 0.92493151 → 0.92493151 | 25.12436 | **−1.45159** | 14.16493 | **3.913e-01** |

  The mass is held to every digit in `V` and drifts 26 % in `V_D`, `μ` is four times smaller
  where it drifts, and the manuscript's reference fits 2.7× better there. Over the full run in
  `V_D`, `μ` continues to `+0.0071`.

- **`H` is conserved at 1.392e-14 relative, 8.535e-16 absolute**, and `S` falls monotonically
  from `1.0112986621` to `6.1014492682e-01` with a worst increment of **−8.559e-11** — every step
  decreased it. `ω` stays admissible throughout, its minimum over the run being its initial
  `9.324e-03`, i.e. the floor's margin is never eroded. The entropy production runs from
  `12.854` to `4.04e-09`. The mass drifts **+33.75 %**.

- **`Δt = 0.02` is resolved**: ten steps against twenty at `Δt/2` land on states differing by
  **2.583e-04** relative in `L²`, so within the admissible range the step is not limiting the
  accuracy either.
