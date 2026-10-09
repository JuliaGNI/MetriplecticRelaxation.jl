# Known issues

Defects that a review found and that are not fixed yet. Each entry names its kind and its
evidence.

### K1 · README.md gives a command that does not run the suite

- location: `README.md:272`
- evidence: README.md:272 gives `julia --project=. --check-bounds=auto test/runtests.jl`. The
  test dependencies are in `test/Project.toml`, so with `--project=.` the command cannot load
  SafeTestsets: `Base.find_package("SafeTestsets")` returns `nothing`. `Pkg.test()` runs the
  suite.
- kind: docs
- found: 2026-09-27

### K3 · The test files copy imports, a seed and a comment they do not need

- location: `test/torus.jl`
- evidence: All 14 files moved from `runtests.jl` keep `using LinearAlgebra`, `using Random`,
  `using SparseArrays`, `Random.seed!(0x5c1e9a3b)` and the five-line "Every resolution below …"
  comment. Only `spline.jl`, `euler_square.jl`, `euler_references.jl`, `gibbs_entropy.jl` and
  `gradshafranov.jl` draw random numbers: `grep -c -E '\brandn?\(' test/*.jl` gives 0 for the
  other nine. Only `spline.jl` uses SparseArrays (`nnz`).
- kind: dead code
- found: 2026-09-27

### K4 · Five test files of `src/euler.jl` have no `test/euler.jl`

- location: `test/euler_problems.jl`
- evidence: `euler_problems.jl`, `euler_square.jl`, `free_state_space.jl`, `gibbs_entropy.jl`
  and `euler_references.jl` test mainly `src/euler.jl`. The D3 convention maps `src/<path>.jl`
  to `test/<path>.jl`. `test-layout.jl --check` accepts the current layout, because `src/` is
  flat.
- kind: defect
- found: 2026-09-27

### K5 · Two mutants survive the Section 5 tests

- location: `src/euler.jl:507`
- evidence: `src/euler.jl:507`, `euler_flow`: `mobility = (x, u) -> u` changed to
  `(x, u) -> -u` survives `test/collision_bracket.jl`. The sign cancels in the recentring.
  `src/diagnostics.jl:241`, `entropy_plateau`: `: 0.0` changed to `: 1.0` survives
  `test/euler_short_runs.jl`. No test asserts the `held` output of the `i == 1` branch. Both
  survive on `origin/main` too, because the test bodies are the same.
- kind: missing test
- found: 2026-09-27

### K6 · Comments and docs say GeometricBrackets is not registered

- location: `Project.toml:21`
- evidence: `Project.toml:21`, `test/Project.toml:23`, `scripts/Project.toml:16`,
  `README.md:281`, `CHANGELOG.md:51` and `CHANGELOG.md:63` (two `[Unreleased]` entries) say
  that GeometricBrackets is unregistered; `README.md:281` and `scripts/Project.toml:16` say the
  same of SimpleSplines, and `README.md:285` says that both resolve from their git remotes at
  `rev = "main"`; SimpleSplines has no `[sources]` entry. General holds GeometricBrackets 0.1.0
  to 0.2.0, and SimpleSplines 0.1.0 to 0.3.2.
- kind: docs
- found: 2026-10-08

### K7 · The `[sources]` entry of GeometricBrackets follows the remote `main`

- location: `Project.toml:40`
- evidence: `GeometricBrackets = {rev = "main", url = "https://github.com/JuliaGNI/GeometricBrackets.jl"}`
  at `Project.toml:40`, `test/Project.toml:26` and `scripts/Project.toml:21`. When that `main`
  moves to a version outside the `"0.2.0"` bound, no environment of this repository resolves.
  `ced4e9e` does not resolve: its `"0.1"` bound excludes the 0.2.0 that `main` carries.
- kind: defect
- found: 2026-10-08

### K8 · `Project.toml` and README call 1.10 the tree's usual Julia floor

- location: `Project.toml:36`
- evidence: `Project.toml:36-38` says the `julia` floor is "1.11 and not the tree's usual 1.10",
  and `README.md:283` says the same. GeometricBase 0.15 requires Julia 1.11, so every package
  that bounds `GeometricBase = "0.15"` has a Julia floor of 1.11 or higher.
- kind: docs
- found: 2026-10-08

### K9 · The floor `FFTW = "1"` cannot be installed on Julia 1.11

- location: `Project.toml:44`
- evidence: FFTW 1.0.x and 1.1.x depend on BinaryProvider, and every BinaryProvider version in
  General bounds `julia = ["0.7", "1.0-1.10"]`. On Julia 1.11.9, `resolve_versions!` with FFTW
  at 1.0.0 gives "Unsatisfiable requirements detected for package BinaryProvider". FFTW 1.2.x
  and 1.3.0 need AbstractFFTs 0.5, which the other bounds of `Project.toml` exclude, so the
  lowest FFTW that resolves here is 1.3.1. With FFTW at 1.3.1, every other direct dependency
  resolves at its floor. So the advisory CI `downgrade` job fails on FFTW before it reaches any
  other floor.
- kind: defect
- found: 2026-10-08
