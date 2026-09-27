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
