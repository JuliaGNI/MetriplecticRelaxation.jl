# Known issues

Defects that a review found and that are not fixed yet. Each entry names its kind and its
evidence.

## KI-1 · docs · README.md gives a command that no longer runs the suite

README.md:272 gives `julia --project=. --check-bounds=auto test/runtests.jl`. The test
dependencies are in `test/Project.toml`, so with `--project=.` the command cannot load
SafeTestsets: `Base.find_package("SafeTestsets")` returns `nothing`. `Pkg.test()` runs the suite.

## KI-2 · docs · the CHANGELOG names the wrong slowest test file

The newest `[Unreleased]` entry says that the slowest file after compilation is
`gradshafranov.jl`, at 6.6 s. A second-run measurement in one process puts `quality/aqua.jl`
first, at 11.1 s to 13.7 s, against 6.6 s to 7.1 s for `gradshafranov.jl`. All files stay under
the 60 s budget.

## KI-3 · dead code · the test files copy imports, a seed and comments they do not need

All 14 files moved from `runtests.jl` keep `using LinearAlgebra`, `using Random`,
`using SparseArrays`, `Random.seed!(0x5c1e9a3b)`, the seed comment and the five-line
"Every resolution below …" comment. Only `spline.jl`, `euler_square.jl`, `euler_references.jl`,
`gibbs_entropy.jl` and `gradshafranov.jl` draw random numbers: `grep -c -E '\brandn?\(' test/*.jl`
gives 0 for the other nine. Only `spline.jl` uses SparseArrays (`nnz`). In the nine other files,
the comment "The runs here excite fields with random degrees of freedom" is false.

## KI-4 · layout · five test files of `src/euler.jl` have no `test/euler.jl`

`euler_problems.jl`, `euler_square.jl`, `free_state_space.jl`, `gibbs_entropy.jl` and
`euler_references.jl` test mainly `src/euler.jl`. The D3 convention maps `src/<path>.jl` to
`test/<path>.jl`. `test-layout.jl --check` accepts the current layout, because `src/` is flat.

## KI-5 · missing test · two mutants survive the Section 5 tests

- `src/euler.jl`, `euler_flow`: `mobility = (x, u) -> u` changed to `(x, u) -> -u` survives
  `test/collision_bracket.jl`. The sign cancels in the recentring.
- `src/diagnostics.jl:241`, `entropy_plateau`: `: 0.0` changed to `: 1.0` survives
  `test/euler_short_runs.jl`. No test asserts the `held` output of the `i == 1` branch.

Both survive on `origin/main` too, because the test bodies are the same.
