#!/usr/bin/env julia
#
# Run every driver script.  Exits nonzero if any of them fails.
#
#     julia --project=scripts scripts/run_all.jl [script ...]
#
# With no arguments it runs the whole suite, in the order below, which follows the sections of
# *Metriplectic relaxation to equilibria*: the A runs are §4 on the periodic torus, the B runs
# §5.4 on the unit square, the C runs §5.5 on the Grad-Shafranov geometries.  Each C run depends
# on the machinery the B runs exercise, so the order is not arbitrary.
#
# Each script runs in its own process, so a failure -- or an `exit(1)` from `summary` -- is
# contained and reported rather than taking the runner down with it.
#
# Every driver takes `--runs-dir` and `--results-dir`; neither is passed here, so each writes to
# the repository's own `runs/` and `results/`, which are not tracked.

const SCRIPTS = String[
# §4  the periodic torus: analytical diffusion, then reduced Euler under the double and the
#     projector bracket
#   "run_a1.jl",
#   "run_a2.jl",
#   "run_a3.jl",
#   "run_a4.jl",
# §5.4  reduced Euler on [0,1]^2: single vortex, perturbed equilibrium, Gibbs entropy
#   "run_b1.jl",
#   "run_b2.jl",
#   "run_b3.jl",
# §5.5  Grad-Shafranov: the rectangle, then the mapped disk
#   "run_c1.jl",
#   "run_c2.jl",
# the convergence-rate fits over the A runs
#   "converge.jl"
]

const RULE = "="^67

function main(args)
    requested = isempty(args) ? SCRIPTS : args
    project = @__DIR__
    failed = String[]

    if isempty(requested)
        println("No driver scripts registered yet.")
        return 0
    end

    for s in requested
        path = joinpath(@__DIR__, s)
        isfile(path) ||
            (push!(failed, "$s (missing)"); @error "no such script" script=s; continue)
        println(RULE)
        println(s)
        println(RULE)
        flush(stdout)
        cmd = `$(Base.julia_cmd()) --project=$project --startup-file=no $path`
        success(pipeline(cmd; stdout, stderr)) || push!(failed, s)
    end

    println()
    if isempty(failed)
        println("All driver scripts passed.")
        return 0
    end
    println("$(length(failed)) SCRIPT(S) FAILED: " * join(failed, ", "))
    return 1
end

exit(main(ARGS))
