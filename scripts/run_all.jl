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
    # The verification scripts first: each settles a property the runs below rest on, and each is
    # seconds rather than minutes. A failure here invalidates everything after it, so running them
    # first is what makes a failed sweep readable.
    "verify_torus_geometry.jl",
    "verify_spectral.jl",
    "verify_spline.jl",
    "verify_diagnostics.jl",
    "verify_projector_factor.jl",
    "verify_projector_rates.jl",
    "verify_euler.jl",
    # §5.5's reference eigenvalue, before anything that compares a run against it: the whole
    # point of `takeda.jl` is that λ = 0.030302 comes from a computation the relaxation has no
    # part in, so its verification runs before the relaxation is built.
    "verify_takeda.jl",
    # The refinement study behind the deliberate deviation: the spline discretisation must converge
    # to the manuscript's spectral one, at the order the spline space has.
    "converge.jl",
    # §4  the periodic torus: analytical diffusion, then reduced Euler under the double and the
    #     projector bracket. A1 is the long one -- the manuscript's Δt = 1e-4 buys it ten times the
    #     step count of the others.
    "run_a1.jl",
    "run_a2.jl",
    "run_a3.jl",
    "run_a4.jl",
    # §5.4  reduced Euler on [0,1]^2 with homogeneous Dirichlet conditions: single vortex,
    #       perturbed equilibrium, Gibbs entropy. Each step is a Newton solve on a NONLOCAL
    #       bracket whose Jacobian is N dense assemblies -- 27 s per step at 26^2 cells -- so
    #       these are by far the most expensive scripts here: an hour or more each, and B3 is a
    #       step-size study and a two-space control on top of its own relaxation.
    "run_b1.jl",
    "run_b2.jl",
    "run_b3.jl",
    # §5.5  Grad-Shafranov on the rectangle. Cheaper than the B runs despite the larger
    #       problem, because its Δt is set by accuracy rather than by cost: 100 steps at
    #       Δt = 0.25 against B1's 200 at Δt = 1, and ~5 s per step at 18×21 cells.
    #
    #       C2, the mapped disk, has no driver: its relaxation is deferred, and the reason is in
    #       `disk_eigenvalue` and in `CHANGELOG.md`. Its reference eigenvalue is not deferred and
    #       is checked in `verify_takeda.jl`.
    "run_c1.jl",
    # The figures, last, because they read what the runs wrote.
    "figures.jl"
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
