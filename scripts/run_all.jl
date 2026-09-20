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
    "verify_euler_grid.jl",
    # The control behind §4.2's 𝔠_η membership bound: an under-relaxed state must be rejected by
    # it. Minutes rather than seconds — it steps A3 to half its final time — and the only one
    # here that is a statement about a TOLERANCE rather than about the mathematics.
    "verify_projector_tolerance.jl",
    # §5.5's reference eigenvalue, before anything that compares a run against it: the whole
    # point of `takeda.jl` is that λ = 0.030302 comes from a computation the relaxation has no
    # part in, so its verification runs before the relaxation is built. `verify_gradshafranov.jl`
    # then settles the measure and the state variable, and it relaxes two spaces as a control, so
    # it is minutes rather than seconds.
    "verify_takeda.jl",
    "verify_gradshafranov.jl",
    # And the resampling §5.5's field maps are drawn on, after the space it resamples from is
    # settled and before `figures.jl` draws with it. Seconds, like its §5.4 counterpart.
    "verify_gradshafranov_grid.jl",
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
    # §5.5  Grad-Shafranov on the rectangle. Comparable in cost to the B runs: 400 steps at
    #       Δt = 0.0625 against B1's 200 at Δt = 1, and ~2.7 s per step at 18×21 cells, plus a
    #       280-step study of the step size itself -- some twenty minutes in all. Δt is small
    #       because it is set by ACCURACY here, unlike §5.4's, where it was set by cost.
    #
    #       C2, the mapped disk, is the expensive one: 375 steps at Δt = 0.004 on 12×24 cubic
    #       cells, ~22 s per step in the transient and far less at the fixed point, plus a
    #       4375-step study of the step size on a coarser mesh at ~1.3 s per step. Well over an
    #       hour in all, and the mesh is set by that cost -- 16×32 measures 674 s per step.
    #       It reproduces `eq:gs-ref`; the state space is what decides that, and its own header
    #       says why.
    "run_c1.jl",
    "run_c2.jl",
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

# `exit` is what this file is for, so the rule does not apply. `run_all.jl` is only ever a
# process: it is run as `julia --project=scripts scripts/run_all.jl` and its exit status is the
# verdict on the sweep. `main` returns 1 when a driver failed, and returning that number to an
# interactive caller instead would make a failed sweep exit 0.
# fatou-ignore discouraged-function
exit(main(ARGS))
