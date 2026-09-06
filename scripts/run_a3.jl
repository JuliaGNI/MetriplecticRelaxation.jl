#!/usr/bin/env julia
#
# Run A3: reduced Euler under the projector bracket, Section 4.2, Figs. 4-6.
#
#     julia --project=scripts scripts/run_a3.jl [--runs-dir DIR] [--results-dir DIR]
#                                               [--spectral N] [--cells N] [--degree P]
#
# Same initial condition, grid and time step as A2 -- the manuscript is explicit that "the
# initial condition as well as the numerical method and the numerical parameters (grid size and
# time steps) are the same as in Fig. 2" -- but the metric bracket is `eq:projector-brackets`
# rather than the double bracket. That single change is what turns A2's incomplete relaxation
# into complete relaxation, so the two runs are a matched pair and are meant to be read against
# each other.
#
# Every check lives in `projector_run.jl`, shared with A4.

using MetriplecticRelaxation
using MetriplecticRelaxation: SECTION4_RUNS
using Printf

# Only this include: `projector_run.jl` brings in `check.jl` and `runner.jl` itself and
# re-exports them. Including either of those again here would make a second `Options` type.
include(joinpath(@__DIR__, "projector_run.jl"))
using .ProjectorRun
using .ProjectorRun: check_summary

const spec = SECTION4_RUNS["a3"]
const opts = parse_options()

println("A3  --  reduced Euler, projector bracket  (Section 4.2, Figs. 4-6)")
@printf("    Δt = %.0e   T = %.1f   steps = %d   spectral %d²   spline %d² cells, degree %d\n",
    spec.Δt, spec.T, round(Int, spec.T / spec.Δt), opts.spectral, opts.cells, opts.degree)
flush(stdout)

section("running")
const res = projector_run(spec, opts)

projector_checks(res, spec, opts)

section("writing")
let (p, c) = save_run(opts, "a3", projector_payload(res, spec, opts, "a3"))
    println("    run     -> ", p)
    println("    trace   -> ", c)
    println("    report  -> ",
        report(opts, "a3", projector_report(res, spec, opts, "a3")))
end

check_summary("run_a3.jl")
