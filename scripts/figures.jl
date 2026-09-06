#!/usr/bin/env julia
#
# Regenerate every Section 4 figure from the saved runs.
#
#     julia --project=scripts scripts/figures.jl [a1 a2 a3 a4] [--runs-dir DIR]
#                                                              [--results-dir DIR]
#
# Reads `<runs-dir>/<name>.jls`, which a `run_aN.jl` wrote, and draws into `<results-dir>`.
# Neither directory is tracked: figures are regenerated, never committed, and the numbers quoted
# in prose come from the drivers' own reports rather than from reading a plot.
#
# A run that has not been done yet is reported and skipped rather than being an error, so this
# is usable while the long ones are still going.

using MetriplecticRelaxation
using MetriplecticRelaxation: SECTION4_RUNS, SECTION4_ORDER, islands_h, relaxation_time,
                              euler_entropy_minimum, contour_average, CENTRAL_ISLANDS,
                              DOMAIN_LENGTH, figure_fields, figure_traces, figure_scatter,
                              figure_cone, figure_rates, figure_tau, fit_rate,
                              initial_condition
using Printf
using Serialization

include(joinpath(@__DIR__, "runner.jl"))
using .Runner

const args = filter(a -> startswith(a, "--"), ARGS)
const wanted = filter(a -> !startswith(a, "--"), ARGS)

# `parse_options` consumes only the flags; the bare run names are this script's own.
const opts = parse_options(vcat(args,
    reduce(vcat, [String[] for _ in 1:0]; init = String[])))
const names = isempty(wanted) ? collect(SECTION4_ORDER) : wanted

"The `N`-by-`N` grid a spectral run was sampled on, as two coordinate vectors."
function grid_axes(N)
    (collect(0:(N - 1)) .* (DOMAIN_LENGTH / N),
        collect(0:(N - 1)) .* (DOMAIN_LENGTH / N))
end

for name in names
    path = joinpath(opts.runs_dir, name * ".jls")
    if !isfile(path)
        @printf("%-4s  not run yet (%s)\n", name, path)
        continue
    end
    p = open(deserialize, path)
    spec = SECTION4_RUNS[name]
    N = p.opts.spectral
    x₁, x₂ = grid_axes(N)
    traces = p.traces
    trg = last(traces)[2]          # the spectral trace, for the field maps
    out(f) = joinpath(opts.results_dir, name * "_" * f)

    # ---------------------------------------------------------------------------------------
    # The field maps, drawn from the spectral run, with the contours of δH/δu over them.
    contours = if name == "a1"
        [islands_h(x₁[i], x₂[j]) for i in 1:N, j in 1:N]
    else
        # The stream function of the FINAL state: the manuscript overlays the contours the
        # relaxed vorticity is constant on, which is the visual form of the equilibrium claim.
        reshape(p.scatter.spectral_final[1], N, N)
    end
    uΩ = p.uΩ.spectral
    println(figure_fields(out("fields.png"), x₁, x₂,
        reshape(trg.initial, N, N) .+ uΩ, reshape(trg.final, N, N) .+ uΩ, contours;
        threshold = name == "a1" ? 1e-4 : nothing,
        label = name == "a1" ? "u" : "ω", colormap = name == "a1" ? :viridis : :balance))

    # ---------------------------------------------------------------------------------------
    Sη = name == "a1" ? trg.S[end] : euler_entropy_minimum(trg.H[1])
    println(figure_traces(out("traces.png"), traces, Sη))

    # ---------------------------------------------------------------------------------------
    # The scatter plot, with the closed-form reference the run is measured against.
    reference = if name == "a1"
        # Fig. 1's blue markers: the average of the initial condition on each contour of h,
        # over the island that holds the maximum.
        hs = collect(range(0.02, 0.98; length = 40))
        u₀ = initial_condition(spec)
        (hs, [contour_average(u₀, h, CENTRAL_ISLANDS[2]; n = 400) for h in hs])
    else
        nothing
    end
    println(figure_scatter(out("scatter.png"),
        p.scatter.spectral_initial, p.scatter.spectral_final;
        reference = reference, xlabel = name == "a1" ? "h" : "φ",
        ylabel = name == "a1" ? "u" : "ω"))

    # ---------------------------------------------------------------------------------------
    if name == "a1"
        # Fig. 1's lower panel, but as a test rather than an illustration: the closed form as a
        # curve, and the rate measured from the run as points over it.
        measured = Float64[]
        for h in p.contours
            τ = relaxation_time(h)
            w = (min(1.5τ / spec.T, 0.5), min(4.5τ / spec.T, 0.9))
            (rate, _, _) = fit_rate(p.deviation.t, p.deviation.y[h]; window = w,
                floor = 1e-8)
            push!(measured, 1 / rate)
        end
        println(figure_tau(out("tau.png"), p.contours,
            [relaxation_time(h) for h in p.contours], measured))
    else
        H₀ = trg.H[1]
        println(figure_cone(out("cone.png"), traces, H₀))
        if hasproperty(p, :distance)
            println(figure_rates(
                out("rates.png"), trg.t, trg.S .- euler_entropy_minimum(H₀),
                p.snapshots.spectral_t, p.distance.spectral))
        end
    end
end
