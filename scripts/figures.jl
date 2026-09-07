#!/usr/bin/env julia
#
# Regenerate every figure from the saved runs.
#
#     julia --project=scripts scripts/figures.jl [a1 … a4 b1 b2 b3] [--runs-dir DIR]
#                                                                   [--results-dir DIR]
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
                              initial_condition,
                              SECTION5_RUNS, SECTION5_ORDER, DIRICHLET_EIGENVALUE,
                              euler_entropy_floor
using Printf
using Serialization

include(joinpath(@__DIR__, "runner.jl"))
using .Runner

"""
    split_args(argv)

Separate the `--flag value` pairs, which belong to [`parse_options`](@ref), from the bare run
names, which are this script's own.

Splitting on `startswith(a, "--")` alone does not work: the *value* of a flag is not itself a
flag, so `--runs-dir /tmp/x` would leave `/tmp/x` looking like a run name. Every option this
script forwards takes exactly one value except `--quiet`, so the value is consumed with it.
"""
function split_args(argv)
    flags = String[]
    names = String[]
    i = 1
    while i <= length(argv)
        a = argv[i]
        if a == "--quiet"
            push!(flags, a)
        elseif startswith(a, "--")
            push!(flags, a)
            i < length(argv) && (push!(flags, argv[i += 1]))
        else
            push!(names, a)
        end
        i += 1
    end
    return flags, names
end

const flags, wanted = split_args(ARGS)
const opts = parse_options(flags)
const requested = isempty(wanted) ? [collect(SECTION4_ORDER); collect(SECTION5_ORDER)] :
                  wanted
const names = filter(in(SECTION4_ORDER), requested)

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
    # A1's payload carries no `scatter` field: its abscissa is the prescribed `h`, which is a
    # function of position and needs no run to produce, so the pair is assembled here from the
    # same `contours` grid the field maps use. The reduced Euler runs store theirs, because
    # there the abscissa is the state-dependent stream function.
    (sc_initial, sc_final) = if name == "a1"
        ((vec(contours), vec(reshape(trg.initial, N, N) .+ uΩ)),
            (vec(contours), vec(reshape(trg.final, N, N) .+ uΩ)))
    else
        (p.scatter.spectral_initial, p.scatter.spectral_final)
    end
    println(figure_scatter(out("scatter.png"), sc_initial, sc_final;
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
        println(figure_tau(out("tau.png"), p.contours, measured))
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

# =============================================================================================
# Section 5.4: the `sv_*`, `pe_*` and `ge_*` panels.
#
# Two figures per run rather than the four §4 gets, and the difference is not laziness. There is
# no second discretisation here, so there is nothing to overlay on the trace panel; there is no
# cone, because `eq:theoretical-limits` is a statement about the periodic torus; and the field
# maps would need the state resampled off the quadrature grid, which is a solver operation and
# does not belong in a figure script. What the manuscript's §5.4 panels actually show — the
# (φ, ω) cloud collapsing onto a curve, and the conservation and dissipation traces — is here.

for name in filter(in(SECTION5_ORDER), requested)
    path = joinpath(opts.runs_dir, name * ".jls")
    if !isfile(path)
        @printf("%-4s  not run yet (%s)\n", name, path)
        continue
    end
    p = open(deserialize, path)
    tr = last(p.traces)[2]
    out(f) = joinpath(opts.results_dir, name * "_" * f)

    # `S_η = λ₁,₁H₀` is the constrained minimum for `s = ω²/2` and has no closed form for
    # `s = ω log ω`, so B3's panel gets no reference line rather than a fabricated one.
    Sη = SECTION5_RUNS[name].entropy === :quadratic ?
         euler_entropy_floor(tr.H[1]; λ = p.λ.discrete) : nothing
    println(figure_traces(out("traces.png"), p.traces, Sη))

    # The reference over the cloud: the straight line ω = λφ for the quadratic entropy, the
    # curve ω = e^{λφ-1} for the Gibbs one. Both are the manuscript's own, drawn from the λ the
    # driver measured and printed.
    φv = p.scatter.final[1]
    xs = range(extrema(φv)...; length = 200)
    reference = if SECTION5_RUNS[name].entropy === :quadratic
        (collect(xs), p.λ.fitted .* xs)
    else
        (collect(xs), exp.(p.λ.formula .* xs .- 1))
    end
    println(figure_scatter(out("scatter.png"), p.scatter.initial, p.scatter.final;
        reference = reference, xlabel = "φ", ylabel = "ω"))
end
