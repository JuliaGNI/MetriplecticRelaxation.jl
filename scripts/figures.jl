#!/usr/bin/env julia
#
# Regenerate every figure from the saved runs.
#
#     julia --project=scripts scripts/figures.jl [a1 … a4 b1 b2 b3 c1] [--runs-dir DIR]
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
                              SECTION5_RUNS, SECTION5_ORDER,
                              euler_entropy_floor, SECTION55_ORDER,
                              EulerSquare, euler_axis, euler_grid,
                              GradShafranovBox, gs_axes, gs_grid, gs_ordinate_grid
using PoissonBrackets: nbasis
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
const requested = if isempty(wanted)
    [collect(SECTION4_ORDER); collect(SECTION5_ORDER); collect(SECTION55_ORDER)]
else
    wanted
end
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
# Three figures per run rather than the four §4 gets, and the one that is missing is missing for
# a reason: there is no second discretisation here, so there is nothing to overlay on the trace
# panel, and there is no cone, because `eq:theoretical-limits` is a statement about the periodic
# torus.
#
# The field maps need the state resampled off the quadrature grid onto a uniform one, which is a
# solver operation rather than a plotting one — so it is `euler_grid` in `src/euler.jl`, checked
# in `verify_euler_grid.jl`, and only called from here. The space it resamples through is fixed
# by `(cells, degree, state)`: the payload records the first two, `SECTION5_RUNS` carries the
# third, so the square is REBUILT here rather than the grids being carried in the payload — one
# sparse Cholesky against re-running an hour of Newton solves per figure revision. `nbasis` is
# checked against the recorded `N`, which is what catches the third, so that a payload written
# against a different space stops the script instead of drawing a figure of the wrong state.

# Samples per axis for the §5.4 field maps. Odd on purpose: all three runs share an initial
# condition centred at (½,½), and an odd node count puts a sample exactly on the peak rather than
# straddling it — which matters because `figure_fields` shares one colour range between its two
# panels. 129 is 16 641 evaluations per field, seconds in all.
const SECTION5_SAMPLES = 129

for name in filter(in(SECTION5_ORDER), requested)
    path = joinpath(opts.runs_dir, name * ".jls")
    if !isfile(path)
        @printf("%-4s  not run yet (%s)\n", name, path)
        continue
    end
    p = open(deserialize, path)
    tr = last(p.traces)[2]
    out(f) = joinpath(opts.results_dir, name * "_" * f)

    # ---------------------------------------------------------------------------------------
    # The field maps: ω at t = 0 and t = T, with the contours of the FINAL stream function over
    # them — the same convention as §4's reduced Euler panels, and the visual form of the
    # equilibrium claim, since the relaxed vorticity is constant on those contours.
    state = SECTION5_RUNS[name].state
    sq = EulerSquare(p.opts.cells, p.opts.degree; state = state)
    nbasis(sq.space) == p.opts.N || error(
        "$(name).jls was written on a space of N = $(p.opts.N) degrees of freedom, but " *
        "EulerSquare($(p.opts.cells), $(p.opts.degree); state = :$(state)) has " *
        "N = $(nbasis(sq.space))")
    xs = euler_axis(SECTION5_SAMPLES)
    ω₀g = euler_grid(sq, tr.initial, SECTION5_SAMPLES)
    ωTg = euler_grid(sq, tr.final, SECTION5_SAMPLES)
    # §4's rule, taken from the data rather than hard-coded per run: `:balance` diverges about
    # its own midpoint, so it is only honest where the shared colour range of `figure_fields`
    # actually straddles zero. B2's added mode makes it straddle; B1's vortex and B3's floored
    # state are both single-signed, and drawing those with a diverging map would put white at
    # 0.47 where a reader takes it for zero. No threshold is needed to separate the cases:
    # measured over both panels, `min/max` is `-0.999` for B2 and exactly `0` for B1 and B3,
    # whose minimum is the Dirichlet edge itself.
    lo, hi = extrema(vcat(vec(ω₀g), vec(ωTg)))
    println(figure_fields(out("fields.png"), xs, xs, ω₀g, ωTg,
        euler_grid(sq, sq.Λ * tr.final, SECTION5_SAMPLES);
        label = "ω", colormap = lo < 0 && hi > 0 ? :balance : :viridis))

    # ---------------------------------------------------------------------------------------
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

# Samples per axis for the §5.5 field maps, radial then axial. Two counts and not one: this
# domain is [1,7]×[-9.5,9.5] and a single count would sample the short axis three times as finely
# as the long one. 61 and 191 make the two spacings EXACTLY equal — 6/60 and 19/190 are both
# 0.1 — and both are odd, which puts a sample on C1's Gaussian centre (4,0) rather than
# straddling it. Exact equality needs Nr = 6k+1 and Nz = 19k+1, and an odd Nz needs k even, so
# the pairs are (49, 153) at k = 8, (61, 191) at k = 10 and (73, 229) at k = 12.
# `verify_gradshafranov_grid.jl` asserts the odd counts and the centre; the matched spacing is a
# choice recorded here.
const SECTION55_SAMPLES = (61, 191)

# §5.5, the Grad-Shafranov runs. The same three panels as §5.4 and for the same reasons, with two
# differences that are the manuscript's own: the entropy floor is `λ_h H₀` for the Grad-Shafranov
# eigenvalue rather than the Dirichlet one, and the scatter ordinate is `u/(Cr²+D)` rather than
# the state — `eq:gs-ref` is a statement about that field, not about `u`.
for name in filter(in(SECTION55_ORDER), requested)
    path = joinpath(opts.runs_dir, name * ".jls")
    if !isfile(path)
        @printf("%-4s  not run yet (%s)\n", name, path)
        continue
    end
    p = open(deserialize, path)
    tr = last(p.traces)[2]
    out(f) = joinpath(opts.results_dir, name * "_" * f)

    # ---------------------------------------------------------------------------------------
    # The field maps, C1 only. C2 relaxes on a `GradShafranovDisk` and not on a box, so the
    # reconstruction below would build the wrong space for it — and its `gsc_*` maps are a
    # separate change that belongs with a C2 figure pass, not with C1's missing family.
    #
    # The colour field is `u/(Cr²+D)` and not the state, and the manuscript says so rather than
    # this reproduction choosing it: "Instead of plotting the state variable u directly, the
    # color plot represents the field u/(Cr²+D), which should be proportional to ψ if the system
    # reaches a state consistent with" `eq:gs-ref`. The figure IS that claim, and it only reads
    # as one in this variable — σ depends on r alone, so the contours of `u` are not the
    # contours of ψ even at the exact equilibrium. The overlay is therefore ψ = Λj, the flux
    # function, which is `eq:gs-ref`'s right-hand side and §5.4's `sq.Λ * tr.final` in the other
    # geometry.
    #
    # `N` is checked against the recorded value for the reason the §5.4 loop gives: a payload
    # written against a different space must stop the script rather than draw a figure of the
    # wrong state.
    if name == "c1"
        box = GradShafranovBox(p.spec.cells, p.spec.degree)
        nbasis(box.space) == p.opts.N || error(
            "$(name).jls was written on a space of N = $(p.opts.N) degrees of freedom, but " *
            "GradShafranovBox($(p.spec.cells), $(p.spec.degree)) has " *
            "N = $(nbasis(box.space))")
        (rs, zs) = gs_axes(SECTION55_SAMPLES...)
        y₀ = gs_ordinate_grid(box, tr.initial, SECTION55_SAMPLES...)
        yT = gs_ordinate_grid(box, tr.final, SECTION55_SAMPLES...)

        # The same question §5.4's loop asks — is a diverging map honest here? — and it needs a
        # TOLERANCE where §5.4 needs none. §5.4 can test `lo < 0 && hi > 0` exactly because its
        # minimum is the Dirichlet edge, which is exactly `0.0`. C1's is not: the interior
        # minimum is the undershoot of a projected Gaussian's tail, measured `lo = -3.34e-08`
        # against `hi = 1.05e-01`, so that test is TRUE and §5.4's rule as written picks
        # `:balance` and puts white at 0.053 — the very failure §5.4's own comment describes,
        # arriving through round-off rather than through a single-signed field.
        #
        # So the rule is the fraction of the range the weaker sign occupies, against 1 %, which
        # is about what a colour map resolves. No margin is being chosen here: C1 measures
        # `3.2e-07`, five orders below the threshold, and §5.4's B2 — the one run that genuinely
        # straddles — measures `0.999`, two orders above it. Nothing lies between.
        lo, hi = extrema(vcat(vec(y₀), vec(yT)))
        signed = min(-lo, hi) > 0.01 * max(-lo, hi)
        println(figure_fields(out("fields.png"), rs, zs, y₀, yT,
            gs_grid(box, box.Λ * tr.final, SECTION55_SAMPLES...);
            label = "u/(Cr²+D)", xlabel = "r", ylabel = "z",
            colormap = signed ? :balance : :viridis))
    end

    println(figure_traces(out("traces.png"), p.traces, p.λ.discrete * tr.H[1]))

    φv = p.scatter.final[1]
    xs = range(extrema(φv)...; length = 200)
    println(figure_scatter(out("scatter.png"), p.scatter.initial, p.scatter.final;
        reference = (collect(xs), p.λ.fitted .* xs), xlabel = "ψ", ylabel = "u/(Cr²+D)"))
end
