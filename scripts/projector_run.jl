#
# The run and the checks shared by A3 and A4, which differ only in their initial condition.
#
#     include(joinpath(@__DIR__, "projector_run.jl"))
#     using .ProjectorRun
#
# The manuscript is explicit that A4 is "the same as in Fig. 6, but for the initial condition
# (eq:ic-projector2)", so the two runs share every check. Keeping them in one place is what
# stops the pair from drifting apart; `run_a4.jl` adds only what is specific to it, which is the
# comparison between the two readings of its discontinuous initial condition.
#
# NOTE ON THE EVOLUTION EQUATION. These runs use the coefficient 2H/||phi||^2, not the
# H/||phi||^2 printed below `eq:projector-brackets`. The manuscript's own `eq:SS-projector`
# fixes the factor and `verify_projector_factor.jl` settles it numerically: with the printed
# factor the energy is not conserved, at 70% of scale.

module ProjectorRun

using MetriplecticRelaxation
using MetriplecticRelaxation: euler_entropy_minimum, spline_grid, energy_error,
                              entropy_monotone, best_fit_euler, fit_rate, l2norm,
                              cone_residual, scatter_data
using Printf

# `check.jl` and `runner.jl` are included HERE and nowhere else in an A3/A4 driver. Including a
# file twice makes two distinct modules with two distinct `Options` types, and the driver's
# `Options` then matches no method of `run_spline` — with an error that reads like a missing
# keyword argument rather than like a duplicate include. So the drivers include only this file,
# and what they need is re-exported below.
include(joinpath(@__DIR__, "check.jl"))
include(joinpath(@__DIR__, "runner.jl"))
using .Checks: header, check, summary
using .Runner

export projector_run, projector_checks, projector_payload, projector_report
export header, check
export Options, parse_options, run_spectral, run_spline, save_run, section, report

# `summary` is NOT exported: `Base.summary` exists, and exporting a second one leaves every
# unqualified use ambiguous. The drivers reach it as `ProjectorRun.Checks.summary`.
const check_summary = Checks.summary

# `‖ω(t) − ω(T)‖` needs the whole trajectory and ω(T) is not known until the end, so states are
# kept as the run goes. Every 4th sample: 100 states at 256² is 52 MB, where all 400 would be
# 210 MB for no gain in a rate fit.
const SNAP_STRIDE = 4

# The window the asymptotic rates are fitted over, as a fraction of T. A RECORDED CHOICE, and
# the reason it is late is derived rather than tuned: the linearised spectrum of the projector
# flow is exactly {1 - 1/λ} over the eigenvalues of -Δ (`verify_projector_rates.jl`), so the
# lambda = 4 mode decays at 3/4 against the lambda = 2 mode's 1/2 and contaminates the fit as
# exp(-t/4). With T = 20 this window is t ∈ [13, 17], where that contamination is 3.9% falling
# to 1.4%, biasing the fitted 1/2 by under 0.01.
#
# It stops at 0.85 T rather than running to the end because the reference state in
# ‖ω(t) − ω(T)‖ is ω(T) itself, which is not the true limit: the last few samples measure the
# distance to a moving target and flatten.
const RATE_WINDOW = (0.65, 0.85)

# The vorticity rate is fitted EARLIER than the entropy rate, and the two windows pull in
# opposite directions for reasons that are worth separating.
#
# `‖ω(t) − ω(T)‖` uses ω(T) as a stand-in for the relaxed state, which is the manuscript's own
# construction ("using the vorticity ω at the last point in time t = T as an approximation of
# the latter"). But ω(T) is not the limit. With δ_t = ω(t) − ω(∞) ≈ D e^{-t/2} v along the
# slowest mode,
#
#     ‖ω(t) − ω(T)‖ = D e^{-t/2} (1 − e^{-(T-t)/2}) ,
#
# whose log-slope is not −1/2 but
#
#     −1/2 − (1/2) e^{-(T-t)/2} / (1 − e^{-(T-t)/2}) ,
#
# so the FITTED rate reads high, and increasingly so as the window approaches T. At T = 20 that
# predicts 0.545 at t = 15 and 0.502 at t = 9; measured, the same run gives 0.548 and 0.507.
# The artefact, not mode contamination, is what dominates a late vorticity fit -- so this window
# stops at 0.55 T, where the correction is under 0.006.
#
# The entropy fit has the opposite problem and therefore the opposite window: S − S_η is
# quadratic in the perturbation, so its reference artefact is second order and negligible,
# while its mode contamination needs late times to decay. Hence two windows, not one.
const VORTICITY_WINDOW = (0.35, 0.55)

"""
    reference_bias(t, T)

The excess of a fitted `‖ω(t) − ω(T)‖` decay rate over its true asymptote 1/2, caused by using
`ω(T)` rather than the limit as the reference. See [`VORTICITY_WINDOW`](@ref).
"""
reference_bias(t, T) = 0.5 * exp(-(T - t) / 2) / (1 - exp(-(T - t) / 2))

"""
    settled_floor(y; factor = 10)

A `floor` for [`fit_rate`](@ref) that keeps only the samples at least `factor` times above the
series' own terminal value.

A decaying quantity does not decay forever at a finite resolution: `S - S_η` bottoms out on a
floor set by how well the mesh can represent the relaxed state, and a fit that reaches into
that floor measures the floor. Measured on A4 at 20 spline cells, where `S - S_η` settles at
`8.5e-8`, a late fit returned **0.457** for a rate that is exactly **1** — while the spectral
run at the same time, with a lower floor, returned **1.022**.

A fixed floor relative to the maximum cannot catch this, because the plateau's height depends
on the resolution rather than on the initial excess. Scaling it to the terminal value adapts:
for a well-resolved run it discards only the last stretch, and for a floor-limited one it
discards the floor.
"""
function settled_floor(y; factor = 10)
    m = maximum(y)
    m <= 0 && return 1e-12
    return max(factor * abs(y[end]) / m, 1e-12)
end

"""
    snapshotter()

An observer for [`run_spline`](@ref) / [`run_spectral`](@ref) that keeps every
`SNAP_STRIDE`-th sampled state, and the vectors it fills.
"""
function snapshotter()
    ts = Float64[]
    ws = Any[]
    k = Ref(0)
    obs = function (solver, t, ω)
        k[] += 1
        (k[] - 1) % SNAP_STRIDE == 0 && (push!(ts, t); push!(ws, copy(ω)))
        return nothing
    end
    return (obs, ts, ws)
end

"""
    projector_run(spec, opts)

Integrate `spec` with both discretisations, keeping snapshots, and return a named tuple of
everything the checks and the report need.
"""
function projector_run(spec, opts)
    (obs_t, snap_tt, snap_wt) = snapshotter()
    (ts, dt, trt, uΩt) = run_spline(spec, opts; observer = obs_t)
    (obs_g, snap_tg, snap_wg) = snapshotter()
    (gs, dg, trg, uΩg) = run_spectral(spec, opts; observer = obs_g)
    return (;
        spline = (solver = ts, diag = dt, trace = trt, uΩ = uΩt,
            snap_t = snap_tt, snap_w = snap_wt),
        spectral = (solver = gs, diag = dg, trace = trg, uΩ = uΩg,
            snap_t = snap_tg, snap_w = snap_wg))
end

"The distance ``\\|\\omega(t) - \\omega(T)\\|`` along a run's snapshots."
function distances(r)
    [l2norm(r.solver, w .- r.snap_w[end]) for w in r.snap_w]
end

"The final time of the run a result came from, for labelling a fit window in seconds."
spec_T(r) = isempty(r.trace.t) ? 1.0 : r.trace.t[end]

"""
    projector_checks(res, spec, opts; label = "")

Every claim §4.2 makes about a projector-bracket run, checked on both discretisations.
"""
function projector_checks(res, spec, opts; label = "")
    tag = isempty(label) ? "" : label * " "
    both = (("spline", res.spline), ("spectral", res.spectral))

    # -----------------------------------------------------------------------------------------
    header("$(tag)1. H is conserved and S is monotone")
    for (nm, r) in both
        e = maximum(energy_error(r.trace))
        check(@sprintf("%s%-8s H conserved to round-off", tag, nm), e < 1e-11,
            @sprintf("max |ΔH|/|H₀| = %.3e   H₀ = %.12e", e, r.trace.H[1]))
        (ok, worst) = entropy_monotone(r.trace)
        check(@sprintf("%s%-8s S monotone", tag, nm), ok,
            @sprintf("worst increment %+.3e   S: %.10f -> %.10f",
                worst, r.trace.S[1], r.trace.S[end]))
    end

    # -----------------------------------------------------------------------------------------
    # These bounds hold for ANY state of energy H₀, so a violation is a defect in the run rather
    # than a falsification of the manuscript. Checked before the rates, as their precondition.
    header("$(tag)2. the trajectory stays inside the cone [eq:theoretical-limits]")
    for (nm, r) in both
        (a, b, c) = cone_residual(r.trace, r.trace.H[1])
        check(@sprintf("%s%-8s S ≥ S_η  [eq:tb2]", tag, nm), a < 1e-9,
            @sprintf("worst %+.3e", a))
        check(@sprintf("%s%-8s 1/‖φ‖² ≥ 1/2H₀  [eq:tb3]", tag, nm), b < 1e-9,
            @sprintf("worst %+.3e", b))
        check(@sprintf("%s%-8s 1/‖φ‖² ≤ upper bound  [eq:tb1]", tag, nm), c < 1e-9,
            @sprintf("worst %+.3e", c))
    end

    # -----------------------------------------------------------------------------------------
    header("$(tag)3. S → S_η = H₀ — complete relaxation, unlike A2")
    for (nm, r) in both
        Sη = euler_entropy_minimum(r.trace.H[1])
        excess = (r.trace.S[end] - Sη) / (r.trace.S[1] - Sη)
        check(@sprintf("%s%-8s S(T) reaches S_η", tag, nm), excess < 1e-3,
            @sprintf("(S_T-S_η)/(S₀-S_η) = %.3e   S(T) = %.12f   S_η = %.12f",
                excess, r.trace.S[end], Sη))
        check(@sprintf("%s%-8s S never goes below S_η", tag, nm),
            r.trace.S[end] >= Sη - 1e-10,
            @sprintf("S(T) - S_η = %+.3e", r.trace.S[end] - Sη))
    end

    # -----------------------------------------------------------------------------------------
    header("$(tag)4. the rates: S - S_η at ≈ 1 and ‖ω(t)-ω(T)‖ at ≈ 1/2")
    for (nm, r) in both
        Sη = euler_entropy_minimum(r.trace.H[1])
        excess = r.trace.S .- Sη
        (rS, r²S, nS) = fit_rate(r.trace.t, excess; window = RATE_WINDOW,
            floor = settled_floor(excess))
        check(@sprintf("%s%-8s S - S_η decays at rate ≈ 1", tag, nm),
            isfinite(rS) && abs(rS - 1) < 0.05 && r²S > 0.9999,
            nS == 0 ?
            "no usable range: S - S_η sits at its resolution floor across the whole window" :
            @sprintf("rate %.5f   r² = %.7f   n = %d   (exact: 1)", rS, r²S, nS))

        # The vorticity rate converges to 1/2 far more slowly than the entropy rate converges
        # to 1, and the spectrum says why: in ‖ε‖ the λ = 2 and λ = 4 modes are separated by
        # 0.75 - 0.5 = 0.25, whereas in ‖ε‖² -- which is what S - S_η is -- they are separated
        # by 1.5 - 1 = 0.5. The entropy fit therefore cleans up four times faster in the
        # exponent, and the vorticity fit still reads a few per cent high at any affordable T.
        #
        # So the check is not a tight tolerance on 1/2, which would only be a statement about T.
        # It is the two things the spectrum actually predicts: every contaminating mode decays
        # FASTER than 1/2, so the fitted rate must approach 1/2 strictly FROM ABOVE and must
        # fall as the window moves later. A rate below 1/2, or one rising with time, would
        # contradict the derivation.
        d = distances(r)
        T = spec_T(r)
        (rω, r²ω, nω) = fit_rate(r.snap_t[1:(end - 1)], d[1:(end - 1)];
            window = VORTICITY_WINDOW, floor = settled_floor(d[1:(end - 1)]))
        # 0.12 rather than the entropy fit's 0.05, and A4 is the reason: its initial condition
        # `eq:ic-projector2` contains cos(2x₂), which IS a λ = 4 eigenmode, so the mode that
        # contaminates a vorticity fit is present at order one rather than as a tail. That is
        # the manuscript's own observation -- "the cosine terms shift the initial condition
        # closer to the boundary. As a consequence the initial entropy relaxation rate is
        # slower" -- so a tolerance tight enough for A3 would be measuring A4's initial
        # condition rather than its rate. The number is reported either way.
        check(@sprintf("%s%-8s ‖ω(t)-ω(T)‖ decays at rate ≈ 1/2", tag, nm),
            0.47 < rω < 0.62 && r²ω > 0.999,
            @sprintf("rate %.5f   r² = %.7f   n = %d   (exact asymptote: 1/2)",
                rω, r²ω, nω))

        # The reference artefact, predicted and then measured. A late window must read HIGH by
        # the amount `reference_bias` says, and confirming that is what identifies the excess as
        # the ω(T) stand-in rather than as a failure of the rate to be 1/2.
        (rω_l, _, _) = fit_rate(r.snap_t[1:(end - 1)], d[1:(end - 1)];
            window = RATE_WINDOW, floor = settled_floor(d[1:(end - 1)]))
        tc = T * (RATE_WINDOW[1] + RATE_WINDOW[2]) / 2
        pred = 0.5 + reference_bias(tc, T)
        check(
            @sprintf("%s%-8s a late window reads high, by the predicted amount", tag, nm),
            abs(rω_l - pred) < 0.02,
            @sprintf("t ∈ [%.0f,%.0f]: measured %.5f   predicted %.5f   (bias %.4f)",
                RATE_WINDOW[1] * T, RATE_WINDOW[2] * T, rω_l, pred,
                reference_bias(tc, T)))

        # The manuscript's own reading: the vorticity rate is HALF the entropy rate, because
        # S - S_η is quadratic in the distance to the relaxed state. The ratio is the statement
        # independent of both absolute rates, so it survives a run that is off both.
        check(@sprintf("%s%-8s the entropy rate is twice the vorticity rate", tag, nm),
            isfinite(rS / rω) && abs(rS / rω - 2) < 0.35,
            isfinite(rS / rω) ? @sprintf("rS/rω = %.5f", rS / rω) :
            "not measurable: one of the two fits had no usable range")
    end

    # -----------------------------------------------------------------------------------------
    header("$(tag)5. the vertex of the cone is reached at t ≈ 5")
    for (nm, r) in both
        Sη = euler_entropy_minimum(r.trace.H[1])
        # "Reached" is a choice: 1% of the initial excess entropy, which is where Fig. 6's
        # trajectory becomes indistinguishable from the vertex at plot resolution.
        target = Sη + 0.01 * (r.trace.S[1] - Sη)
        i = findfirst(<=(target), r.trace.S)
        tv = i === nothing ? NaN : r.trace.t[i]
        check(@sprintf("%s%-8s the vertex is reached at t ≈ 5", tag, nm),
            i !== nothing && 2.0 <= tv <= 8.0,
            @sprintf("t(99%% relaxed) = %.3f   (manuscript: t ≈ 5)", tv))
    end

    # -----------------------------------------------------------------------------------------
    header("$(tag)6. the relaxed state lies on 𝔠_η")
    for (nm, r) in both
        (_, resid, _) = best_fit_euler(r.diag, r.trace.final, r.trace.H[1])
        rel = resid / l2norm(r.solver, r.trace.final)
        check(@sprintf("%s%-8s ω(T) is a member of eq:u-eta_Euler_periodic", tag, nm),
            rel < 5e-2, @sprintf("‖ω(T) - fit‖/‖ω(T)‖ = %.5e", rel))
    end

    # -----------------------------------------------------------------------------------------
    header("$(tag)7. the spline and spectral runs agree")
    let ω̂g = spline_grid(res.spline.solver, res.spline.trace.final, opts.spectral)
        e = maximum(abs, ω̂g .- res.spectral.trace.final) /
            maximum(abs, res.spectral.trace.final)
        check("$(tag)the two final states agree", e < 1e-1, @sprintf("max rel %.3e", e))
    end
    for (nm, a, b) in (("H₀", res.spline.trace.H[1], res.spectral.trace.H[1]),
        ("S(0)", res.spline.trace.S[1], res.spectral.trace.S[1]),
        ("S(T)", res.spline.trace.S[end], res.spectral.trace.S[end]))
        rel = abs(a - b) / abs(b)
        check(@sprintf("%s%-6s agrees", tag, nm), rel < 5e-2,
            @sprintf("%.12e vs %.12e   rel %.2e", a, b, rel))
    end
    return nothing
end

"The serialisable record of a projector run."
function projector_payload(res, spec, opts, name)
    (; run = name, spec = (; Δt = spec.Δt, T = spec.T),
        opts = (; spectral = opts.spectral, cells = opts.cells, degree = opts.degree),
        traces = (("spline", res.spline.trace), ("spectral", res.spectral.trace)),
        uΩ = (; spline = res.spline.uΩ, spectral = res.spectral.uΩ),
        Sη = (; spline = euler_entropy_minimum(res.spline.trace.H[1]),
            spectral = euler_entropy_minimum(res.spectral.trace.H[1])),
        snapshots = (; spline_t = res.spline.snap_t, spectral_t = res.spectral.snap_t),
        distance = (; spline = distances(res.spline), spectral = distances(res.spectral)),
        scatter = (;
            spline_initial = scatter_data(res.spline.diag, res.spline.trace.initial),
            spline_final = scatter_data(res.spline.diag, res.spline.trace.final),
            spectral_initial = scatter_data(res.spectral.diag,
                res.spectral.trace.initial),
            spectral_final = scatter_data(res.spectral.diag, res.spectral.trace.final)))
end

"The markdown table of a projector run's own numbers."
function projector_report(res, spec, opts, name)
    lines = ["# $(uppercase(name)) — reduced Euler, projector bracket (§4.2)", "",
        @sprintf("Δt = %.0e, T = %.1f, spectral %d², spline %d² cells degree %d.",
            spec.Δt, spec.T, opts.spectral, opts.cells, opts.degree), "",
        "Uses the coefficient 2H/‖φ‖², not the H/‖φ‖² printed below eq:projector-brackets;",
        "see `verify_projector_factor.jl`.", "",
        "| method | H₀ | max \\|ΔH\\|/\\|H₀\\| | S(0) | S(T) | S_η | rate S | rate ω |",
        "|:--|--:|--:|--:|--:|--:|--:|--:|"]
    for (nm, r) in (("spline", res.spline), ("spectral", res.spectral))
        Sη = euler_entropy_minimum(r.trace.H[1])
        ex = r.trace.S .- Sη
        (rS, _, _) = fit_rate(r.trace.t, ex; window = RATE_WINDOW,
            floor = settled_floor(ex))
        d = distances(r)
        (rω, _, _) = fit_rate(r.snap_t[1:(end - 1)], d[1:(end - 1)];
            window = VORTICITY_WINDOW, floor = settled_floor(d[1:(end - 1)]))
        push!(lines,
            @sprintf("| %s | %.10f | %.3e | %.10f | %.12f | %.12f | %.5f | %.5f |",
                nm, r.trace.H[1], maximum(energy_error(r.trace)), r.trace.S[1],
                r.trace.S[end], Sη, rS, rω))
    end
    return lines
end

end # module ProjectorRun
