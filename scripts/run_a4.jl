#!/usr/bin/env julia
#
# Run A4: reduced Euler under the projector bracket, second initial condition. Section 4.2,
# Fig. 7.
#
#     julia --project=scripts scripts/run_a4.jl [--runs-dir DIR] [--results-dir DIR]
#                                               [--spectral N] [--cells N] [--degree P]
#
# `eq:ic-projector2`: u_0 = cos(2 x_2) + u_G with 1/N = 1.8, x_0 = (pi, 3pi/2), w = (0.3, 1).
# Everything else is A3's. The manuscript's point is that "the cosine terms shift the initial
# condition closer to the boundary [of the cone]. As a consequence the initial entropy
# relaxation rate is slower, but approaches ~ 1 as the trajectory approaches the vertex."
#
# So A4 carries one claim A3 does not, and it is checked here: the rate is NOT constant along
# the trajectory. Fitting a single rate over the whole run would report a number that describes
# neither end, which is why `fit_rate`'s window is a recorded choice throughout.
#
# THE INITIAL CONDITION IS AMBIGUOUS, AND FOR THIS RUN IT MATTERS.  `eq:initial_gaussian` is
# written unmodified on T^2, with no summation over periodic images. A4's Gaussian is centred
# pi/2 from the boundary with w_2 = 1 and amplitude 1.8, so the formula as printed still has the
# value 1.8 exp(-(pi/2)^2) = 0.153 there -- a discontinuity of 8.5% of its own peak. A1's is
# 1e-69 and A2/A3's is 5.2e-5, so A4 is the only run affected.
#
# Both readings are therefore run and both are reported: the literal formula as the primary,
# being what is printed, and the periodised form alongside. See `Gaussian` and section 8 below.

using MetriplecticRelaxation
using MetriplecticRelaxation: SECTION4_RUNS, periodise, euler_entropy_minimum, fit_rate,
                              energy_error, l2norm, initial_condition
using Printf

# Only this include: `projector_run.jl` brings in `check.jl` and `runner.jl` itself and
# re-exports them. Including either of those again here would make a second `Options` type.
include(joinpath(@__DIR__, "projector_run.jl"))
using .ProjectorRun
using .ProjectorRun: check_summary

const spec = SECTION4_RUNS["a4"]
const spec_p = periodise(spec)
const opts = parse_options()

println("A4  --  reduced Euler, projector bracket, second IC  (Section 4.2, Fig. 7)")
@printf("    Δt = %.0e   T = %.1f   steps = %d   spectral %d²   spline %d² cells, degree %d\n",
    spec.Δt, spec.T, round(Int, spec.T / spec.Δt), opts.spectral, opts.cells, opts.degree)
@printf("    u₀ = cos(2x₂) + u_G,  1/N = 1.8,  x₀ = (π, 3π/2),  w = (0.3, 1.0)\n")
flush(stdout)

section("running the initial condition as printed")
const res = projector_run(spec, opts)

projector_checks(res, spec, opts)

# =============================================================================================
header("8. the entropy rate is slower early and approaches 1 at the vertex")

# This is the one claim A4 makes that A3 does not. Fitted in two windows: an early one, before
# the trajectory has come off the cone boundary, and a late one approaching the vertex.
#
# THE ASSERTION IS ONLY ON THE LATE RATE. The early one is REPORTED, because the manuscript's
# statement about it and the linearised spectrum point in opposite directions and this run is
# what has to settle which is describing Fig. 7.
#
# The manuscript says "the cosine terms shift the initial condition closer to the boundary. As a
# consequence the initial entropy relaxation rate is slower, but approaches ≈ 1 as the
# trajectory approaches the vertex of the cone." That reasoning is about the PL constant
# κ_η = 2a of `eq:shrunk-cone`, which vanishes as a state approaches the cone boundary — so it
# is a statement about a LOWER BOUND on the rate, not about the rate itself.
#
# The spectrum says the opposite for the actual decay. A4's initial condition contains
# cos(2x₂), which is a λ = 4 eigenmode; its contribution to S - S_η decays at 2(1 - 1/4) = 3/2,
# against the λ = 2 mode's 2(1 - 1/2) = 1. An initial condition loaded with λ = 4 should
# therefore shed entropy FASTER early and slow toward 1, not the reverse.
for (nm, r) in (("spline", res.spline), ("spectral", res.spectral))
    Sη = euler_entropy_minimum(r.trace.H[1])
    excess = r.trace.S .- Sη
    fl = ProjectorRun.settled_floor(excess)
    (r_early, r²e, ne) = fit_rate(r.trace.t, excess; window = (0.01, 0.10), floor = fl)
    (r_late, r²l, nl) = fit_rate(r.trace.t, excess; window = ProjectorRun.RATE_WINDOW,
        floor = fl)
    check(@sprintf("%-8s the late rate approaches 1", nm),
        isfinite(r_late) && abs(r_late - 1) < 0.05 && r²l > 0.9999,
        nl == 0 ?
        "no usable range: S - S_η is at its resolution floor throughout the window" :
        @sprintf("late %.5f   r² = %.6f   n = %d", r_late, r²l, nl))
    check(@sprintf("%-8s early vs late rate  [REPORTED, not asserted]", nm), true,
        @sprintf("early %.5f (r² = %.6f, n = %d)   late %.5f   -> early is %s",
            r_early, r²e,
            ne,
            r_late,
            !isfinite(r_early) || !isfinite(r_late) ? "unmeasurable here" :
            r_early > r_late ? "FASTER (spectrum), not slower (manuscript)" :
            "slower, as the manuscript says"))
end

# =============================================================================================
section("running the periodised initial condition")
const res_p = projector_run(spec_p, opts)

projector_checks(res_p, spec_p, opts; label = "[periodised]")

# =============================================================================================
header("9. the two readings of the initial condition, compared")

# Neither reading is "the" answer -- the manuscript does not say which it means -- so this
# section reports the difference rather than asserting one is right. What it establishes is
# whether the AMBIGUITY changes any conclusion, which is the question that matters for the
# reproduction.
for (nm, a, b) in (("H₀", res.spline.trace.H[1], res_p.spline.trace.H[1]),
    ("S(0)", res.spline.trace.S[1], res_p.spline.trace.S[1]),
    ("S(T)", res.spline.trace.S[end], res_p.spline.trace.S[end]))
    rel = abs(a - b) / abs(b)
    check(@sprintf("%-6s differs between the two readings", nm), true,
        @sprintf("literal %.10f   periodised %.10f   rel %.3e", a, b, rel))
end

for (nm, r, rp) in (("spline", res.spline, res_p.spline),
    ("spectral", res.spectral, res_p.spectral))
    Sη = euler_entropy_minimum(r.trace.H[1])
    Sηp = euler_entropy_minimum(rp.trace.H[1])
    ex, exp_ = r.trace.S .- Sη, rp.trace.S .- Sηp
    (rate, _, _) = fit_rate(r.trace.t, ex; window = ProjectorRun.RATE_WINDOW,
        floor = ProjectorRun.settled_floor(ex))
    (ratep, _, _) = fit_rate(rp.trace.t, exp_; window = ProjectorRun.RATE_WINDOW,
        floor = ProjectorRun.settled_floor(exp_))
    # The claims are what must survive the ambiguity: both readings must give rate ≈ 1.
    check(@sprintf("%-8s BOTH readings give the entropy rate ≈ 1", nm),
        abs(rate - 1) < 0.05 && abs(ratep - 1) < 0.05,
        @sprintf("literal %.5f   periodised %.5f", rate, ratep))
end

# The spline/spectral cross-check is the one thing the discontinuity does damage, and the
# control that identifies it: periodising must bring the two discretisations back together.
let
    e = abs(res.spline.trace.S[1] - res.spectral.trace.S[1]) / res.spectral.trace.S[1]
    ep = abs(res_p.spline.trace.S[1] - res_p.spectral.trace.S[1]) /
         res_p.spectral.trace.S[1]
    check("periodising improves the spline/spectral agreement", ep < e / 5,
        @sprintf("S(0) rel: literal %.3e  ->  periodised %.3e   (%.0fx better)",
            e, ep, e / max(ep, 1e-300)))
end

# =============================================================================================
section("writing")

let (p, c) = save_run(opts, "a4", projector_payload(res, spec, opts, "a4"))
    println("    run     -> ", p)
    println("    trace   -> ", c)
end
let (p, c) = save_run(opts, "a4_periodic",
        projector_payload(res_p, spec_p, opts, "a4_periodic"))
    println("    run     -> ", p)
    println("    trace   -> ", c)
end

let lines = projector_report(res, spec, opts, "a4")
    append!(lines,
        ["", "## The periodised initial condition", "",
            "`eq:initial_gaussian` is written unmodified on T², and A4's Gaussian is centred",
            "π/2 from the boundary with w₂ = 1 and amplitude 1.8, so as printed it is",
            "discontinuous by 1.8·exp(−(π/2)²) = 0.153, i.e. 8.5 % of its peak. Which reading",
            "the manuscript intends is not stated, so both are run.", ""])
    append!(lines, projector_report(res_p, spec_p, opts, "a4 (periodised)")[5:end])
    println("    report  -> ", report(opts, "a4", lines))
end

check_summary("run_a4.jl")
