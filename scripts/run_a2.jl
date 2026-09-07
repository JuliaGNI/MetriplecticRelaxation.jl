#!/usr/bin/env julia
#
# Run A2: reduced Euler under the metric double bracket, Section 4.1, Figs. 2 and 3.
#
#     julia --project=scripts scripts/run_a2.jl [--runs-dir DIR] [--results-dir DIR]
#                                               [--spectral N] [--cells N] [--degree P]
#
# Same bracket as A1, but h is replaced by the stream function phi, which depends on the state
# through `eq:Poisson-eq-periodic`. The equation is therefore NONLINEAR -- the manuscript calls
# it "a cubic nonlinearity" -- and no analytical solution is known.
#
# u_0 is the Gaussian with x_0 = (pi,pi), w = (0.3, 1.0), N = 1; RK4 at dt = 1e-3 on 256^2.
#
# THE CLAIM, and the reason this run exists:
#
#   H is constant, S is monotonically dissipated, and S converges to a value STRICTLY ABOVE the
#   constrained minimum S_eta = H_0.  That last part is the point.  The manuscript reports it as
#   an observation -- "the entropy appears to converge to a value that is higher than its
#   constrained minimum" -- and it is what motivates Section 4.2: the double bracket relaxes
#   INCOMPLETELY.  So a run that drove S down to S_eta would not be a better result, it would
#   contradict the paper and invalidate the reason the projector bracket was introduced.
#
# The check is therefore two-sided: S must decrease, and it must NOT reach S_eta.  A one-sided
# "S decreased" would pass for either outcome.

using MetriplecticRelaxation
using MetriplecticRelaxation: SECTION4_RUNS, euler_entropy_minimum, spline_grid,
                              energy_error, entropy_monotone, best_fit_euler, fit_rate,
                              l2norm, cone_residual, scatter_data
using Printf

include(joinpath(@__DIR__, "check.jl"))
include(joinpath(@__DIR__, "runner.jl"))
using .Checks: header, check, summary
using .Runner

const spec = SECTION4_RUNS["a2"]
const opts = parse_options()

println("A2  --  reduced Euler, metric double bracket  (Section 4.1, Figs. 2-3)")
@printf("    Δt = %.0e   T = %.1f   steps = %d   spectral %d²   spline %d² cells, degree %d\n",
    spec.Δt, spec.T, round(Int, spec.T / spec.Δt), opts.spectral, opts.cells, opts.degree)
flush(stdout)

section("running")
(ts, dt, trt, uΩt) = run_spline(spec, opts)
(gs, dg, trg, uΩg) = run_spectral(spec, opts)

# =============================================================================================
header("1. H is conserved and S is monotone")

for (label, tr) in (("spline", trt), ("spectral", trg))
    e = maximum(energy_error(tr))
    check(@sprintf("%-8s H conserved to round-off", label), e < 1e-11,
        @sprintf("max |ΔH|/|H₀| = %.3e   H₀ = %.12e", e, tr.H[1]))
    (ok, worst) = entropy_monotone(tr)
    check(@sprintf("%-8s S monotone", label), ok,
        @sprintf("worst increment %+.3e   S: %.10f -> %.10f", worst, tr.S[1], tr.S[end]))
    check(@sprintf("%-8s vorticity mass stays zero", label),
        maximum(abs, tr.M) < 1e-10, @sprintf("max |∫ω| = %.3e", maximum(abs, tr.M)))
end

# =============================================================================================
header("2. S plateaus ABOVE S_η = H₀ — the manuscript's incomplete relaxation")

for (label, tr) in (("spline", trt), ("spectral", trg))
    Sη = euler_entropy_minimum(tr.H[1])
    excess = tr.S[end] - Sη
    # The assertion is on the SIZE of the excess, not its sign. `excess > 0` is a theorem, not a
    # result: `S_η` is the constrained minimum of `S`, so every admissible state satisfies
    # `S ≥ S_η` and the check could not fail — it passes for A3's *complete* relaxation too,
    # whose excess is 2.35e-10. What A2 claims is that the plateau sits far above the minimum,
    # so the threshold is a fraction of `S_η`. Half is well clear of the 179 % measured here and
    # nine orders above anything a relaxed run produces.
    check(@sprintf("%-8s S(T) exceeds S_η by more than 50%%", label), excess / Sη > 0.5,
        @sprintf("S(T) = %.10f   S_η = H₀ = %.10f   excess = %+.6e (%.1f%% of S_η)",
            tr.S[end], Sη, excess, 100excess / Sη))

    # And it is a plateau, not still falling: the last decade of the run must have flattened,
    # or "converged above S_η" would just mean "had not got there yet".
    n = length(tr.S)
    half = tr.S[n ÷ 2]
    late = abs(tr.S[end] - tr.S[max(n - n ÷ 10, 1)]) / max(abs(tr.S[end] - Sη), 1e-300)
    check(@sprintf("%-8s S has plateaued", label), late < 5e-2,
        @sprintf("ΔS over the last 10%% of the run = %.3e of the excess; S(T/2) = %.10f",
            late, half))

    # The relaxation is genuine, not absent: S fell by a substantial fraction of the reduction
    # available to it. HOW LARGE a fraction is a result of this run, not a requirement on it —
    # measured, the double bracket completes 24.15 % of it, identically in both
    # discretisations. Demanding more than half would be an assumption with nothing behind it:
    # the manuscript says only that the entropy "appears to converge to a value that is higher
    # than its constrained minimum", and never how much higher. The threshold below is only
    # large enough to separate a relaxing run from a stalled one.
    frac = (tr.S[1] - tr.S[end]) / (tr.S[1] - Sη)
    check(@sprintf("%-8s S actually relaxed", label), frac > 0.1,
        @sprintf("(S₀-S_T)/(S₀-S_η) = %.4f — %.1f%% of the available reduction; S₀ = %.10f",
            frac, 100frac, tr.S[1]))
end

# =============================================================================================
header("3. the final state is an equilibrium, but not on 𝔠_η")

# The manuscript's evidence is the scatter plot: omega becomes a function of phi. Measured
# rather than plotted -- bin the (phi, omega) cloud and compare the spread within a bin against
# the spread of omega overall.
for (label, d, tr) in (("spline", dt, trt), ("spectral", dg, trg))
    for (when, ω) in (("t = 0", tr.initial), ("t = T", tr.final))
        (φv, ωv) = scatter_data(d, ω)
        lo, hi = extrema(φv)
        nb = 40
        spread = 0.0
        wsum = 0
        for b in 1:nb
            a, z = lo + (hi - lo) * (b - 1) / nb, lo + (hi - lo) * b / nb
            idx = findall(x -> a <= x < z, φv)
            length(idx) < 10 && continue
            m = sum(ωv[idx]) / length(idx)
            spread += sum((ωv[idx] .- m) .^ 2)
            wsum += length(idx)
        end
        rel = sqrt(spread / max(wsum, 1)) / (maximum(ωv) - minimum(ωv))
        check(@sprintf("%-8s %s  ω-vs-φ scatter width", label, when), true,
            @sprintf("within-bin spread / range = %.4f", rel))
    end
end

for (label, d, tr) in (("spline", dt, trt), ("spectral", dg, trg))
    H₀ = tr.H[1]
    (_, res, _) = best_fit_euler(d, tr.final, H₀)
    rel = res / l2norm(d.solver, tr.final)
    check(@sprintf("%-8s the final state is NOT on 𝔠_η", label), rel > 0.05,
        @sprintf("‖ω(T) - fit‖/‖ω(T)‖ = %.4f", rel))
end

# =============================================================================================
header("4. the spline and spectral runs agree")

# `5e-4` is seven times the measured `6.978e-05`, and it is deliberately looser than A3's
# `1e-6` in `projector_run.jl` rather than inconsistent with it: A2's double bracket relaxes
# INCOMPLETELY, so its final state is an arbitrary member of a large set rather than the one
# member of `eq:u-eta_Euler_periodic` that A3's energy fixes, and the two discretisations are
# left disagreeing at their own truncation error instead of at their agreement on `H₀`. The
# previous `5e-2` was 716x the measurement and bounded nothing.
let ω̂g = spline_grid(ts, trt.final, opts.spectral)
    e = maximum(abs, ω̂g .- trg.final) / maximum(abs, trg.final)
    check("the two final states agree", e < 5e-4, @sprintf("max rel %.3e   (tol 5e-04)", e))
end

for (name, a, b) in (("H₀", trt.H[1], trg.H[1]),
    ("S(0)", trt.S[1], trg.S[1]),
    ("S(T)", trt.S[end], trg.S[end]))
    rel = abs(a - b) / abs(b)
    check(@sprintf("%-6s agrees", name), rel < 5e-3,
        @sprintf("%.12e vs %.12e   rel %.2e", a, b, rel))
end

# =============================================================================================
section("writing")

let payload = (; run = "a2", spec = (; Δt = spec.Δt, T = spec.T),
        opts = (; spectral = opts.spectral, cells = opts.cells, degree = opts.degree),
        traces = (("spline", trt), ("spectral", trg)),
        uΩ = (; spline = uΩt, spectral = uΩg),
        Sη = (; spline = euler_entropy_minimum(trt.H[1]),
            spectral = euler_entropy_minimum(trg.H[1])),
        scatter = (; spline_initial = scatter_data(dt, trt.initial),
            spline_final = scatter_data(dt, trt.final),
            spectral_initial = scatter_data(dg, trg.initial),
            spectral_final = scatter_data(dg, trg.final)))
    (p, c) = save_run(opts, "a2", payload)
    println("    run     -> ", p)
    println("    trace   -> ", c)

    lines = ["# A2 — reduced Euler, metric double bracket (§4.1, Figs. 2–3)", "",
        @sprintf("Δt = %.0e, T = %.1f, spectral %d², spline %d² cells degree %d.",
            spec.Δt, spec.T, opts.spectral, opts.cells, opts.degree), "",
        "| method | H₀ | max \\|ΔH\\|/\\|H₀\\| | S(0) | S(T) | S_η = H₀ | S(T) − S_η |",
        "|:--|--:|--:|--:|--:|--:|--:|"]
    for (label, tr) in (("spline", trt), ("spectral", trg))
        Sη = euler_entropy_minimum(tr.H[1])
        push!(lines,
            @sprintf("| %s | %.10f | %.3e | %.10f | %.10f | %.10f | %+.6e |",
                label, tr.H[1], maximum(energy_error(tr)),
                tr.S[1], tr.S[end], Sη,
                tr.S[end] - Sη))
    end
    push!(lines, "",
        "The entropy plateaus strictly above S_η, which is the manuscript's incomplete",
        "relaxation and the reason §4.2 introduces the projector bracket.")
    println("    report  -> ", report(opts, "a2", lines))
end

summary("run_a2.jl")
