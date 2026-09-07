#!/usr/bin/env julia
#
# Run B1: the single vortex, Section 5.4, figures `sv_*`.
#
#     julia --project=scripts scripts/run_b1.jl [--runs-dir DIR] [--results-dir DIR]
#                                               [--cells N] [--degree P] [--samples N]
#
# Reduced Euler on [0,1]^2 with homogeneous Dirichlet boundary conditions, under the
# collision-like div-grad bracket of eq:collision-bracket with s = omega^2/2, hence mobility
# M = 1 by eq:M-condition.  omega_0 is the Gaussian with x_0 = (1/2,1/2), w_1^2 = 0.01,
# w_2^2 = 0.07, N = 1 -- SQUARED widths, see `gaussian_w2`.
#
# THE CLAIM:
#
#   omega relaxes to the lowest Dirichlet eigenmode, omega = lambda_11 phi with
#   lambda_11 = 2 pi^2 = 19.7392; H is conserved to machine precision; S decreases monotonically
#   to S_eta = lambda_11 H_0.
#
# The three are not one claim.  H is conserved because the bracket is DEGENERATE on it, S falls
# because the bracket is POSITIVE SEMI-DEFINITE, and the limit is the lowest eigenmode because
# that is the constrained minimiser -- a run could satisfy any two and fail the third.  So the
# check on the final state is TWO-SIDED: the fitted lambda must reach the eigenvalue AND the
# residual ||omega - lambda phi|| must vanish.  The fitted lambda is a projection coefficient
# whose error is SECOND order in that residual, so it is within a per cent of the eigenvalue for
# states that are nothing like the eigenmode -- 22.58 against 19.74 at t = 0 -- and asserting
# only that would pass for a run that had barely moved.

using MetriplecticRelaxation
using MetriplecticRelaxation: SECTION5_RUNS, EulerSquare, euler_entropy_floor,
                              eigenmode_fit,
                              dirichlet_eigenvalue, euler_state, euler_flow, state_extrema,
                              energy_error, entropy_monotone, entropy_plateau, l2norm,
                              l2inner, integrate, scatter_data, DIRICHLET_EIGENVALUE
using PoissonBrackets: Integrator, ImplicitMidpoint, integrate_step!, entropy_production,
                       default_f_abstol, nbasis, degree, ncells
using Printf

include(joinpath(@__DIR__, "check.jl"))
include(joinpath(@__DIR__, "runner.jl"))
using .Checks: header, check, summary
using .Runner

const spec = SECTION5_RUNS["b1"]
# 26 cells and degree 2 are recorded choices, not the manuscript's 64² P2. What rules out the
# manuscript's mesh is the nonlocal bracket's Jacobian — N dense assemblies per Newton matrix,
# see `EulerSquare`. What fixes the number at 26 is the narrow direction, w₁ = 0.1: it is the
# coarsest mesh on which the L² projection of the initial Gaussian stops oscillating below zero,
# measured in `verify_euler.jl`. All three runs share it.
const opts = parse_options(; cells = 26, degree = 2, samples = 400)

println("B1  --  single vortex, collision bracket, s = ω²/2  (§5.4, figs sv_*)")
@printf("    Δt = %.4g   T = %.1f   steps = %d   %d² cells, degree %d\n",
    spec.Δt, spec.T, round(Int, spec.T / spec.Δt), opts.cells, opts.degree)
flush(stdout)

section("running")
const production = Tuple{Float64, Float64}[]
(sq, d, flow, tr) = run_euler(spec, opts;
    observer = (s, f, t, ω̂) -> push!(production, (t, entropy_production(f, ω̂))))

const H₀ = tr.H[1]
const λh = dirichlet_eigenvalue(sq)
const Sη = euler_entropy_floor(H₀)
const Sηh = euler_entropy_floor(H₀; λ = λh)

# =============================================================================================
header("1. H is conserved to machine precision and S is monotone")

# What bounds the drift here is the NEWTON RESIDUAL TOLERANCE, not the method, and the tolerance
# is absolute: `default_f_abstol` is 4 max(8,N) eps ‖ω̂₀‖_∞, which at N = 676 is 6.0e-13 and
# knows nothing about H₀ = 5.6e-4. So the *relative* energy error floors near
# f_abstol × scale / H₀ rather than near eps. Measured, the absolute drift is 1.4e-15 — machine
# precision on a quantity of this size — and the relative one 2.4e-12. The threshold is set from
# that measurement and reports both numbers, because "to machine precision" on a small H₀ is a
# claim about the absolute drift.
let e = maximum(energy_error(tr)),
    tol = default_f_abstol(Float64, nbasis(sq.space), tr.initial)

    check("H conserved to the Newton residual tolerance", e < 1e-11,
        @sprintf("max |ΔH|/|H₀| = %.3e   max |ΔH| = %.3e   H₀ = %.12e   f_abstol = %.2e",
            e, e * abs(H₀), H₀, tol))
end

let (ok, worst) = entropy_monotone(tr)
    check("S monotone", ok,
        @sprintf("worst increment %+.3e   S: %.10e → %.10e", worst, tr.S[1], tr.S[end]))
end

check("(S,S) > 0 at every sample", all(p -> p[2] > 0, production),
    @sprintf("min (S,S) = %.4e at t = %.2f   max = %.4e",
        minimum(p[2] for p in production),
        production[argmin([p[2] for p in production])][1],
        maximum(p[2] for p in production)))

# The mass is a Casimir of the continuous bracket and is NOT conserved discretely, because the
# constant function is not in the homogeneous-Dirichlet space. That absence is exactly what makes
# the closed-form reference below exact, so the drift is reported and not asserted on.
check("the mass drift is reported, not asserted", true,
    @sprintf("∫ω: %.10f → %.10f   relative change %+.3e",
        tr.M[1], tr.M[end], (tr.M[end] - tr.M[1]) / tr.M[1]))

# B1 dissipates from the first step, and this number is the CONTROL for B2's plateau: the same
# diagnostic on the same mesh, on a run that starts nowhere near an equilibrium.
let (tb, ib, held) = entropy_plateau(tr)
    check("S starts falling immediately — the control for B2's plateau", ib <= 2,
        @sprintf("1%% of the total fall reached at t = %.4g (sample %d); held = %.3e",
            tb, ib, held))
end

# =============================================================================================
header("2. S decreases to S_η = λ₁,₁ H₀")

check("S ≥ S_η at every sample — the Poincaré floor", minimum(tr.S) >= Sη,
    @sprintf("min S - S_η = %+.4e   S_η = λ₁,₁H₀ = %.10e", minimum(tr.S) - Sη, Sη))

let excess = (tr.S[end] - Sηh) / Sηh
    check("S(T) is within 1e-5 of the space's own floor λ_h H₀", excess < 1e-5,
        @sprintf("(S(T) - λ_h H₀)/λ_h H₀ = %+.4e   S(T) = %.12e   λ_h H₀ = %.12e",
            excess, tr.S[end], Sηh))
end

let frac = (tr.S[1] - tr.S[end]) / (tr.S[1] - Sηh)
    check("S completed the reduction available to it", frac > 0.999,
        @sprintf("(S₀-S_T)/(S₀-S_η) = %.8f   S₀ = %.10e", frac, tr.S[1]))
end

# =============================================================================================
header("3. the relaxed state IS ω = λ₁,₁ φ")

let (λ₀, _, rel₀) = eigenmode_fit(sq, tr.initial), (λ, _, rel) = eigenmode_fit(sq, tr.final)
    check("the FITTED λ reaches the space's eigenvalue λ_h",
        abs(λ - λh) / λh < 1e-6,
        @sprintf("λ(T) = %.10f   λ_h = %.10f   2π² = %.10f   rel %.3e",
            λ, λh, DIRICHLET_EIGENVALUE, abs(λ - λh) / λh))
    check("the fitted λ reaches 2π² to the space's own accuracy",
        abs(λ - DIRICHLET_EIGENVALUE) / DIRICHLET_EIGENVALUE < 1e-4,
        @sprintf("λ(T) - 2π² = %+.3e   λ_h - 2π² = %+.3e", λ - DIRICHLET_EIGENVALUE,
            λh - DIRICHLET_EIGENVALUE))
    check("the RESIDUAL ‖ω - λφ‖/‖ω‖ vanishes — the two-sided half", rel < 1e-3,
        @sprintf("t = T: %.4e   t = 0: %.4e (λ₀ = %.6f)", rel, rel₀, λ₀))
end

let (lo, hi) = state_extrema(sq, tr.final)
    check("the relaxed vorticity is single-signed, as sin(πx₁)sin(πx₂) is", lo >= 0,
        @sprintf("ω(T) ∈ [%+.4e, %.6f]", lo, hi))
end

# =============================================================================================
header("4. Δt is a cost choice, not an accuracy one")

# The measurement behind the recorded Δt = 1: five steps at Δt against ten at Δt/2, on the run's
# own space. If those disagree, Δt was chosen wrongly and every number above is suspect.
let n = 5
    function endpoint(Δt, steps)
        f = euler_flow(sq, spec)
        ω̂ = euler_state(sq, spec)
        integ = Integrator(f, ImplicitMidpoint(), Δt; û₀ = ω̂)
        for _ in 1:steps
            integrate_step!(ω̂, integ)
        end
        return ω̂
    end
    a, b = endpoint(spec.Δt, n), endpoint(spec.Δt / 2, 2n)
    e = l2norm(sq, a .- b) / l2norm(sq, b)
    check(@sprintf("Δt and Δt/2 agree at t = %.4g", n * spec.Δt), e < 1e-3,
        @sprintf("relative L² difference %.3e", e))
end

# =============================================================================================
section("writing")

let payload = (; run = "b1", spec = (; Δt = spec.Δt, T = spec.T, entropy = spec.entropy),
        opts = (; cells = opts.cells, degree = opts.degree, N = nbasis(sq.space)),
        traces = (("euler", tr),),
        λ = (; fitted = eigenmode_fit(sq, tr.final)[1], discrete = λh,
            continuum = DIRICHLET_EIGENVALUE),
        Sη = (; continuum = Sη, discrete = Sηh),
        production = production,
        scatter = (; initial = scatter_data(d, tr.initial),
            final = scatter_data(d, tr.final)))
    (p, c) = save_run(opts, "b1", payload)
    println("    run     -> ", p)
    println("    trace   -> ", c)

    (λ, _, rel) = eigenmode_fit(sq, tr.final)
    lines = ["# B1 — single vortex, collision bracket, s = ω²/2 (§5.4, figs sv_*)", "",
        @sprintf("Δt = %.4g, T = %.1f, %d² cells degree %d, N = %d.",
            spec.Δt, spec.T, opts.cells, opts.degree, nbasis(sq.space)), "",
        "| quantity | value |", "|:--|--:|",
        @sprintf("| H₀ | %.12e |", H₀),
        @sprintf("| max \\|ΔH\\|/\\|H₀\\| | %.3e |", maximum(energy_error(tr))),
        @sprintf("| S(0) | %.10e |", tr.S[1]),
        @sprintf("| S(T) | %.10e |", tr.S[end]),
        @sprintf("| S_η = 2π² H₀ | %.10e |", Sη),
        @sprintf("| S_η = λ_h H₀ | %.10e |", Sηh),
        @sprintf("| (S(T) − λ_h H₀)/λ_h H₀ | %+.4e |", (tr.S[end] - Sηh) / Sηh),
        @sprintf("| fitted λ at t = T | %.10f |", λ),
        @sprintf("| discrete λ_h | %.10f |", λh),
        @sprintf("| continuum λ₁,₁ = 2π² | %.10f |", DIRICHLET_EIGENVALUE),
        @sprintf("| ‖ω(T) − λφ(T)‖/‖ω(T)‖ | %.4e |", rel),
        @sprintf("| ∫ω, t = 0 → T | %.10f → %.10f |", tr.M[1], tr.M[end]), "",
        "The mass drifts because the constant function is not in the homogeneous-Dirichlet",
        "space — the same absence that makes ω = λ₁,₁φ the exact relaxed state."]
    println("    report  -> ", report(opts, "b1", lines))
end

summary("run_b1.jl")
