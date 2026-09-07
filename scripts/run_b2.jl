#!/usr/bin/env julia
#
# Run B2: the perturbed equilibrium, Section 5.4, figures `pe_*`.
#
#     julia --project=scripts scripts/run_b2.jl [--runs-dir DIR] [--results-dir DIR]
#                                               [--cells N] [--degree P] [--samples N]
#
# Same domain, same bracket and same entropy as B1.  What differs is the initial condition:
# omega_0 = sin(6 pi x_1) sin(4 pi x_2) + omega_G with N = 100, i.e. an EIGENMODE of -Delta on
# the square plus a 1 % Gaussian bump.
#
# An eigenmode satisfies omega = lambda phi exactly, so it is a stationary state of the
# relaxation -- and, not being the lowest one, an unstable one.  That is what makes this a
# perturbed equilibrium rather than a second arbitrary initial condition.
#
# THE CLAIM, and it is TWO-PHASE:
#
#   S is FLAT first, and only then monotonically decreasing.
#
# Both halves have to be asserted, and neither implies the other.  Monotonicity alone is
# satisfied by a run that fell from the first step -- which is exactly what B1 does, and
# `run_b1.jl` reports the same diagnostic on the same mesh as the control.  Flatness alone is
# satisfied by a run that never moved.  So this driver measures WHEN the entropy starts falling,
# how flat it was before that, and that it never rises anywhere.

using MetriplecticRelaxation
using MetriplecticRelaxation: SECTION5_RUNS, EulerSquare, euler_entropy_floor,
                              eigenmode_fit,
                              dirichlet_eigenvalue, euler_state, euler_flow, state_extrema,
                              energy_error, entropy_monotone, entropy_plateau, l2norm,
                              scatter_data, DIRICHLET_EIGENVALUE
using PoissonBrackets: Integrator, ImplicitMidpoint, integrate_step!, entropy_production,
                       default_f_abstol, nbasis
using Printf

include(joinpath(@__DIR__, "check.jl"))
include(joinpath(@__DIR__, "runner.jl"))
using .Checks: header, check, summary
using .Runner

const spec = SECTION5_RUNS["b2"]
const opts = parse_options(; cells = 26, degree = 2, samples = 400)

println("B2  --  perturbed equilibrium, collision bracket, s = ω²/2  (§5.4, figs pe_*)")
@printf("    Δt = %.4g   T = %.1f   steps = %d   %d² cells, degree %d\n",
    spec.Δt, spec.T, round(Int, spec.T / spec.Δt), opts.cells, opts.degree)
flush(stdout)

section("running")
const production = Tuple{Float64, Float64}[]
(sq, d, flow, tr) = run_euler(spec, opts;
    observer = (s, f, t, ω̂) -> push!(production, (t, entropy_production(f, ω̂))))

const H₀ = tr.H[1]
const λh = dirichlet_eigenvalue(sq)
const Sηh = euler_entropy_floor(H₀; λ = λh)

# =============================================================================================
header("1. H is conserved to machine precision and S never rises")

# The relative energy error floors at the NEWTON RESIDUAL TOLERANCE over H₀, not at eps:
# `default_f_abstol` is absolute and knows nothing about H₀. See `run_b1.jl`, which measures
# both numbers on the same mesh.
let e = maximum(energy_error(tr)),
    tol = default_f_abstol(Float64, nbasis(sq.space), tr.initial)

    check("H conserved to the Newton residual tolerance", e < 1e-11,
        @sprintf("max |ΔH|/|H₀| = %.3e   max |ΔH| = %.3e   H₀ = %.12e   f_abstol = %.2e",
            e, e * abs(H₀), H₀, tol))
end

let (ok, worst) = entropy_monotone(tr)
    check("S monotone over the whole run", ok,
        @sprintf("worst increment %+.3e   S: %.10e → %.10e", worst, tr.S[1], tr.S[end]))
end

check("S ≥ λ_h H₀ at every sample — the Poincaré floor", minimum(tr.S) >= Sηh,
    @sprintf("min S - λ_h H₀ = %+.4e   λ_h H₀ = %.10e", minimum(tr.S) - Sηh, Sηh))

# =============================================================================================
header("2. PHASE ONE: S is flat while the instability grows")

const (t_break, i_break, held) = entropy_plateau(tr)

check("S holds still before it falls", held < 1e-2,
    @sprintf("max |S(t)-S₀|/(S₀-S_T) over t < %.4g is %.3e", t_break, held))

check("the plateau lasts several time units, not one step", t_break >= 4spec.Δt,
    @sprintf("1%% of the total fall reached at t = %.4g (sample %d of %d), i.e. %.1f%% of T",
        t_break, i_break, length(tr.t), 100t_break / spec.T))

# The mechanism, measured rather than asserted: what ends the plateau is the entropy production
# growing, and it grows exponentially. If it were flat the plateau would be the whole run.
let p₀ = production[1][2], pmax = maximum(x -> x[2], production)
    check("(S,S) grows by orders of magnitude — the instability", pmax / p₀ > 1e2,
        @sprintf("(S,S): %.4e at t = 0 → %.4e at the peak (t = %.4g), a factor %.3g",
            p₀, pmax, production[argmax([x[2] for x in production])][1], pmax / p₀))
end

# The initial state IS close to an equilibrium, which is why the plateau exists at all — and how
# close is a property of the mesh, not of the manuscript. The Gaussian is a 1 % perturbation;
# the projection error of a mode with three wavelengths across the domain is the other part, and
# whichever is larger is what actually seeds the instability. Reported, not asserted.
let (λ₀, _, rel₀) = eigenmode_fit(sq, tr.initial)
    check("ω₀ is near an eigenmode of -Δ, and this is how near", rel₀ < 0.2,
        @sprintf("‖ω₀ - λφ₀‖/‖ω₀‖ = %.4e with λ₀ = %.4f; 52π² = %.4f; the Gaussian is 1%%",
            rel₀, λ₀, 52π^2))
end

# =============================================================================================
header("3. PHASE TWO: S decreases to the ground state")

let frac = (tr.S[1] - tr.S[end]) / (tr.S[1] - Sηh)
    check("S completed almost all of the reduction available to it", frac > 0.99,
        @sprintf("(S₀-S_T)/(S₀-λ_h H₀) = %.8f   S₀ = %.10e   S_T = %.10e",
            frac, tr.S[1], tr.S[end]))
end

# The tolerance on λ is the SQUARE of the state residual, not a constant, and that is a measured
# law rather than a loosened threshold: λ is the projection coefficient (ω,φ)/(φ,φ), whose error
# is second order in the state's, and |Δλ|/λ_h = 0.2474 rel² holds to five digits on B1 and B2
# alike — three orders of magnitude apart in rel. So `rel²` is the accuracy the state itself
# permits, a sharper statement than any absolute floor, and one that tightens on its own as a
# run relaxes. B2 does not reach B1's 7e-05 residual within an affordable wall clock: it relaxes
# about four times more slowly, and T is bounded by the step cost.
let (λ, _, rel) = eigenmode_fit(sq, tr.final)
    check("the relaxed state is the LOWEST eigenmode, not the one it started on",
        abs(λ - λh) / λh < rel^2 && rel < 5e-2,
        @sprintf("λ(T) = %.8f   λ_h = %.8f   2π² = %.8f   |Δλ|/λ_h = %.3e against ‖ω-λφ‖/‖ω‖² = %.3e",
            λ, λh, DIRICHLET_EIGENVALUE, abs(λ - λh) / λh, rel^2))
    check("...and it is single-signed, where ω₀ was not",
        state_extrema(sq, tr.final)[1] >= 0 && state_extrema(sq, tr.initial)[1] < 0,
        @sprintf("ω(T) ∈ [%+.4e, %.6f]   ω₀ ∈ [%+.4f, %.4f]",
            state_extrema(sq, tr.final)..., state_extrema(sq, tr.initial)...))
end

# =============================================================================================
header("4. Δt is a cost choice, not an accuracy one")

let n = 10
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

let payload = (; run = "b2", spec = (; Δt = spec.Δt, T = spec.T, entropy = spec.entropy),
        opts = (; cells = opts.cells, degree = opts.degree, N = nbasis(sq.space)),
        traces = (("euler", tr),),
        plateau = (; t_break = t_break, i_break = i_break, held = held),
        λ = (; initial = eigenmode_fit(sq, tr.initial)[1],
            fitted = eigenmode_fit(sq, tr.final)[1], discrete = λh,
            continuum = DIRICHLET_EIGENVALUE, mode = 52π^2),
        Sη = (; discrete = Sηh, continuum = euler_entropy_floor(H₀)),
        production = production,
        scatter = (; initial = scatter_data(d, tr.initial),
            final = scatter_data(d, tr.final)))
    (p, c) = save_run(opts, "b2", payload)
    println("    run     -> ", p)
    println("    trace   -> ", c)

    (λ, _, rel) = eigenmode_fit(sq, tr.final)
    (λ₀, _, rel₀) = eigenmode_fit(sq, tr.initial)
    lines = [
        "# B2 — perturbed equilibrium, collision bracket, s = ω²/2 (§5.4, figs pe_*)", "",
        @sprintf("Δt = %.4g, T = %.1f, %d² cells degree %d, N = %d.",
            spec.Δt, spec.T, opts.cells, opts.degree, nbasis(sq.space)), "",
        "| quantity | value |", "|:--|--:|",
        @sprintf("| H₀ | %.12e |", H₀),
        @sprintf("| max \\|ΔH\\|/\\|H₀\\| | %.3e |", maximum(energy_error(tr))),
        @sprintf("| S(0) | %.10e |", tr.S[1]),
        @sprintf("| S(T) | %.10e |", tr.S[end]),
        @sprintf("| S_η = λ_h H₀ | %.10e |", Sηh),
        @sprintf("| plateau ends at t | %.4g (sample %d) |", t_break, i_break),
        @sprintf("| held before then | %.3e of the total fall |", held),
        @sprintf("| (S,S) at t = 0 | %.4e |", production[1][2]),
        @sprintf("| (S,S) at the peak | %.4e |", maximum(x -> x[2], production)),
        @sprintf("| λ, ‖ω-λφ‖/‖ω‖ at t = 0 | %.4f, %.4e |", λ₀, rel₀),
        @sprintf("| λ, ‖ω-λφ‖/‖ω‖ at t = T | %.8f, %.4e |", λ, rel),
        @sprintf("| the initial mode's eigenvalue 52π² | %.6f |", 52π^2),
        @sprintf("| discrete λ_h | %.8f |", λh),
        @sprintf("| ∫ω, t = 0 → T | %.10f → %.10f |", tr.M[1], tr.M[end]), "",
        "The entropy holds still and then falls, which is the perturbed-equilibrium",
        "signature. `run_b1.jl` reports the same plateau diagnostic on the same mesh for a",
        "run that starts far from any equilibrium, and that is the control."]
    println("    report  -> ", report(opts, "b2", lines))
end

summary("run_b2.jl")
