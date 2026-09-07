#!/usr/bin/env julia
#
# Run C1: Grad-Shafranov relaxation on the rectangle, Section 5.5, figures `gsr_*`.
#
#     julia --project=scripts scripts/run_c1.jl [--runs-dir DIR] [--results-dir DIR]
#                                               [--samples N]
#
# The collision-like div-grad bracket on Omega = [1,7] x [-9.5,9.5] with the measure
# dmu = dr dz / r, the Herrnegger-Maschke entropy s = y^2/2(Cr^2+D) with C = 0.6, D = 0.2, and
# u_0 the Gaussian with x_0 = (4,0), w_1^2 = 0.5, w_2^2 = 3.2, N = 1 -- SQUARED widths, see
# `gaussian_w2`.  The state variable is j = u/r; `GradShafranovBox` says why.
#
# `--cells` and `--degree` are deliberately NOT read: Section 5.5's mesh is anisotropic, so a
# single number cannot express it, and the pair is a recorded choice in `SECTION55_RUNS`.
#
# THE CLAIM:
#
#   the scatter of u/(Cr^2+D) against psi collapses onto a LINE through the origin whose slope is
#   the eigenvalue lambda of the classical Grad-Shafranov solver; H is conserved to machine
#   precision; S decreases monotonically.
#
# The three are independent.  H is conserved because the bracket is DEGENERATE on it, S falls
# because the bracket is POSITIVE SEMI-DEFINITE, and the slope is lambda because the relaxed
# state is the constrained minimiser -- a run could satisfy any two and fail the third.
#
# The scatter check is therefore TWO-SIDED, as B1's was: the fitted lambda must reach the space's
# own eigenvalue AND the residual ||u/(Cr^2+D) - lambda psi|| must fall to the floor set by the L2
# projection of that field.  lambda alone is a projection coefficient whose error is SECOND order
# in the residual, so it is within a per cent for states that are nothing like an equilibrium --
# 0.0514 against 0.0302 at t = 0 -- and asserting only on it would pass a run that had barely
# moved.
#
# What the residual floor IS matters, and it is not zero.  The discrete relaxed state satisfies
# Pi(sigma j) = lambda psi exactly -- the L2 PROJECTION of the ordinate is proportional to psi --
# so the pointwise residual `verify_gradshafranov.jl` measures is the projection error of
# sigma j, which falls with the mesh (2.2e-3, 4.3e-4, 1.5e-4, 7.0e-5 at 6x7, 10x12, 14x16 and
# 18x21 cells).  The run is asked to REACH that floor, not to beat it.

using MetriplecticRelaxation
using MetriplecticRelaxation: SECTION55_RUNS, gs_state, gs_flow, gs_fit,
                              gs_rayleigh, gs_eigenvalue,
                              GS_LAMBDA_RECTANGLE, GS_LAMBDA_CONTINUUM, TakedaGrid,
                              takeda_iterate, energy_error, entropy_monotone, l2norm,
                              scatter_data
using PoissonBrackets: Integrator, ImplicitMidpoint, integrate_step!, entropy_production,
                       default_f_abstol, nbasis
using Printf

include(joinpath(@__DIR__, "check.jl"))
include(joinpath(@__DIR__, "runner.jl"))
using .Checks: header, check, summary
using .Runner

const spec = SECTION55_RUNS["c1"]
const opts = parse_options(; samples = 200)

println("C1  --  Grad-Shafranov on the rectangle, collision bracket, s = y²/2(Cr²+D)  (§5.5, figs gsr_*)")
@printf("    Δt = %.4g   T = %.1f   steps = %d   %d×%d cells, degree %d   dμ = dr dz / r\n",
    spec.Δt, spec.T, round(Int, spec.T / spec.Δt), spec.cells..., spec.degree)
flush(stdout)

# The reference eigenvalue, computed here by the CLASSICAL solver on the manuscript's own
# 64×64 node grid — no relaxation involved. `verify_takeda.jl` is what establishes that this
# number is right; recomputing it here is what makes the comparison below self-contained.
section("the reference eigenvalue, from the classical iteration")
const takeda = takeda_iterate(TakedaGrid(64, 64))
@printf("    Takeda on 64×64 nodes: λ = %.10f in %d sweeps   printed %.6f   continuum %.10f\n",
    takeda.λ, takeda.iterations, GS_LAMBDA_RECTANGLE, GS_LAMBDA_CONTINUUM)
flush(stdout)

section("running")
const production = Tuple{Float64, Float64}[]
(box, d, flow, tr) = run_gs(spec, opts;
    observer = (b, f, t, ĵ) -> push!(production, (t, entropy_production(f, ĵ))))

const H₀ = tr.H[1]
const λh = gs_eigenvalue(box)
const Sη = λh * H₀

# =============================================================================================
header("1. H is conserved to machine precision and S is monotone")

# The bound is the NEWTON RESIDUAL TOLERANCE, not the method, and it has a derivation rather
# than a number. `H` is quadratic, so the midpoint increment satisfies
# ΔH = ∂H/∂ĵ(j̄)·(ĵⁿ⁺¹−ĵⁿ) exactly; the increment is −Δt 𝔾(j̄) ∂S/∂ĵ(j̄) plus the Newton residual
# ρ, and the first term is annihilated by the degeneracy — so ΔH = ∂H/∂ĵ(j̄)·ρ per step, bounded
# by ‖∂H/∂ĵ‖ times `default_f_abstol` and INDEPENDENT of Δt. Over `n` steps the worst case is
# `n` times that and the random walk is `√n` times it, which is the bound asserted here.
#
# `default_f_abstol` is 4 max(8,N) eps ‖ĵ₀‖_∞ and knows nothing about H₀, so the RELATIVE energy
# error floors near f_abstol/H₀ rather than near eps. Both numbers are reported, because "to
# machine precision" on a given H₀ is a claim about the absolute drift.
let e = maximum(energy_error(tr)), tol = default_f_abstol(Float64, nbasis(box), tr.initial),
    n = round(Int, spec.T / spec.Δt), bound = sqrt(n) * tol / abs(H₀)

    check("H conserved at the Newton residual tolerance, over √n steps", e < bound,
        @sprintf("max |ΔH|/|H₀| = %.3e   max |ΔH| = %.3e   H₀ = %.12e   f_abstol/H₀ = %.2e   √n bound %.2e   ratio to one step %.2f",
            e, e * abs(H₀), H₀, tol / abs(H₀), bound, e / (tol / abs(H₀))))
end

# For a QUADRATIC entropy the midpoint rule dissipates exactly, at any step size: the same
# identity as above gives S(y) − S(x) = −Δt gᵀ𝔾g ≤ 0 by positive semi-definiteness alone. Unlike
# §5.4's B3, monotonicity here is not a statement about Δt — which is why Δt is an accuracy
# choice throughout.
#
# What IS a statement about the run is where the identity stops being measurable. Once S has
# reached λ_h H₀ to twelve digits the evaluated difference of two nearly equal numbers floats at
# round-off and its SIGN is not meaningful. So the claim is split: strictly falling while it
# falls, and no rise beyond round-off after that.
let falling = findall(i -> (tr.S[i] - Sη) / Sη > 1e-11, eachindex(tr.S)),
    (ok, worst) = entropy_monotone(tr; tol = 1e-14)

    check("S falls strictly at every step at which it is still falling",
        all(diff(tr.S[falling]) .< 0),
        @sprintf("%d of %d samples above the floor by more than 1e-11; worst increment among them %+.3e",
            length(falling), length(tr.S),
            length(falling) > 1 ? maximum(diff(tr.S[falling])) : 0.0))
    check("and never rises beyond round-off once it has arrived", ok,
        @sprintf("worst increment over the whole trace %+.3e (relative to S₀)   S: %.10e → %.10e",
            worst, tr.S[1], tr.S[end]))
end

# The same round-off floor, on the same cause: (S,S) = gᵀ𝔾g with g = ∂S/∂ĵ, and at the
# equilibrium g lands in the kernel of 𝔾, so the evaluated form is a difference of large terms
# that cancel. Positivity is asserted where the quantity is above round-off and non-negativity
# relative to its own scale everywhere — the same convention `ispositive_semidefinite` uses, and
# for the same reason: what is being tested is a sign.
let ps = [p[2] for p in production], sc = maximum(ps),
    active = findall(i -> (tr.S[i] - Sη) / Sη > 1e-11, eachindex(tr.S))

    check("(S,S) > 0 at every sample where the state has not yet arrived",
        all(>(0), ps[active]),
        @sprintf("smallest over those %d samples %.4e   largest %.4e — %.1f orders",
            length(active), minimum(ps[active]), sc,
            log10(sc / minimum(ps[active]))))
    check("and (S,S) ≥ 0 to round-off against its own scale, everywhere",
        minimum(ps) ≥ -1e-13 * sc,
        @sprintf("min (S,S) = %+.4e at t = %.3f   max = %.4e   min/max = %+.3e",
            minimum(ps), production[argmin(ps)][1], sc, minimum(ps) / sc))
end

# The mass Casimir ∫u dμ = ∫j dx is NOT conserved discretely, because the constant function is
# not in the homogeneous-Dirichlet space — and that absence is exactly what forces the
# equilibrium multipliers μ and c to zero, which is what makes eq:gs-ref exact. Reported, not
# asserted on.
check("the mass drift is reported, not asserted  [REPORTED]", true,
    @sprintf("∫u dμ: %.10f → %.10f   relative change %+.3e",
        tr.M[1], tr.M[end], (tr.M[end] - tr.M[1]) / tr.M[1]))

# =============================================================================================
header("2. S decreases to S_η = λ_h H₀, the Poincaré floor")

check("S/H ≥ λ_h at every sample — the floor holds off the flow as well as on it",
    minimum(tr.S ./ tr.H) ≥ λh * (1 - 1e-12),
    @sprintf("min S/H = %.12f   λ_h = %.12f   excess %+.3e",
        minimum(tr.S ./ tr.H), λh, minimum(tr.S ./ tr.H) - λh))

let excess = (tr.S[end] - Sη) / Sη
    check("S(T) reaches λ_h H₀ to round-off", abs(excess) < 1e-11,
        @sprintf("(S(T) − λ_h H₀)/λ_h H₀ = %+.4e   S(T) = %.12e   λ_h H₀ = %.12e",
            excess, tr.S[end], Sη))
end

let frac = (tr.S[1] - tr.S[end]) / (tr.S[1] - Sη)
    check("S completed the reduction available to it", frac > 0.99999,
        @sprintf("(S₀−S_T)/(S₀−S_η) = %.12f   S₀ = %.10e", frac, tr.S[1]))
end

# =============================================================================================
header("3. the scatter collapses onto u/(Cr²+D) = λψ")

let (λ₀, _, rel₀) = gs_fit(box, tr.initial), (λ, res, rel) = gs_fit(box, tr.final),
    q = gs_rayleigh(box, tr.final)

    check(
        "the RESIDUAL ‖u/(Cr²+D) − λψ‖/‖u/(Cr²+D)‖ reaches its floor — the two-sided half",
        rel < 1e-4,
        @sprintf("t = T: %.4e (absolute %.4e)   t = 0: %.4e (λ₀ = %.8f)",
            rel, res, rel₀, λ₀))

    # The Rayleigh quotient is stationary at the discrete equilibrium and its value THERE is
    # λ_h exactly, both being the smallest eigenvalue of the pencil (MΛ, W). So this is the
    # sharp test of whether the run arrived, and it is sharper than the fit below.
    check("S/H reaches λ_h — the sharp statement that the run arrived",
        abs(q - λh) / λh < 1e-10,
        @sprintf("S/H = %.14f   λ_h = %.14f   relative %+.3e", q, λh, (q - λh) / λh))

    # The pointwise fit is looser, and by a known amount: the discrete equilibrium satisfies
    # Π(σj) = λψ, so the fit sees the projection error of σj rather than a state error, and its
    # deviation from λ_h is second order in `rel`. The ratio is reported because it is the
    # informative number — a constant threshold here would only record which one was chosen.
    check("and the FITTED λ agrees with λ_h at second order in the residual",
        abs(λ - λh) / λh < 1e-6,
        @sprintf("λ(T) = %.12f   λ_h = %.12f   relative %+.3e   |Δλ|/(λ_h rel²) = %.4f",
            λ, λh, (λ - λh) / λh, abs(λ - λh) / (λh * rel^2)))
end

# =============================================================================================
header("4. against the classical Grad-Shafranov solver")

# THE claim of §5.5: the relaxation and the classical iteration agree. They are two different
# computations of the same eigenvalue on two different discretisations, so what they can agree to
# is the larger of the two discretisation errors — which here is the finite-volume grid's.
let e = abs(λh - takeda.λ) / takeda.λ
    check("the relaxed λ agrees with Takeda's iteration to its own grid error", e < 1e-3,
        @sprintf("relaxation λ_h = %.10f   Takeda 64×64 = %.10f   relative %+.3e",
            λh, takeda.λ, e))
end

let e = abs(λh - GS_LAMBDA_CONTINUUM) / GS_LAMBDA_CONTINUUM
    check("and with the continuum eigenvalue to the spline space's own accuracy", e < 1e-5,
        @sprintf("λ_h = %.12f   continuum %.10f   relative %+.3e",
            λh, GS_LAMBDA_CONTINUUM, e))
end

# The manuscript's printed 0.030302 is 0.22 % above the continuum value, which is its own
# solver's grid error — see `GS_LAMBDA_RECTANGLE`. So this is the loosest of the three
# comparisons by construction, and reporting it loosely is the honest thing to do.
let e = abs(λh - GS_LAMBDA_RECTANGLE) / GS_LAMBDA_RECTANGLE
    check("the printed λ = 0.030302 is reproduced to 0.3 %", e < 3e-3,
        @sprintf("λ_h = %.10f   printed %.6f   relative %+.3e   (the printed value is %+.3e from the continuum one)",
            λh, GS_LAMBDA_RECTANGLE,
            e,
            (GS_LAMBDA_RECTANGLE - GS_LAMBDA_CONTINUUM) / GS_LAMBDA_CONTINUUM))
end

# =============================================================================================
header("5. Δt is an accuracy choice, and this is the measurement behind it")

# Δt, Δt/2 and Δt/4 to the same physical time, on the run's own space, and the RATE at which the
# successive differences fall rather than a threshold on one of them. Crank-Nicolson damps a
# stiff mode by only 4/(λΔt) per step, so too large a Δt does not blow up — it silently slows the
# apparent relaxation as Δt², and nothing else in this script would notice. The observed factor
# must be near 4, the order of implicit midpoint; a LARGER factor means the difference is still
# dominated by the stiff-mode error and the step is above the asymptotic regime. Measured here:
# 9.8, 13.8 and 3.6 for the halvings from Δt = 0.5, which is why the recorded step is 0.0625.
#
# The comparison is at t = 2.5 and not after five steps. The relaxation is essentially over by
# t = 5, so a five-step comparison at this Δt measures the first step's error on a rough initial
# condition rather than the trajectory — and it reads an order of magnitude larger for that
# reason. 40 + 80 + 160 steps, some nine minutes at this resolution.
let n = round(Int, 2.5 / spec.Δt)
    function endpoint(Δt, steps)
        f = gs_flow(box)
        ĵ = gs_state(box, spec)
        integ = Integrator(f, ImplicitMidpoint(), Δt; û₀ = ĵ)
        for _ in 1:steps
            integrate_step!(ĵ, integ)
        end
        return ĵ
    end
    a, b, c = endpoint(spec.Δt, n), endpoint(spec.Δt / 2, 2n), endpoint(spec.Δt / 4, 4n)
    e1 = l2norm(box, a .- b) / l2norm(box, c)
    e2 = l2norm(box, b .- c) / l2norm(box, c)
    check(
        @sprintf("the step-halving difference falls at second order at t = %.4g",
            n * spec.Δt),
        3.0 < e1 / e2 < 5.0,
        @sprintf("‖Δt − Δt/2‖ = %.3e   ‖Δt/2 − Δt/4‖ = %.3e   ratio %.3f",
            e1, e2, e1 / e2))
    check("and it is small at the recorded Δt", e1 < 2e-3,
        @sprintf("relative L² difference %.3e at Δt = %.4g", e1, spec.Δt))
end

# =============================================================================================
section("writing")

let (λ, _, rel) = gs_fit(box, tr.final),
    payload = (; run = "c1",
        spec = (; Δt = spec.Δt, T = spec.T, cells = spec.cells,
            degree = spec.degree),
        opts = (; N = nbasis(box.space)),
        traces = (("gradshafranov", tr),),
        λ = (; fitted = λ, rayleigh = gs_rayleigh(box, tr.final), discrete = λh,
            takeda = takeda.λ, printed = GS_LAMBDA_RECTANGLE,
            continuum = GS_LAMBDA_CONTINUUM),
        Sη = (; discrete = Sη),
        production = production,
        scatter = (; initial = scatter_data(d, tr.initial),
            final = scatter_data(d, tr.final)))

    (p, c) = save_run(opts, "c1", payload)
    println("    run     -> ", p)
    println("    trace   -> ", c)

    lines = ["# C1 — Grad-Shafranov on the rectangle (§5.5, figs `gsr_*`)", "",
        @sprintf("Δt = %.4g, T = %.1f, %d×%d cells degree %d, N = %d, state = j = u/r.",
            spec.Δt, spec.T, spec.cells..., spec.degree, nbasis(box.space)), "",
        "| quantity | value |", "|:--|--:|",
        @sprintf("| H₀ | %.12e |", H₀),
        @sprintf("| max \\|ΔH\\|/\\|H₀\\| | %.3e |", maximum(energy_error(tr))),
        @sprintf("| max \\|ΔH\\| | %.3e |", maximum(energy_error(tr)) * abs(H₀)),
        @sprintf("| S(0) | %.10e |", tr.S[1]),
        @sprintf("| S(T) | %.10e |", tr.S[end]),
        @sprintf("| S_η = λ_h H₀ | %.10e |", Sη),
        @sprintf("| (S(T) − λ_h H₀)/λ_h H₀ | %+.4e |", (tr.S[end] - Sη) / Sη),
        @sprintf("| fitted λ at t = T | %.12f |", λ),
        @sprintf("| Rayleigh S/H at t = T | %.12f |", gs_rayleigh(box, tr.final)),
        @sprintf("| the space's own λ_h | %.12f |", λh),
        @sprintf("| Takeda, 64×64 nodes | %.10f |", takeda.λ),
        @sprintf("| the manuscript's printed λ | %.6f |", GS_LAMBDA_RECTANGLE),
        @sprintf("| the continuum λ | %.10f |", GS_LAMBDA_CONTINUUM),
        @sprintf("| ‖u/(Cr²+D) − λψ‖/‖u/(Cr²+D)‖ at t = T | %.4e |", rel),
        @sprintf("| the same at t = 0 | %.4e |", gs_fit(box, tr.initial)[3]),
        @sprintf("| ∫u dμ, t = 0 → T | %.10f → %.10f |", tr.M[1], tr.M[end]), "",
        "The mass drifts because the constant function is not in the homogeneous-Dirichlet",
        "space — the same absence that forces the equilibrium multipliers to zero and so makes",
        "`u/(Cr²+D) = λψ` exact.", "",
        "The printed λ = 0.030302 is 0.22 % above the continuum eigenvalue 0.0302346260, which",
        "is the discretisation error of the authors' own reference solver; see",
        "`verify_takeda.jl`."]
    println("    report  -> ", report(opts, "c1", lines))
end

summary("run_c1.jl")
