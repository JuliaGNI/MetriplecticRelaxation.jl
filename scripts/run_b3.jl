#!/usr/bin/env julia
#
# Run B3: the Gibbs entropy, Section 5.4, figures `ge_*`.
#
#     julia --project=scripts scripts/run_b3.jl [--runs-dir DIR] [--results-dir DIR]
#                                               [--cells N] [--degree P] [--samples N]
#
# Same domain and same bracket as B1 and B2.  What differs is the entropy: s = y log y rather
# than y^2/2, so by eq:M-condition the mobility is M = y rather than 1, and omega_0 carries
# 1/N = 10.
#
# THE CLAIM:
#
#   omega relaxes to omega = exp(lambda phi - 1) with lambda = (M(omega) + S(omega)) / 2 H_0.
#
# THE TIME STEP IS THE VARIABLE UNDER TEST, and that is the whole difference from B1 and B2.
# For s = y^2/2 the entropy is a quadratic form and Crank-Nicolson dissipates it monotonically
# whatever the step; for s = y log y it is not, and the manuscript says only that "sufficiently
# small time steps must be used" without saying what that means.  So section 0 below does not
# pick a step and assert -- it sweeps, and reports where monotonicity holds and where it fails.
#
# THE INITIAL CONDITION IS NOT THE PRINTED ONE, and that is a FINDING rather than a convenience.
# s = y log y is undefined at y <= 0 and eq:M-condition gives it the mobility M = y, which the
# same equation requires to be positive.  The printed Gaussian is strictly positive as a
# function -- but with w_1^2 = 0.01 it decays to 1.4e-11 of its peak at the far corner, which a
# Galerkin scheme cannot tell from zero: measured, the projected state starts at 6e-13 and the
# FIRST step drives the iterate to -6.6e-5, at which point the entropy raises a DomainError --
# no step completes.  B3 therefore carries a positive background of 1 % of its peak, `B3_FLOOR`,
# which is four orders above the undershoot it has to absorb.  See `B3_FLOOR` for what it costs.
#
# THE STATE SPACE IS THE SAME AS B1's AND B2's, and the free alternative runs here as a CONTROL.
# Both spaces are viable once the floor makes the state admissible, and they disagree exactly as
# the equilibrium condition predicts: in V_D the mass is not conserved, so the multiplier mu is
# forced to zero and the manuscript's lambda = (M+S)/2H_0 is an identity; in V the mass IS a
# discrete Casimir, mu is fixed by it, and the formula acquires an extra term.  Section 4 runs
# both on a coarse mesh and reports the two mu's side by side.

const CONTROL_CELLS = 16
const CONTROL_STEPS = 25
const CONTROL_DT = 0.02

using MetriplecticRelaxation
using MetriplecticRelaxation: SECTION5_RUNS, EulerSpec, EulerSquare, euler_state,
                              euler_flow,
                              state_extrema, energy_error, entropy_monotone, l2norm,
                              gibbs_lambda, gibbs_fit, gibbs_residual,
                              scatter_data, entropy, energy, vorticity_mass, Diagnostics,
                              B3_FLOOR
using PoissonBrackets: Integrator, ImplicitMidpoint, integrate_step!, entropy_production,
                       default_f_abstol, nbasis
using Printf

include(joinpath(@__DIR__, "check.jl"))
include(joinpath(@__DIR__, "runner.jl"))
using .Checks: header, check, summary
using .Runner

const spec = SECTION5_RUNS["b3"]
const opts = parse_options(; cells = 26, degree = 2, samples = 400)
const h = 1.0 / opts.cells

println("B3  --  Gibbs entropy s = ω log ω, collision bracket  (§5.4, figs ge_*)")
@printf("    Δt = %.4g   T = %.1f   steps = %d   %d² cells, degree %d   state = :%s\n",
    spec.Δt, spec.T, round(Int, spec.T / spec.Δt), opts.cells, opts.degree, spec.state)
flush(stdout)

# =============================================================================================
header("0. the initial condition has to be admissible before anything else is said")

const sq0 = EulerSquare(opts.cells, opts.degree; state = spec.state)
const bare = EulerSpec("b3-printed", spec.section, spec.entropy, spec.state,
    spec.gaussian, nothing, spec.Δt, spec.T)

# The printed initial condition, first, because the floor is a departure from it and has to be
# justified by measurement rather than by assertion.
let (lo, _) = state_extrema(sq0, euler_state(sq0, bare))
    check("the PRINTED ω₀ has no admissible margin at all", lo < 1e-9,
        @sprintf("min ω₀ = %+.4e — at or below the Gaussian's own far-corner value, w₁ = 0.1",
            lo))
    # And "no margin" means "fails", not "is tight".
    fb = euler_flow(sq0, bare)
    ω̂b = euler_state(sq0, bare)
    integ = Integrator(fb, ImplicitMidpoint(), spec.Δt; û₀ = ω̂b)
    (broke, why) = try
        for _ in 1:5
            integrate_step!(ω̂b, integ)
        end
        (state_extrema(sq0, ω̂b)[1] <= 0,
            @sprintf("min ω = %+.4e after 5 steps", state_extrema(sq0, ω̂b)[1]))
    catch err
        (err isa DomainError,
            string(typeof(err)) * ": " * sprint(showerror, err)[1:min(60,
                end)])
    end
    check("CONTROL: the printed ω₀ leaves the admissible set within five steps", broke, why)
end

let (lo, hi) = state_extrema(sq0, euler_state(sq0, spec))
    ok = check(
        @sprintf("the FLOORED ω₀ is admissible with a margin on the %d² space", opts.cells),
        lo > 1e-3, @sprintf("min ω₀ = %+.4e   max = %.4f   floor = %.3g", lo, hi, B3_FLOOR))
    ok || (println("\nω₀ is not admissible: y log y is undefined. Refine and re-run.");
        summary("run_b3.jl"))
end

header("0b. the time step is the variable under test")

# The sweep. Ten steps at each step size, on the run's own mesh: monotone dissipation is a
# per-step property, so a short horizon settles it and a long one only costs.
const STEPS = 10
const SWEEP = (0.005, 0.01, 0.02, 0.04, 0.08, 0.16, 0.5)
const sweep = NamedTuple[]

for Δt in SWEEP
    f = euler_flow(sq0, spec)
    d = Diagnostics(sq0, spec)
    ω̂ = euler_state(sq0, spec)
    S₀ = entropy(d, ω̂)
    H₀ = energy(d, ω̂)
    integ = Integrator(f, ImplicitMidpoint(), Δt; û₀ = ω̂)
    worst = -Inf
    drift = 0.0
    admissible = true
    Sprev = S₀
    for _ in 1:STEPS
        try
            integrate_step!(ω̂, integ)
        catch err
            err isa DomainError || rethrow()
            admissible = false
            break
        end
        state_extrema(sq0, ω̂)[1] > 0 || (admissible = false; break)
        S = entropy(d, ω̂)
        worst = max(worst, (S - Sprev) / abs(S₀))
        Sprev = S
        drift = max(drift, abs(energy(d, ω̂) - H₀) / abs(H₀))
    end
    push!(sweep, (; Δt, worst, drift, admissible, S = Sprev, S₀))
    # A row is a MEASUREMENT, not a claim, and a row that fails is the point of the sweep rather
    # than a defect — so it is printed with `true` and the outcome in the detail, the same way
    # `run_a2.jl` reports its scatter widths. It carries `[REPORTED]` because without it a row
    # reading `[PASS] Δt = 0.16  LEFT THE ADMISSIBLE SET` inverts the harness's own marker. The
    # claims are the three assertions below.
    check(
        @sprintf("Δt = %-7.4g  %s  [REPORTED]", Δt,
            admissible && worst <= 0 ? "monotone" : "LEFT THE ADMISSIBLE SET"),
        true,
        @sprintf("worst ΔS/S₀ = %+.3e   max |ΔH|/|H₀| = %.3e   over %d steps",
            worst, drift, STEPS))
end

# Three claims about the sweep, and the first is what makes the other two mean anything: if the
# outcome were not monotone in Δt there would be no threshold to locate, only a scatter.
let ok = map(r -> r.admissible && r.worst <= 0, sweep),
    bad = filter(r -> !r.admissible || r.worst > 0, sweep),
    good = filter(r -> r.admissible && r.worst <= 0, sweep)

    # `ok` must be a run of trues followed by a run of falses — no interleaving.
    threshold = findfirst(!, ok)
    check("the outcome is monotone in Δt, so there IS a threshold",
        threshold === nothing || !any(ok[threshold:end]),
        @sprintf("outcomes over %s: %s", string(SWEEP),
            join(map(o -> o ? "ok" : "FAIL", ok), " ")))

    check("the sweep brackets that threshold: some Δt work and some do not",
        !isempty(bad) && !isempty(good),
        @sprintf("monotone up to Δt = %.4g; first failure at Δt = %s",
            isempty(good) ? NaN : maximum(r -> r.Δt, good),
            isempty(bad) ? "none in the sweep" : @sprintf("%.4g", minimum(r -> r.Δt, bad))))

    check(@sprintf("the run's own Δt = %.4g is on the monotone side", spec.Δt),
        isempty(bad) || spec.Δt < minimum(r -> r.Δt, bad),
        @sprintf("Δt = %.4g against the first failure at %s", spec.Δt,
            isempty(bad) ? "none" : @sprintf("%.4g", minimum(r -> r.Δt, bad))))
end

# =============================================================================================
section("running")
const production = Tuple{Float64, Float64}[]
const admissibility = Tuple{Float64, Float64}[]
(sq, d, flow, tr) = run_euler(spec, opts;
    observer = (s, f, t, ω̂) -> (push!(production, (t, entropy_production(f, ω̂)));
        push!(admissibility, (t, state_extrema(s, ω̂)[1]))))

const H₀ = tr.H[1]
const Mfinal = tr.M[end]
const Sfinal = tr.S[end]
const λ = gibbs_lambda(Mfinal, Sfinal, H₀)

# =============================================================================================
header("1. H is conserved, S is monotone, and ω stays admissible")

let e = maximum(energy_error(tr)),
    tol = default_f_abstol(Float64, nbasis(sq.space), tr.initial)

    check("H conserved to the Newton residual tolerance", e < 1e-10,
        @sprintf("max |ΔH|/|H₀| = %.3e   max |ΔH| = %.3e   H₀ = %.12e   f_abstol = %.2e",
            e, e * abs(H₀), H₀, tol))
end

let (ok, worst) = entropy_monotone(tr)
    check("S monotone at the chosen Δt", ok,
        @sprintf("worst increment %+.3e   S: %.10e → %.10e", worst, tr.S[1], tr.S[end]))
end

check("ω > 0 at every sample — y log y is defined throughout",
    all(p -> p[2] > 0, admissibility),
    @sprintf("min over the run = %+.4e at t = %.4g", minimum(p[2] for p in admissibility),
        admissibility[argmin([p[2] for p in admissibility])][1]))

check("(S,S) > 0 at every sample", all(p -> p[2] > 0, production),
    @sprintf("(S,S): %.4e at t = 0 → %.4e at t = T", production[1][2], production[end][2]))

# The constants are NOT in this space, so the mass is not a discrete Casimir and drifts — the
# same absence that forces μ = 0 and makes the manuscript's λ formula an identity. Reported, as
# in B1 and B2.
check("the mass drift is reported, not asserted  [REPORTED]", true,
    @sprintf("∫ω: %.10f → %.10f   relative change %+.3e",
        tr.M[1], tr.M[end], (tr.M[end] - tr.M[1]) / tr.M[1]))

# =============================================================================================
header("2. λ = (M(ω) + S(ω)) / 2H₀, and the multiplier μ it assumes away")

const (λfit, μfit, rfit) = gibbs_fit(sq, tr.final; margin = 2h)

# What is actually an identity is the THREE-quantity relation. Substituting
# log ω = λφ + μ - 1 into S = ∫ω log ω gives S = 2λH + (μ-1)M for ANY μ; the manuscript's
# λ = (M+S)/2H₀ is that at μ = 0. So the identity is checked with the fitted μ, and the
# manuscript's special case is then measured against it rather than assumed.
# The tolerance is the FIT's own residual and not a constant: the identity holds exactly for a
# state that satisfies log ω = λφ + μ - 1, and `rfit` is how far this state is from doing so.
# Demanding better than the fit would be demanding the state be something it is not.
let predicted = 2λfit * H₀ + (μfit - 1) * Mfinal,
    rel = abs(predicted - Sfinal) / abs(Sfinal)

    check("S = 2λH₀ + (μ-1)M closes to the accuracy of the fit itself",
        rel < max(2rfit, 1e-6),
        @sprintf("S = %.10f   2λH₀+(μ-1)M = %.10f   rel %.3e   against a fit residual %.3e",
            Sfinal, predicted, rel, rfit))
end

check("some member of the family fits at all", rfit < 0.1,
    @sprintf("weighted rms of log ω - λφ - μ + 1, over the spread of log ω = %.4e", rfit))

# The manuscript's formula is the μ = 0 case, and in this space μ is not free: the constants are
# absent, so nothing constrains the mass and the equilibrium carries no mass multiplier. That is
# the prediction, and this is the measurement of it.
check("the fitted mass multiplier μ has relaxed toward zero", abs(μfit / λfit) < 5e-2,
    @sprintf("μ = %+.8f   λ_fit = %.8f   μ/λ = %+.3e", μfit, λfit, μfit / λfit))

check("λ from the formula agrees with λ from the fit", abs(λfit - λ) / abs(λfit) < 5e-2,
    @sprintf("(M+S)/2H₀ = %.8f   fitted λ = %.8f   rel %.3e", λ, λfit,
        abs(λfit - λ) / abs(λfit)))

check("λ is reported with the quantities it was built from  [REPORTED]", true,
    @sprintf("M(ω) = %.10f   S(ω) = %.10f   H₀ = %.10f   →   λ = %.8f",
        Mfinal, Sfinal, H₀, λ))

# =============================================================================================
header("3. ω = exp(λφ + μ − 1)")

# The manuscript's reference is e^{λφ-1} with λ from the formula; the fitted pair is the control
# that says whether any member of the family fits better. Both are reported at four margins,
# because e^{λφ-1} is e^{-1} on ∂Ω while every ω_h of this space is zero there: the outermost
# cell is where no state of V_D can match the reference, and excluding it is what separates that
# boundary layer from a disagreement in the interior.
const residuals = [(m, gibbs_residual(sq, tr.final, λ; margin = m),
                       gibbs_residual(sq, tr.final, λfit, μfit; margin = m))
                   for m in (0.0, h, 2h, 4h)]

check("the interior state IS the manuscript's e^{λφ-1}", residuals[3][2] < 0.1,
    @sprintf("‖ω - e^{λφ-1}‖/‖ω‖ at margins 0, h, 2h, 4h = %.4e, %.4e, %.4e, %.4e",
        (r[2] for r in residuals)...))

check("the fitted pair does no better, which is what μ ≈ 0 means",
    residuals[3][3] < 2 * residuals[3][2],
    @sprintf("‖ω - e^{λφ+μ-1}‖/‖ω‖ at margins 0, h, 2h, 4h = %.4e, %.4e, %.4e, %.4e",
        (r[3] for r in residuals)...))

check("the disagreement is the boundary layer — the residual falls with the margin",
    residuals[1][2] > residuals[2][2] > residuals[3][2],
    @sprintf("0 → h → 2h: %.4e → %.4e → %.4e; e^{λφ-1} is %.6f on ∂Ω and every ω_h is 0 there",
        residuals[1][2], residuals[2][2], residuals[3][2], exp(-1)))

check("the relaxation moved the state toward the reference",
    residuals[3][2] < gibbs_residual(sq, tr.initial, λ; margin = 2h),
    @sprintf("interior residual, t = 0: %.4e → t = T: %.4e",
        gibbs_residual(sq, tr.initial, λ; margin = 2h), residuals[3][2]))

# =============================================================================================
header("4. CONTROL: the same run in the free space, where the mass IS a Casimir")

# The choice of space is a claim, so it gets a control rather than a paragraph. On a coarse mesh
# and a short horizon — this is about which multiplier survives, not about the relaxed state —
# the two spaces are stepped side by side. In V_D the constants are absent, nothing constrains
# the mass, and μ must relax to zero; in V the mass is a genuine discrete Casimir, μ is fixed by
# it, and the manuscript's λ = (M+S)/2H₀ cannot hold.
const control = NamedTuple[]

for st in (:dirichlet, :free)
    sqc = EulerSquare(CONTROL_CELLS, opts.degree; state = st)
    dc = Diagnostics(sqc, spec)
    fc = euler_flow(sqc, spec)
    ω̂c = euler_state(sqc, spec)
    Hc = energy(dc, ω̂c)
    Mstart = vorticity_mass(dc, ω̂c)
    integ = Integrator(fc, ImplicitMidpoint(), CONTROL_DT; û₀ = ω̂c)
    for _ in 1:CONTROL_STEPS
        integrate_step!(ω̂c, integ)
    end
    hc = 1.0 / CONTROL_CELLS
    (λc, μc, _) = gibbs_fit(sqc, ω̂c; margin = 2hc)
    λfc = gibbs_lambda(vorticity_mass(dc, ω̂c), entropy(dc, ω̂c), Hc)
    push!(control,
        (; state = st, λfit = λc, μ = μc, λformula = λfc,
            mass = (Mstart, vorticity_mass(dc, ω̂c)),
            residual = gibbs_residual(sqc, ω̂c, λfc; margin = 2hc)))
    check(
        @sprintf("%-10s the mass %s, and μ = %+.5f", string(st),
            st === :dirichlet ? "drifts" : "is held", μc),
        st === :dirichlet ?
        abs(vorticity_mass(dc, ω̂c) - Mstart) / Mstart > 1e-3 :
        abs(vorticity_mass(dc, ω̂c) - Mstart) / Mstart < 1e-9,
        @sprintf("∫ω %.8f → %.8f   λ_fit = %.5f   (M+S)/2H₀ = %.5f   residual %.3e",
            Mstart, vorticity_mass(dc, ω̂c), λc, λfc,
            gibbs_residual(sqc, ω̂c, λfc; margin = 2hc)))
end

let d = control[1], fr = control[2]
    check("CONTROL: μ is smaller where the mass is not conserved",
        abs(d.μ) < abs(fr.μ),
        @sprintf("μ: %+.5f in V_D against %+.5f in V, after the same %d steps at Δt = %.3g",
            d.μ, fr.μ, CONTROL_STEPS, CONTROL_DT))
    check("CONTROL: and the manuscript's reference fits better there too",
        d.residual < fr.residual,
        @sprintf("‖ω - e^{λφ-1}‖/‖ω‖: %.4e in V_D against %.4e in V", d.residual,
            fr.residual))
end

# =============================================================================================
header("5. Δt is resolved, measured on this run's own state")

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
    check(@sprintf("Δt and Δt/2 agree at t = %.4g", n * spec.Δt), e < 1e-2,
        @sprintf("relative L² difference %.3e", e))
end

# =============================================================================================
section("writing")

let payload = (; run = "b3",
        spec = (; Δt = spec.Δt, T = spec.T, entropy = spec.entropy, state = spec.state),
        opts = (; cells = opts.cells, degree = opts.degree, N = nbasis(sq.space)),
        traces = (("euler", tr),),
        sweep = sweep,
        control = control,
        floor = B3_FLOOR,
        λ = (; formula = λ, fit = (λfit, μfit, rfit)),
        residuals = residuals,
        production = production,
        admissibility = admissibility,
        scatter = (; initial = scatter_data(d, tr.initial),
            final = scatter_data(d, tr.final)))
    (p, c) = save_run(opts, "b3", payload)
    println("    run     -> ", p)
    println("    trace   -> ", c)

    lines = ["# B3 — Gibbs entropy s = ω log ω, collision bracket (§5.4, figs ge_*)", "",
        @sprintf("Δt = %.4g, T = %.1f, %d² cells degree %d, N = %d, state space :%s.",
            spec.Δt, spec.T, opts.cells, opts.degree, nbasis(sq.space), spec.state), "",
        "## The step-size sweep — the variable under test", "",
        "| Δt | worst ΔS/S₀ over 10 steps | max \\|ΔH\\|/\\|H₀\\| | ω admissible |",
        "|--:|--:|--:|:--|"]
    for r in sweep
        push!(lines,
            @sprintf("| %.4g | %+.3e | %.3e | %s |", r.Δt, r.worst, r.drift,
                r.admissible ? "yes" : "**no**"))
    end
    append!(lines,
        ["", "## The relaxation at Δt = " * @sprintf("%.4g", spec.Δt), "",
            "| quantity | value |", "|:--|--:|",
            @sprintf("| H₀ | %.12e |", H₀),
            @sprintf("| max \\|ΔH\\|/\\|H₀\\| | %.3e |", maximum(energy_error(tr))),
            @sprintf("| S(0) | %.10e |", tr.S[1]),
            @sprintf("| S(T) | %.10e |", Sfinal),
            @sprintf("| M(ω(0)) → M(ω(T)) | %.10f → %.10f |", tr.M[1], Mfinal),
            @sprintf("| λ = (M+S)/2H₀, the manuscript's | %.8f |", λ),
            @sprintf("| fitted λ | %.8f |", λfit),
            @sprintf("| fitted μ | %+.8f |", μfit),
            @sprintf("| μ/λ | %+.4e |", μfit / λfit),
            @sprintf("| relative fit residual | %.4e |", rfit),
            @sprintf("| min ω over the run | %+.4e |",
                minimum(p[2] for p in admissibility))])
    for (m, r0, rf) in residuals
        push!(lines,
            @sprintf("| ‖ω - e^{λφ+μ-1}‖/‖ω‖ at margin %.4f — μ=0 / fitted μ | %.4e / %.4e |",
                m, r0, rf))
    end
    append!(lines,
        ["", "## The state-space control", "",
            @sprintf("%d steps at Δt = %.3g on a %d² mesh, in each space.",
                CONTROL_STEPS, CONTROL_DT, CONTROL_CELLS), "",
            "| space | ∫ω | fitted λ | fitted μ | (M+S)/2H₀ | ‖ω − e^{λφ−1}‖/‖ω‖ |",
            "|:--|:--|--:|--:|--:|--:|"])
    for c in control
        push!(lines,
            @sprintf("| `:%s` | %.8f → %.8f | %.5f | %+.5f | %.5f | %.4e |",
                c.state, c.mass[1], c.mass[2], c.λfit, c.μ, c.λformula, c.residual))
    end
    append!(lines,
        ["",
            "Two findings rather than results, and both are about the manuscript.", "",
            @sprintf("**The printed initial condition is not admissible.** `s = y log y` needs"),
            "`ω > 0` and `eq:M-condition` gives it the mobility `M = y`, which the same equation",
            "requires to be positive. The printed Gaussian is strictly positive as a function,",
            "but with `w₁² = 0.01` it decays to 1.4e-11 of its peak at the far corner, which a",
            "Galerkin scheme cannot tell from zero: the projected state starts at ~6e-13 and the",
            "first step drives the iterate to -6.6e-5, where the entropy raises a `DomainError`",
            @sprintf("before any step completes. B3 therefore carries a %.3g",
                B3_FLOOR),
            "background — 1 % of its peak — and the run above measures the printed condition",
            "failing before it uses the floored one.", "",
            "**Which space the vorticity lives in decides whether `λ = (M+S)/2H₀` holds.** The",
            "formula is the μ = 0 case of `S = 2λH + (μ−1)M`, and μ is the multiplier of the",
            "mass Casimir. In the homogeneous-Dirichlet space the constants are absent, nothing",
            "constrains the mass, and μ relaxes to zero; in the free space the mass is a genuine",
            "discrete Casimir and μ does not. The control table above is that comparison."])
    println("    report  -> ", report(opts, "b3", lines))
end

summary("run_b3.jl")
