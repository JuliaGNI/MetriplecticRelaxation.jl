#!/usr/bin/env julia
#
# The diagnostics, checked before any run quotes a number through them.
#
#     julia --project=scripts scripts/verify_diagnostics.jl
#
# The one that matters is section 2. `best_fit_euler` replaces the manuscript's "best fit ...
# varying the three phases" with a closed-form projection, on the grounds that the three phases
# parameterise the unit sphere in four orthogonal mode coefficients. If that reduction is wrong
# the A3 and A4 residuals are wrong, and nothing else in the reproduction would notice. It is
# therefore checked against a brute-force scan over the phases, which is the thing it claims to
# replace.

using MetriplecticRelaxation
using MetriplecticRelaxation: SpectralTorus, SplineTorus, Diagnostics, torus_field,
                              spectral_state, spline_state, energy, entropy, potential,
                              potential_norm², best_fit_euler, fit_rate, Trace, record!,
                              cone_residual, energy_error, entropy_monotone,
                              euler_minimiser, euler_entropy_minimum, l2inner, l2norm,
                              mean_value, integrate, poisson_periodic, SECTION4_RUNS
using PoissonBrackets: project
using Printf
using Random

include(joinpath(@__DIR__, "check.jl"))
using .Checks: header, check, summary

Random.seed!(0x5c1e9a3b)

const g = SpectralTorus(48)

# =============================================================================================
header("1. the energies and the entropy")

let spec = SECTION4_RUNS["a2"], d = Diagnostics(g, spec)
    ω, uΩ = spectral_state(g, spec)
    φ = potential(d, ω)
    check("H = ½(φ,ω) matches ½‖∇φ‖² by parts",
        abs(energy(d, ω) - l2inner(g, φ, ω) / 2) < 1e-13,
        @sprintf("H = %.10f", energy(d, ω)))
    check("S = ½‖ω‖²", abs(entropy(d, ω) - l2inner(g, ω, ω) / 2) < 1e-13,
        @sprintf("S = %.10f", entropy(d, ω)))
    check("the vorticity has zero mass", abs(integrate(g, ω)) < 1e-13,
        @sprintf("%.2e", abs(integrate(g, ω))))
    # eq:tb2 and eq:tb3, which hold for any state of this energy.
    check("S ≥ H₀  [eq:tb2]", entropy(d, ω) >= energy(d, ω),
        @sprintf("S = %.6f ≥ H₀ = %.6f", entropy(d, ω), energy(d, ω)))
    check("‖φ‖² ≤ 2H₀  [eq:tb3]", potential_norm²(d, ω) <= 2energy(d, ω) + 1e-12,
        @sprintf("‖φ‖² = %.6f ≤ 2H₀ = %.6f", potential_norm²(d, ω), 2energy(d, ω)))
end

let spec = SECTION4_RUNS["a1"], d = Diagnostics(g, spec)
    ω, _ = spectral_state(g, spec)
    check("A1's energy is linear in ω",
        abs(energy(d, 2 .* ω) - 2energy(d, ω)) / abs(energy(d, ω)) < 1e-12,
        @sprintf("H(2ω)/H(ω) = %.14f", energy(d, 2 .* ω) / energy(d, ω)))
    check("A1's potential is h - h_Ω, independent of ω",
        potential(d, ω) === potential(d, 3 .* ω), "")
end

# =============================================================================================
header("2. best_fit_euler is the minimiser over the three phases")

let spec = SECTION4_RUNS["a3"], d = Diagnostics(g, spec)
    ω, _ = spectral_state(g, spec)
    H₀ = energy(d, ω)
    fit, res, coef = best_fit_euler(d, ω, H₀)

    # (a) the fit is IN the family: it is a member of eq:u-eta_Euler_periodic, so its own
    #     entropy is S_η = H₀ exactly.
    Sfit = l2inner(g, fit, fit) / 2
    check("the fit lies on 𝔠_η: S(fit) = S_η = H₀",
        abs(Sfit - euler_entropy_minimum(H₀)) / H₀ < 1e-12,
        @sprintf("S(fit) = %.12f   H₀ = %.12f", Sfit, H₀))

    # (b) the brute-force scan over the three phases cannot beat it.
    best = Inf
    bestθ = (0.0, 0.0, 0.0)
    n = 60
    for i in 0:(n - 1), j in 0:(n - 1), k in 0:(n - 1)
        θ₀, θ₁, θ₂ = 2π * i / n, 2π * j / n, 2π * k / n
        w = torus_field(g, euler_minimiser(H₀, θ₀, θ₁, θ₂))
        r = l2inner(g, ω .- w, ω .- w)
        if r < best
            best = r
            bestθ = (θ₀, θ₁, θ₂)
        end
    end
    scan = sqrt(best)
    check("closed form ≤ a 60³ scan over (θ₀,θ₁,θ₂)", res <= scan + 1e-12,
        @sprintf("closed form %.10f   scan %.10f   (scan is worse by %.2e)",
            res, scan, scan - res))
    # And the scan gets close, or the closed form would be minimising something else.
    check("the scan approaches the closed form", (scan - res) / res < 1e-3,
        @sprintf("rel gap %.2e   at θ = (%.3f, %.3f, %.3f)",
            (scan - res) / res, bestθ...))

    # (c) refining the scan reduces the gap, which is what says the gap is the scan's own
    #     resolution rather than a defect in the closed form.
    best2 = Inf
    n2 = 120
    for i in 0:(n2 - 1), j in 0:(n2 - 1), k in 0:(n2 - 1)
        w = torus_field(g,
            euler_minimiser(H₀, 2π * i / n2, 2π * j / n2, 2π * k / n2))
        r = l2inner(g, ω .- w, ω .- w)
        r < best2 && (best2 = r)
    end
    check("doubling the scan resolution halves the gap",
        (sqrt(best2) - res) < (scan - res) / 1.5,
        @sprintf("gap %.3e -> %.3e", scan - res, sqrt(best2) - res))

    # (d) a state already in the family is fitted exactly.
    let w = torus_field(g, euler_minimiser(H₀, 0.7, 1.3, 2.9))
        _, r2, _ = best_fit_euler(d, w, H₀)
        check("a state on 𝔠_η is fitted to round-off", r2 / l2norm(g, w) < 1e-13,
            @sprintf("residual/‖ω‖ = %.2e", r2 / l2norm(g, w)))
    end
end

# =============================================================================================
header("3. fit_rate recovers a known exponential")

let t = collect(0.0:0.01:10.0)
    for λ in (0.5, 1.0, 2.0)
        y = 3.7 .* exp.(-λ .* t)
        (r, r², n) = fit_rate(t, y)
        check(@sprintf("rate %.1f recovered", λ), abs(r - λ) < 1e-10,
            @sprintf("%.12f   r² = %.10f   n = %d", r, r², n))
    end
    # A non-exponential decay must show up in r², or a quoted rate would mean nothing.
    let y = 1 ./ (1 .+ t) .^ 2
        (r, r², _) = fit_rate(t, y)
        check("an algebraic decay has poor r²", r² < 0.999,
            @sprintf("rate %.4f   r² = %.6f", r, r²))
    end
    # The floor drops round-off tails rather than fitting them.
    let y = max.(3.7 .* exp.(-1.0 .* t), 1e-17)
        (r, _, _) = fit_rate(t, y)
        check("a round-off floor does not drag the slope", abs(r - 1.0) < 1e-8,
            @sprintf("%.12f", r))
    end
end

# =============================================================================================
header("4. the cone and the monotonicity test")

let spec = SECTION4_RUNS["a3"], d = Diagnostics(g, spec)
    ω, _ = spectral_state(g, spec)
    H₀ = energy(d, ω)
    tr = Trace(ω)
    # A synthetic trajectory that relaxes onto the family: it must lie inside the cone.
    fit, _, _ = best_fit_euler(d, ω, H₀)
    for (i, s) in enumerate(range(0, 1; length = 40))
        w = (1 - s) .* ω .+ s .* fit
        # Rescale to hold the energy fixed, which is what the cone assumes.
        w .*= sqrt(H₀ / energy(d, w))
        record!(tr, d, 10.0 * (i - 1) / 39, w)
    end
    (a, b, c) = cone_residual(tr, H₀)
    check("S ≥ S_η along the path", a < 1e-12, @sprintf("worst %.2e", a))
    check("1/‖φ‖² ≥ 1/2H₀  [eq:tb3]", b < 1e-12, @sprintf("worst %.2e", b))
    check("1/‖φ‖² ≤ upper bound  [eq:tb1]", c < 1e-12, @sprintf("worst %.2e", c))

    (ok, worst) = entropy_monotone(tr)
    check("the synthetic path dissipates entropy monotonically", ok,
        @sprintf("worst increment %.2e", worst))
    check("energy_error is at round-off on an energy-fixed path",
        maximum(energy_error(tr)) < 1e-12,
        @sprintf("max %.2e", maximum(energy_error(tr))))
end

# A control: a path that INCREASES entropy must be caught.
let spec = SECTION4_RUNS["a3"], d = Diagnostics(g, spec)
    ω, _ = spectral_state(g, spec)
    tr = Trace(ω)
    for i in 1:10
        record!(tr, d, 0.1i, (1 + 0.01i) .* ω)
    end
    (ok, worst) = entropy_monotone(tr)
    check("CONTROL: a rising entropy is caught", !ok,
        @sprintf("worst increment %.4f", worst))
end

# =============================================================================================
header("5. the spline diagnostics agree with the spectral ones")

let ts = SplineTorus(64, 3)
    for name in ("a1", "a2", "a3", "a4")
        spec = SECTION4_RUNS[name]
        dg, dt = Diagnostics(g, spec), Diagnostics(ts, spec)
        ωg, _ = spectral_state(SpectralTorus(64), spec)
        ω̂, _ = spline_state(ts, spec)
        dg64 = Diagnostics(SpectralTorus(64), spec)

        eH = abs(energy(dg64, ωg) - energy(dt, ω̂)) / abs(energy(dg64, ωg))
        eS = abs(entropy(dg64, ωg) - entropy(dt, ω̂)) / abs(entropy(dg64, ωg))
        tol = name == "a4" ? 5e-2 : 5e-3
        check(@sprintf("%s: H agrees", name), eH < tol, @sprintf("rel %.2e", eH))
        check(@sprintf("%s: S agrees", name), eS < tol, @sprintf("rel %.2e", eS))
    end
end

summary("verify_diagnostics.jl")
