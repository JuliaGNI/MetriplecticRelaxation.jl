#!/usr/bin/env julia
#
# A discrepancy inside Section 4.2 of the manuscript, settled numerically.
#
#     julia --project=scripts scripts/verify_projector_factor.jl
#
# THE CLAIM UNDER TEST.  Just below `eq:projector-brackets` the manuscript writes the
# projector-bracket evolution of the reduced Euler case as
#
#     d_t u = -[ u - u_Omega - H(u)/||phi||^2 phi ]                        (as printed)
#
# but its own definitions give
#
#     d_t u = -[ u - u_Omega - 2 H(u)/||phi||^2 phi ]                      (as derived)
#
# because `eq:L2-projector` sets c(u,v) = ||phi||^{-2} (phi, v) and `eq:Euler_H_periodic`
# defines H = (1/2)(phi,u), so with v = delta S/delta u = omega,
#
#     c = ||phi||^{-2} (phi, omega) = 2 H / ||phi||^2 .
#
# This is INTERNAL to Section 4.2: `eq:SS-projector` states (S,S) = 2S - 4 H_0^2/||phi||^2, and
# only the factor 2 reproduces it.  The two printed equations cannot both be right.
#
# The analytic test case one page earlier is NOT affected, and section 3 below confirms that:
# `eq:analytical_H` defines H = (h - h_Omega, u) with no factor 1/2, so there c = H/||h||^2 and
# the printed equation is correct.  That asymmetry is what makes this look like a transcription
# slip rather than a different convention -- a different convention would have moved both.
#
# WHAT WOULD GO WRONG UNQUESTIONED.  The factor-1 field is still symmetric and still dissipates
# entropy, so a run using it looks healthy in the two plots the manuscript shows.  What it does
# NOT do is conserve the energy: section 2 measures dH/dt at order one rather than at round-off.
# It would relax to a state of the wrong energy, and hence to the wrong member of the family
# `eq:u-eta_Euler_periodic`.

using MetriplecticRelaxation
using MetriplecticRelaxation: SpectralTorus, torus_field, poisson_periodic,
                              l2inner, l2norm, mean_value
using Printf

include(joinpath(@__DIR__, "check.jl"))
using .Checks: header, check, relerr, summary

const g = SpectralTorus(48)

# A field with no special relation to the Laplacian's eigenspaces: if omega were an
# eigenfunction then phi would be parallel to omega, the projector would annihilate it, and
# every quantity below would be zero for a reason unrelated to the factor.
const Ω = let u = torus_field(g,
        (a, b) -> 0.9cos(a) + 0.5sin(2b) + 0.4cos(a - b) + 0.2sin(3a + b) + 0.3cos(2a + 3b))
    u .- mean_value(g, u)
end
const Φ = poisson_periodic(g, Ω)
const H₀ = l2inner(g, Φ, Ω) / 2
const S₀ = l2inner(g, Ω, Ω) / 2

"The field as printed in the manuscript, with coefficient `κ H/‖φ‖²`."
field(κ) = .-(Ω .- (κ * H₀ / l2inner(g, Φ, Φ)) .* Φ)

# =============================================================================================
header("1. the projector is a projector, and phi is not parallel to omega")

check("ω is not an eigenfunction of -Δ (so the test is not vacuous)",
    abs(l2inner(g, Φ, Ω)) / (l2norm(g, Φ) * l2norm(g, Ω)) < 0.99,
    @sprintf("cos∠(φ,ω) = %.4f", l2inner(g, Φ, Ω) / (l2norm(g, Φ) * l2norm(g, Ω))))

"``\\Pi_H v = v - \\|\\phi\\|^{-2}(\\phi, v)\\phi``, `eq:L2-projector`."
Π(v) = v .- (l2inner(g, Φ, v) / l2inner(g, Φ, Φ)) .* Φ

let Πω = Π(Ω)
    # `Π_H φ = 0` is REPORTED and not asserted, and APPLYING the projector to `φ` rather than
    # writing the formula out at `v = φ` does not make it assertable: `Π`'s body IS that
    # formula, so at `v = φ` the coefficient is `x/x`, exactly `1.0` in floating point, and
    # `φ .- 1.0 .* φ` is exactly zero elementwise. Measured, `all(iszero, Π(ψ))` holds for
    # `ψ = φ`, for a random field, and for `φ` scaled by `1e±9`: no field could make this fail.
    # The content is in the two rows below -- only the right coefficient makes `Π_H ω`
    # orthogonal to `φ`, and only a projector is idempotent.
    println(@sprintf("      Π_H φ = 0 holds by the formula's own algebra   ‖Π_H φ‖/‖φ‖ = %.2e",
        l2norm(g, Π(Φ)) / l2norm(g, Φ)))
    # The projected field is orthogonal to φ, which is what "projector onto the orthogonal
    # complement" means and is not true of an arbitrary v.
    check("Π_H ω ⊥ φ", abs(l2inner(g, Φ, Πω)) / (l2norm(g, Φ) * l2norm(g, Πω)) < 1e-13,
        @sprintf("cos∠ = %.2e", abs(l2inner(g, Φ, Πω)) / (l2norm(g, Φ) * l2norm(g, Πω))))
    check("Π_H is idempotent", l2norm(g, Π(Πω) .- Πω) / l2norm(g, Πω) < 1e-13,
        @sprintf("%.2e", l2norm(g, Π(Πω) .- Πω) / l2norm(g, Πω)))
    # And it is not the identity, or the three lines above would hold of nothing.
    check("Π_H is not the identity", l2norm(g, Πω .- Ω) / l2norm(g, Ω) > 0.1,
        @sprintf("‖Π_H ω - ω‖/‖ω‖ = %.4f", l2norm(g, Πω .- Ω) / l2norm(g, Ω)))
end

# =============================================================================================
header("2. only the factor 2 conserves the energy")

for κ in (1.0, 2.0)
    f = field(κ)
    dH = l2inner(g, Φ, f)
    rel = abs(dH) / (l2norm(g, Φ) * l2norm(g, f))
    label = κ == 2 ? "κ = 2 (derived): dH/dt = 0" : "κ = 1 (as printed): dH/dt ≠ 0"
    check(label, κ == 2 ? rel < 1e-14 : rel > 1e-2,
        @sprintf("dH/dt = %+.6e   normalised %.2e", dH, rel))
end

# =============================================================================================
header("3. only the factor 2 reproduces eq:SS-projector")

let SS = 2S₀ - 4H₀^2 / l2inner(g, Φ, Φ)
    for κ in (1.0, 2.0)
        dS = l2inner(g, Ω, field(κ))
        rel = abs(-dS - SS) / abs(SS)
        label = κ == 2 ? "κ = 2 (derived): -dS/dt = (S,S)" :
                "κ = 1 (as printed): -dS/dt ≠ (S,S)"
        check(label, κ == 2 ? rel < 1e-13 : rel > 1e-2,
            @sprintf("-dS/dt = %.10e   (S,S) = %.10e   rel %.2e", -dS, SS, rel))
    end
end

# Both fields dissipate entropy, which is why the error is not visible in an entropy trace.
# That is REPORTED and not asserted, because as an assertion it could not fail. With
# `H₀ = (φ,ω)/2` and `S₀ = (ω,ω)/2`, `field(κ)` gives the closed form
#
#     -dS/dt = 2S₀ - 2κH₀²/‖φ‖² ,
#
# and Cauchy-Schwarz bounds `4H₀² = (φ,ω)² ≤ ‖φ‖²·2S₀`, i.e. `H₀²/‖φ‖² ≤ S₀/2`, so
# `-dS/dt ≥ (2-κ)S₀`. Both κ here satisfy `κ ≤ 2`, and section 1 has already shown `φ ∦ ω`,
# which makes the inequality strict: no field could have made either row fail. Measured,
# `H₀²/(‖φ‖²S₀) = 0.3967` against the bound `1/2`. The sign is not unconditional -- on this
# same field `κ = 3` and `κ = 4` give `-dS/dt = -5.067` and `-15.638` -- it is unconditional
# over the two κ the manuscript's discrepancy is between.
#
# What IS asserted is the closed form, which a sign or factor slip in `field` breaks at order
# one and which is what makes the bound above a statement about this code rather than about
# the algebra alone.
for κ in (1.0, 2.0)
    dS = l2inner(g, Ω, field(κ))
    closed = 2S₀ - 2κ * H₀^2 / l2inner(g, Φ, Φ)
    check(@sprintf("κ = %g: -dS/dt = 2S₀ - 2κH₀²/‖φ‖²", κ), relerr(-dS, closed) < 1e-13,
        @sprintf("-dS/dt = %+.10e   closed form %+.10e   rel %.2e", -dS, closed,
            relerr(-dS, closed)))
    println(@sprintf("      κ = %g dissipates, invisibly in Fig. 4   dS/dt = %+.6e   floor %+.6e",
        κ, dS, (2 - κ) * S₀))
end

# =============================================================================================
header("4. the analytic test case is unaffected: there the printed factor is right")

# `eq:analytical_H`: H(u) = (h - h_Ω, u), linear, no factor 1/2. So the projector coefficient
# is c = ||h-h_Ω||^{-2}(h-h_Ω, ω) = H/||h-h_Ω||^2 — exactly as printed.
# The test field here is NOT `Ω`. `h = cos²x₁ sin²x₂` expands to
# ¼(1 + cos2x₁ - cos2x₂ - ½cos(2x₁+2x₂) - ½cos(2x₁-2x₂)), whose modes are disjoint from `Ω`'s —
# so `(h - h_Ω, Ω)` would be exactly zero by orthogonality, `dH/dt = (κ-1)H` would vanish for
# EVERY κ, and the control below would pass without testing anything. `Ω₄` shares three of `h`'s
# four modes, which is what gives the κ = 2 line something to fail on.
let hf = torus_field(g, (a, b) -> cos(a)^2 * sin(b)^2), hz = hf .- mean_value(g, hf),
    Ω₄ = let u = Ω .+ torus_field(g,
            (a, b) -> 0.6cos(2a) - 0.35cos(2b) + 0.25cos(2a + 2b))
        u .- mean_value(g, u)
    end, Hlin = l2inner(g, hz, Ω₄), nh² = l2inner(g, hz, hz)

    check("the test field overlaps h's spectrum (or this section is vacuous)",
        abs(Hlin) / (l2norm(g, hz) * l2norm(g, Ω₄)) > 0.05,
        @sprintf("H = %.6f   normalised %.4f", Hlin,
            abs(Hlin) / (l2norm(g, hz) * l2norm(g, Ω₄))))

    for κ in (1.0, 2.0)
        f = .-(Ω₄ .- (κ * Hlin / nh²) .* hz)
        dH = l2inner(g, hz, f)
        rel = abs(dH) / (l2norm(g, hz) * l2norm(g, f))
        label = κ == 1 ? "κ = 1 (as printed): dH/dt = 0" : "κ = 2: dH/dt ≠ 0"
        check(label, κ == 1 ? rel < 1e-14 : rel > 1e-2,
            @sprintf("dH/dt = %+.6e   normalised %.2e", dH, rel))
    end
end

summary("verify_projector_factor.jl")
