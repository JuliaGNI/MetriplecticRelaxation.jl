#!/usr/bin/env julia
#
# The Fourier spectral layer, checked against PoissonBrackets' own independent implementation
# of the same operators.
#
#     julia --project=scripts scripts/verify_spectral.jl
#
# The point of section 1 is that this file's FFT-based derivatives and PoissonBrackets'
# `spectral_grid`, which assembles a dense differentiation MATRIX D = Re(V diag(ik) W), are two
# genuinely different computations of the same operator. Agreement to round-off is what lets
# the rest of the reproduction treat this layer as the reference it claims to be, rather than
# as a second copy of the manuscript's code checked against itself.
#
# Sections 2-4 check the properties the Section 4 claims rest on: the Poisson solve, the
# identity between the nested-bracket and divergence forms of `eq:parallel-diffusion`, and the
# semi-discrete conservation laws -- energy exactly, entropy monotonically.

using MetriplecticRelaxation
using MetriplecticRelaxation: SpectralTorus, torus_field, ∂₁, ∂₂, laplacian,
                              poisson_periodic, canonical_bracket, hamiltonian_field,
                              double_bracket_field, projector_bracket_field,
                              integrate, l2inner, l2norm, mean_value, islands_h,
                              SECTION4_RUNS, spectral_state
using PoissonBrackets: spectral_grid, ∂x, ∂y
using LinearAlgebra
using Printf

include(joinpath(@__DIR__, "check.jl"))
using .Checks: header, check, summary

const N = 32
const g = SpectralTorus(N)
const gp = spectral_grid(N)

# Fixed, band-limited test fields, sharing Fourier modes so that no inner product between them
# is zero by orthogonality alone.
u(a, b) = 0.7cos(a) + 0.4sin(2b) + 0.3cos(a - b) + 0.15sin(3a + b) + 1.1
v(a, b) = 0.5cos(a) + 0.3sin(2b) + 0.25cos(a - b) - 0.2sin(a + b) + 0.6
hf(a, b) = cos(a) + 0.6sin(b) + 0.25cos(a + b)

const U = torus_field(g, u)
const V = torus_field(g, v)
const H = torus_field(g, hf)

# =============================================================================================
header("1. the FFT operators agree with PoissonBrackets' differentiation matrix")

# `spectral_grid` samples with the same [i,j] -> (nodes[i], nodes[j]) convention, so the two
# grids hold the same numbers and the fields can be compared elementwise.
check("the two grids sample the same points",
    maximum(abs, g.x₁ .- [gp.nodes[i] for i in 1:N, j in 1:N]) < 1e-14 &&
        maximum(abs, g.x₂ .- [gp.nodes[j] for i in 1:N, j in 1:N]) < 1e-14,
    @sprintf("max |Δx| = %.2e",
        max(maximum(abs, g.x₁ .- [gp.nodes[i] for i in 1:N, j in 1:N]),
        maximum(abs, g.x₂ .- [gp.nodes[j] for i in 1:N, j in 1:N]))))

for (name, F) in (("u", U), ("h", H))
    e₁ = maximum(abs, ∂₁(g, F) .- ∂x(gp, F))
    e₂ = maximum(abs, ∂₂(g, F) .- ∂y(gp, F))
    check("∂₁$name: FFT vs matrix", e₁ < 1e-12, @sprintf("max abs diff %.2e", e₁))
    check("∂₂$name: FFT vs matrix", e₂ < 1e-12, @sprintf("max abs diff %.2e", e₂))
end

# The derivatives are exact, not merely consistent: these fields are band-limited well below
# Nyquist, so the spectral derivative is the analytic one.
let e = maximum(abs, ∂₁(g, U) .- torus_field(g,
        (a, b) -> -0.7sin(a) - 0.3sin(a - b) + 0.45cos(3a + b)))
    check("∂₁u is the analytic derivative", e < 1e-12, @sprintf("max abs diff %.2e", e))
end

let e = maximum(abs, canonical_bracket(g, U, V) .-
                     (∂x(gp, U) .* ∂y(gp, V) .- ∂y(gp, U) .* ∂x(gp, V)))
    check("canonical_bracket agrees with the package's", e < 1e-12,
        @sprintf("max abs diff %.2e", e))
end

# =============================================================================================
header("2. the Poisson solve of eq:Poisson-eq-periodic")

let ω = torus_field(g, (a, b) -> cos(a) + 0.5sin(2b) + 0.25cos(a + b)),
    exact = torus_field(g, (a, b) -> cos(a) + 0.5sin(2b) / 4 + 0.25cos(a + b) / 2)

    φ = poisson_periodic(g, ω)
    check("-Δφ = ω on eigenmodes", maximum(abs, φ .- exact) < 1e-13,
        @sprintf("max abs diff %.2e", maximum(abs, φ .- exact)))
    check("φ has zero mean", abs(mean_value(g, φ)) < 1e-14,
        @sprintf("%.2e", abs(mean_value(g, φ))))
    check("-Δ(poisson(ω)) = ω", maximum(abs, .-laplacian(g, φ) .- ω) < 1e-12,
        @sprintf("max abs diff %.2e", maximum(abs, .-laplacian(g, φ) .- ω)))
end

# The zero mode is discarded rather than raising: a field with a mean solves for its
# fluctuation, which is what `eq:Poisson-eq-periodic` intends.
let ω = torus_field(g, (a, b) -> cos(a) + 3.0)
    φ = poisson_periodic(g, ω)
    check("a nonzero mean is discarded, not amplified",
        maximum(abs, φ .- torus_field(g, (a, b) -> cos(a))) < 1e-13,
        @sprintf("max abs diff %.2e",
            maximum(abs, φ .- torus_field(g, (a, b) -> cos(a)))))
end

# =============================================================================================
header("3. eq:parallel-diffusion: [h,[h,u]] = div(X_h ⊗ X_h ∇u)")

let X = hamiltonian_field(g, H), ∂₁u = ∂₁(g, U), ∂₂u = ∂₂(g, U)

    # The divergence form, assembled component by component.
    F₁ = X[1] .* (X[1] .* ∂₁u .+ X[2] .* ∂₂u)
    F₂ = X[2] .* (X[1] .* ∂₁u .+ X[2] .* ∂₂u)
    div = ∂₁(g, F₁) .+ ∂₂(g, F₂)
    nested = double_bracket_field(g, U, H)
    e = maximum(abs, div .- nested) / maximum(abs, nested)
    check("nested bracket = divergence form", e < 1e-12, @sprintf("rel %.2e", e))

    check("∇·X_h = 0", maximum(abs, ∂₁(g, X[1]) .+ ∂₂(g, X[2])) < 1e-12,
        @sprintf("%.2e", maximum(abs, ∂₁(g, X[1]) .+ ∂₂(g, X[2]))))
    check("X_h · ∇h = 0 pointwise",
        maximum(abs, X[1] .* ∂₁(g, H) .+ X[2] .* ∂₂(g, H)) < 1e-13,
        @sprintf("%.2e", maximum(abs, X[1] .* ∂₁(g, H) .+ X[2] .* ∂₂(g, H))))
end

# =============================================================================================
header("4. the semi-discrete conservation laws")

# Double bracket, analytic test case: H(u) = (h - h_Ω, u) is conserved because the bracket is
# degenerate on it, and S = ‖ω‖²/2 decreases because it is positive semi-definite.
let ω = U .- mean_value(g, U), hΩ = mean_value(g, H), hz = H .- hΩ
    f = double_bracket_field(g, ω, H)
    dH = l2inner(g, hz, f)
    dS = l2inner(g, ω, f)
    scale = l2norm(g, hz) * l2norm(g, f)
    check("double bracket: dH/dt = 0", abs(dH) / scale < 1e-13,
        @sprintf("%.2e   (normalised)", abs(dH) / scale))
    check("double bracket: dS/dt < 0", dS < 0, @sprintf("dS/dt = %.6e", dS))
    # dS/dt = -∫|X_h·∇ω|², which is the entropy production written out.
    X = hamiltonian_field(g, H)
    prod = -integrate(g, (X[1] .* ∂₁(g, ω) .+ X[2] .* ∂₂(g, ω)) .^ 2)
    check("dS/dt = -∫|X_h·∇ω|²", abs(dS - prod) / abs(prod) < 1e-12,
        @sprintf("%.10e vs %.10e", dS, prod))
end

# Projector bracket, reduced Euler: H = ½(φ,ω) conserved, S monotone, and (S,S) matches
# eq:SS-projector.
let ω = U .- mean_value(g, U)
    φ = poisson_periodic(g, ω)
    f = projector_bracket_field(g, ω, φ)
    H₀ = l2inner(g, φ, ω) / 2
    dH = l2inner(g, φ, f)
    dS = l2inner(g, ω, f)
    scale = l2norm(g, φ) * l2norm(g, f)
    check("projector: dH/dt = 0", abs(dH) / scale < 1e-13,
        @sprintf("%.2e   (normalised)", abs(dH) / scale))
    check("projector: dS/dt < 0", dS < 0, @sprintf("dS/dt = %.6e", dS))

    # eq:SS-projector: (S,S) = 2S - 4H₀²/‖φ‖², and dS/dt = -(S,S).
    S = l2inner(g, ω, ω) / 2
    SS = 2S - 4H₀^2 / l2inner(g, φ, φ)
    check("(S,S) = 2S - 4H₀²/‖φ‖²  [eq:SS-projector]", abs(-dS - SS) / abs(SS) < 1e-12,
        @sprintf("%.12e vs %.12e   rel %.2e", -dS, SS, abs(-dS - SS) / abs(SS)))
end

# =============================================================================================
header("5. the state is mean-free and stays that way")

for name in ("a1", "a2", "a3", "a4")
    spec = SECTION4_RUNS[name]
    ω, uΩ = spectral_state(g, spec)
    check(@sprintf("%s: initial ω has zero mean", name), abs(mean_value(g, ω)) < 1e-14,
        @sprintf("mean = %.2e   u_Ω = %.8f", mean_value(g, ω), uΩ))

    ĥ = spec.h === nothing ? nothing : torus_field(g, spec.h)
    f = if spec.bracket === :double
        double_bracket_field(g, ω, ĥ === nothing ? poisson_periodic(g, ω) : ĥ)
    else
        projector_bracket_field(g, ω, poisson_periodic(g, ω))
    end
    check(@sprintf("%s: the vector field is mean-free", name),
        abs(mean_value(g, f)) / max(l2norm(g, f), 1e-300) < 1e-13,
        @sprintf("%.2e   (normalised)", abs(mean_value(g, f)) / max(l2norm(g, f), 1e-300)))
end

summary("verify_spectral.jl")
