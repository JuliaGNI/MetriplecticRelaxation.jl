#!/usr/bin/env julia
#
# The B-spline Galerkin layer, checked against the closed forms and against the spectral layer.
#
#     julia --project=scripts scripts/verify_spline.jl
#
# The claims here are the ones the Section 4 runs rest on, in the order a defect would surface:
#
#   1. the space integrates and projects what it should;
#   2. the bordered Poisson solve reproduces -Delta phi = omega, and the multiplier it carries
#      is at round-off (a nonzero multiplier would mean the right-hand side was not mean-free,
#      which would be a defect in the state convention rather than in the solve);
#   3. M*Lambda is symmetric, which is what EllipticEnergy asserts without checking, since
#      QuadraticHamiltonian's own check would assemble the dense matrix PoissonMap exists to
#      avoid;
#   4. the metric brackets are symmetric, positive semi-definite and degenerate on the energy
#      the flow conserves -- degeneracy_residual(flow, u), not the bracket's own two-argument
#      form, which is clean by construction;
#   5. A1's assembled fast path is the same operator as the generic vectorfield;
#   6. the spline and spectral discretisations agree on the vector field itself, which is the
#      cross-check the whole reproduction turns on, at the one place it can be made pointwise.

using MetriplecticRelaxation
using MetriplecticRelaxation: SplineTorus, PoissonMap, EllipticEnergy, LinearHamiltonian,
                              spline_state, spline_flow, spline_rhs, spline_grid,
                              fixed_double_operator, SpectralTorus, spectral_state,
                              torus_field, poisson_periodic, double_bracket_field,
                              projector_bracket_field, islands_h, SECTION4_RUNS,
                              integrate, l2inner, l2norm, mean_value,
                              initial_condition, periodise
using PoissonBrackets: nbasis, project, evaluate, mass_matrix, vectorfield, hamiltonian,
                       gradient, entropy, entropy_gradient, issymmetric,
                       ispositive_semidefinite, degeneracy_residual, domainvolume
using LinearAlgebra
using Printf
using Random
using SparseArrays

include(joinpath(@__DIR__, "check.jl"))
using .Checks: header, check, summary

Random.seed!(0x5c1e9a3b)

const t = SplineTorus(24, 3)
const N = nbasis(t.space)

# =============================================================================================
header("1. the space")

check("nbasis = n² for a periodic basis", N == 24^2, @sprintf("N = %d", N))
check("∫1 dx = |Ω| = 4π²", abs(integrate(t, ones(N)) - 4π^2) < 1e-12,
    @sprintf("%.14f vs %.14f", integrate(t, ones(N)), 4π^2))
check("domainvolume agrees", abs(domainvolume(t.space) - 4π^2) < 1e-12,
    @sprintf("%.14f", domainvolume(t.space)))

# L² projection reproduces a field already in the space: a trigonometric polynomial is not, so
# this is checked on the constant and on the projection being idempotent.
let û = project(t.space, x -> cos(x[1]) + 0.4sin(2x[2]) + 0.3cos(x[1] - x[2]))
    v̂ = project(t.space, x -> evaluate(t.space, û, (x[1], x[2])))
    check("L² projection is idempotent", norm(û - v̂) / norm(û) < 1e-12,
        @sprintf("rel %.2e", norm(û - v̂) / norm(û)))
end

let û = project(t.space, x -> 2.5)
    check("a constant projects onto the constant coefficient vector",
        norm(û .- 2.5) / 2.5 < 1e-12, @sprintf("rel %.2e", norm(û .- 2.5) / 2.5))
end

# =============================================================================================
header("2. the bordered Poisson solve")

for (name, f, λ) in (("cos x₁", x -> cos(x[1]), 1.0),
    ("sin x₂", x -> sin(x[2]), 1.0),
    ("cos(x₁+x₂)", x -> cos(x[1] + x[2]), 2.0),
    ("cos 2x₁", x -> cos(2x[1]), 4.0))
    ω̂ = project(t.space, f)
    φ̂ = t.Λ * ω̂
    # -Δ has these as eigenfunctions with eigenvalue λ, so φ = ω/λ. The Galerkin solve
    # reproduces that exactly, because the eigenfunction relation is inherited by the discrete
    # operators whenever the field is in the space -- which it is not, so this is a
    # discretisation-error check rather than an identity.
    e = norm(φ̂ .- ω̂ ./ λ) / norm(ω̂ ./ λ)
    check(@sprintf("-Δφ = ω for %-11s (λ = %g)", name, λ), e < 5e-4, @sprintf("rel %.2e",
        e))
end

let ω̂ = project(t.space, x -> cos(x[1]) + 0.4sin(2x[2]) + 0.3cos(x[1] - x[2]))
    φ̂ = t.Λ * ω̂
    check("φ has zero mean", abs(mean_value(t, φ̂)) < 1e-12,
        @sprintf("%.2e", abs(mean_value(t, φ̂))))
end

# The Lagrange multiplier of the bordered system is at round-off exactly when the right-hand
# side is mean-free. If the state convention ever stopped subtracting the mean, this is the
# check that would notice.
let ω̂ = project(t.space, x -> cos(x[1]) + 0.4sin(2x[2])), N = nbasis(t.space)
    rhs = vcat(t.M * ω̂, 0.0)
    sol = t.Λ.F \ rhs
    check("the bordered multiplier is at round-off", abs(sol[N + 1]) < 1e-12,
        @sprintf("α = %.2e", abs(sol[N + 1])))
end

# =============================================================================================
header("3. M Λ is symmetric, as EllipticEnergy assumes")

# Assembled densely on this small space only. On a run-sized space this is exactly the matrix
# PoissonMap exists not to form.
let MΛ = t.M * Matrix(t.Λ)
    asym = norm(MΛ - MΛ') / norm(MΛ)
    check("‖MΛ - (MΛ)ᵀ‖/‖MΛ‖ = 0", asym < 1e-11, @sprintf("%.2e", asym))
end

let û = randn(N), Hen = EllipticEnergy(t.Λ, t.M)
    û .-= mean_value(t, û)
    # H = ½(φ,ω) and the gradient is Mφ̂: check the gradient against a directional difference.
    v̂ = randn(N)
    v̂ .-= mean_value(t, v̂)
    ε = 1e-6
    fd = (hamiltonian(Hen, t.space, û .+ ε .* v̂) -
          hamiltonian(Hen, t.space, û .- ε .* v̂)) / 2ε
    an = dot(gradient(Hen, t.space, û), v̂)
    check("EllipticEnergy: gradient matches a central difference",
        abs(fd - an) / abs(an) < 1e-8, @sprintf("%.10f vs %.10f", fd, an))
    check("H = ½(φ,ω) > 0 for ω ≠ 0", hamiltonian(Hen, t.space, û) > 0,
        @sprintf("H = %.6f", hamiltonian(Hen, t.space, û)))
end

# =============================================================================================
header("4. the metric brackets, on the flow's own energy")

for name in ("a1", "a2", "a3", "a4")
    spec = SECTION4_RUNS[name]
    f = spline_flow(t, spec)
    ω̂, _ = spline_state(t, spec)
    # A state with no special structure, so that no property below holds by accident.
    ŵ = ω̂ .+ 0.1 .* randn(N)
    ŵ .-= mean_value(t, ŵ)

    check(@sprintf("%s: the bracket is symmetric", name), issymmetric(f.metric, ŵ), "")
    check(@sprintf("%s: the bracket is positive semi-definite", name),
        ispositive_semidefinite(f.metric, ŵ), "")
    r = degeneracy_residual(f, ŵ)
    check(@sprintf("%s: 𝔾 ∂H/∂û = 0 on the FLOW's Hamiltonian", name), r < 1e-11,
        @sprintf("%.2e", r))

    # The two conservation laws, semi-discretely.
    v = vectorfield(f, ŵ)
    dH = dot(gradient(f, ŵ), v)
    dS = dot(entropy_gradient(f, ŵ), v)
    scale = norm(gradient(f, ŵ)) * norm(v)
    check(@sprintf("%s: dH/dt = 0", name), abs(dH) / scale < 1e-11,
        @sprintf("%.2e   (normalised)", abs(dH) / scale))
    check(@sprintf("%s: dS/dt < 0", name), dS < 0, @sprintf("dS/dt = %.6e", dS))
end

# =============================================================================================
header("5. A1's assembled operator is the generic vector field")

let spec = SECTION4_RUNS["a1"]
    f = spline_flow(t, spec)
    A = fixed_double_operator(t, spec)
    ω̂, _ = spline_state(t, spec)
    ŵ = ω̂ .+ 0.05 .* randn(N)
    ŵ .-= mean_value(t, ŵ)

    fast = -(t.Mfac \ (A * ŵ))
    generic = vectorfield(f, ŵ)
    e = norm(fast - generic) / norm(generic)
    check("assembled path = vectorfield(MetriplecticFlow)", e < 1e-11,
        @sprintf("rel %.2e", e))
    check("the assembled operator is sparse", nnz(A) < N^2 / 4,
        @sprintf("nnz = %d of N² = %d  (%.1f%%)", nnz(A), N^2, 100nnz(A) / N^2))
    check("A is symmetric", norm(A - A') / norm(A) < 1e-12,
        @sprintf("%.2e", norm(A - A') / norm(A)))
    # A·e = 0: constants are in the kernel, which is why S may be written in u or in ω.
    check("A annihilates constants", norm(A * ones(N)) / norm(A) < 1e-12,
        @sprintf("%.2e", norm(A * ones(N)) / norm(A)))
end

# =============================================================================================
header("6. spline and spectral agree on the vector field")

# The two discretisations are compared where they can be compared pointwise: the right-hand side
# evaluated on a common initial condition, sampled on the spectral grid. They do not agree to
# round-off and must not be expected to -- one is 4th-order accurate and the other spectral --
# so this is a RELATIVE agreement at the level of the coarser method's own error, and section 7
# of `run_a1.jl` is what turns it into a refinement statement.
let Ns = 64, tf = SplineTorus(64, 3), g = SpectralTorus(Ns)
    for name in ("a1", "a2", "a3", "a4")
        spec = SECTION4_RUNS[name]
        ω̂, _ = spline_state(tf, spec)
        ωs, _ = spectral_state(g, spec)

        # The initial conditions themselves, on the common grid. A4 is the exception and it is
        # expected: its Gaussian sits π/2 from the periodic boundary with w₂ = 1 and amplitude
        # 1.8, so the formula as printed is discontinuous on T² by 0.153, and the two methods
        # resolve that jump differently. See `Gaussian`.
        ω̂g = spline_grid(tf, ω̂, Ns)
        e₀ = maximum(abs, ω̂g .- ωs) / maximum(abs, ωs)
        tol = name == "a4" ? 1e-1 : 5e-3
        check(@sprintf("%s: the two initial conditions agree", name), e₀ < tol,
            @sprintf("max rel %.2e   (tol %.0e)", e₀, tol))

        rhs = spline_rhs(tf, spec)
        vs = if spec.bracket === :double
            ĥ = spec.h === nothing ? poisson_periodic(g, ωs) : torus_field(g, spec.h)
            double_bracket_field(g, ωs, ĥ)
        else
            projector_bracket_field(g, ωs, poisson_periodic(g, ωs))
        end
        vg = spline_grid(tf, rhs(ω̂), Ns)
        e = maximum(abs, vg .- vs) / maximum(abs, vs)
        check(@sprintf("%s: the two vector fields agree", name), e < 1e-1,
            @sprintf("max rel %.2e", e))
    end
end

# =============================================================================================
header("7. A4's disagreement is the initial condition, not the solver")

# The control for section 6. If periodising the Gaussian — the only change — brings A4 back to
# the agreement A2 gets with an otherwise identical Gaussian, then the 7e-2 above is the stated
# initial condition's own discontinuity and not a defect in either discretisation. If it did
# NOT, the cause would be somewhere in the solvers and section 6 would be hiding it.
let Ns = 64, tf = SplineTorus(64, 3), g = SpectralTorus(Ns)
    a4 = SECTION4_RUNS["a4"]

    ω̂, _ = spline_state(tf, a4)
    ωs, _ = spectral_state(g, a4)
    e_literal = maximum(abs, spline_grid(tf, ω̂, Ns) .- ωs) / maximum(abs, ωs)

    a4p = periodise(a4)
    ω̂p, _ = spline_state(tf, a4p)
    ωsp, _ = spectral_state(g, a4p)
    e_periodic = maximum(abs, spline_grid(tf, ω̂p, Ns) .- ωsp) / maximum(abs, ωsp)

    check("periodising A4 restores spline/spectral agreement", e_periodic < 5e-3,
        @sprintf("literal %.2e  ->  periodised %.2e   (%.0fx better)",
            e_literal, e_periodic, e_literal / e_periodic))

    # And the jump itself, measured rather than asserted: the value the printed formula still
    # has at the far edge of the domain.
    jump = a4.gaussian(π, 2π)
    check("the jump across x₂ = 0 is 0.153", abs(jump - 0.153) < 5e-3,
        @sprintf("u_G(π, 2π) = %.6f   = %.1f%% of the peak %.3f",
            jump, 100jump / a4.gaussian(a4.gaussian.x₀...), a4.gaussian(a4.gaussian.x₀...)))

    # A2's Gaussian differs only in centre and amplitude, and is four orders of magnitude
    # smaller at the wrap — which is why only A4 is affected.
    a2 = SECTION4_RUNS["a2"]
    check("A2's Gaussian is negligible at the wrap", a2.gaussian(π, 0.0) < 1e-4,
        @sprintf("u_G(π, 0) = %.3e", a2.gaussian(π, 0.0)))

    # Periodising changes nothing measurable for A1-A3, which is why it is a choice only for A4.
    for nm in ("a1", "a2", "a3")
        s = SECTION4_RUNS[nm]
        d = maximum(abs,
            torus_field(g, initial_condition(s)) .-
            torus_field(g, initial_condition(periodise(s))))
        check(@sprintf("%s: periodising changes nothing", nm), d < 1e-4,
            @sprintf("max abs diff %.2e", d))
    end
end

summary("verify_spline.jl")
