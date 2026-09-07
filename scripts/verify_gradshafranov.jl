#!/usr/bin/env julia
#
# Verify `src/gradshafranov.jl`: the metriplectic discretisation of the Grad-Shafranov problem
# of Section 5.5, before any relaxation is run on it.
#
#     julia --project=scripts scripts/verify_gradshafranov.jl
#
# THE CLAIM under test is the one thing Section 5.5 adds to Section 5.4: the measure
# dmu = dr dz / r.  It has to appear in THREE places -- the inner moment quadrature of the
# collision bracket, the outer div_mu assembly, and the Delta-star operator -- and the checks
# below show, for each of the three separately, that a stray factor of r there is an O(1) error
# and not a tolerance-scale one.  A check that cannot tell the two readings apart is vacuous,
# so each control reports the factor by which it is wrong.
#
# Section 3 also settles the manuscript's own dx-versus-dmu inconsistency by MEASUREMENT rather
# than by reading: eq:entropy-2D writes S(u) = int s(x,u) dx while eq:GradShafranov-S-H writes
# int s(r,u) dmu, and the two differ by a factor of r.  The dmu reading reproduces eq:gs-ref;
# the dx reading does not, and misses it by a factor of 4.

using MetriplecticRelaxation
using MetriplecticRelaxation: GradShafranovBox, SECTION55_RUNS,
                              GS_RADIAL, GS_AXIAL, GS_LAMBDA_CONTINUUM,
                              HERRNEGGER_C, HERRNEGGER_D, gs_density,
                              herrnegger_mobility, gs_stiffness, gs_eigenvalue, gs_state,
                              gs_flow, gs_fit, gs_rayleigh, gs_current, gs_ordinate,
                              initial_condition, integrate,
                              Diagnostics, scatter_data
using PoissonBrackets: CollisionBracket, MetriplecticFlow, QuadraticHamiltonian, Integrator,
                       ImplicitMidpoint, integrate_step!, metric_operator, metric_matrix,
                       metric_apply, degeneracy_residual, ispositive_semidefinite,
                       issymmetric, entropy_production, entropy,
                       hamiltonian, vectorfield, mass_matrix, weighted_matrix,
                       quadrature_nodes, quadrature_weights, basis_values, field,
                       nbasis, default_f_abstol
using LinearAlgebra
using Random
using SparseArrays
using Printf

include(joinpath(@__DIR__, "check.jl"))
using .Checks: header, check, summary

const C, D = HERRNEGGER_C, HERRNEGGER_D
const spec = SECTION55_RUNS["c1"]

println("verify_gradshafranov.jl  --  the §5.5 discretisation and its measure")
@printf("    Ω = [%.1f, %.1f] × [%.1f, %.1f]   C = %.1f   D = %.1f   dμ = dr dz / r\n",
    GS_RADIAL..., GS_AXIAL..., C, D)

const box = GradShafranovBox((10, 12), 2)
const flow = gs_flow(box)
const ĵ₀ = gs_state(box, spec)
const s = box.space
const λh = gs_eigenvalue(box)

@printf("    %s   λ_h = %.10f   continuum %.10f\n", repr(box), λh, GS_LAMBDA_CONTINUUM)

# =================================================================================================
header("1. the initial condition is the manuscript's")

# `w²`, not `w`: §5 prints the SQUARED widths. `gaussian_w2` is the only place the square root
# is taken, and C1's values differ from B1's, so the transcription is re-checked here.
let u₀ = initial_condition(spec), (r₀, z₀) = (4.0, 0.0)
    check("u₀ peaks at 1 on the axis, N = 1", abs(u₀(r₀, z₀) - 1) < 1e-15,
        @sprintf("u₀(%.1f, %.1f) = %.16f", r₀, z₀, u₀(r₀, z₀)))
    for (k, w², shift) in ((1, 0.5, (sqrt(0.5), 0.0)), (2, 3.2, (0.0, sqrt(3.2))))
        v = u₀(r₀ + shift[1], z₀ + shift[2])
        check(
            @sprintf("the e-folding distance along axis %d is √(w²) = %.6f, not w² = %.1f",
                k, sqrt(w²), w²),
            abs(v - exp(-1)) < 1e-14,
            @sprintf("u₀ at the offset = %.16f   1/e = %.16f", v, exp(-1)))
    end
end

# The mesh criterion. The relaxed state is the lowest Δ* eigenmode, which the space resolves to
# 1e-5 at 6 cells, so what needs resolution is the initial Gaussian's narrower direction,
# w₁ = √0.5 = 0.707 on a domain of width 6. These rows are the measurement behind the recorded
# 18×21 rather than an argument for it; 22×26 is there to show the trend continues.
function gaussian_error(cells)
    b = GradShafranovBox(cells, 2)
    ĵ = gs_state(b, spec)
    xs, ws = quadrature_nodes(b.space), quadrature_weights(b.space)
    exact = [initial_condition(spec)(p[1], p[2]) for p in xs]
    return (b, sqrt(dot(ws, (gs_current(b, ĵ) .- exact) .^ 2) / dot(ws, exact .^ 2)))
end

const mesh_rows = [(cells, gaussian_error(cells)...)
                   for cells in ((10, 12), (14, 16), (18, 21), (22, 26))]

for (cells, b, e) in mesh_rows
    println(@sprintf("      %2d×%2d cells   N = %4d   ‖r j_h − u₀‖/‖u₀‖ = %.4e   λ_h = %.10f",
        cells..., nbasis(b), e, gs_eigenvalue(b)))
end

# `r · j_h` has to be the printed Gaussian: the state is j = u/r, so the division happens once
# in `gs_state` and multiplying back is the check that it happened exactly once. A factor of r
# out would show as an O(1) error rather than as the projection error the mesh rows measure.
let (_, _, e) = only(row for row in mesh_rows if row[1] == spec.cells)
    check("r · j_h recovers the printed Gaussian u₀ on the run's own mesh", e < 1e-2,
        @sprintf("‖r j_h − u₀‖/‖u₀‖ = %.4e at %d×%d cells", e, spec.cells...))
end

let es = [e for (_, _, e) in mesh_rows], ns = [c[1] for (c, _, _) in mesh_rows],
    q = log(es[1] / es[end]) / log(ns[end] / ns[1])

    check(
        "and the residual is a projection error, falling at the space's own order", q > 2.5,
        @sprintf("fitted order %.2f between 10 and 22 radial cells (degree 2 gives 3)", q))
end

# =================================================================================================
header("2. the state variable j = u/r is what makes the structure exact")

# Λ = (K^μ)⁻¹ M, hence Λᵀ K^μ = M, hence ∂H/∂ĵ = M ψ̂ — which is EXACTLY the condition the
# bracket's M⁻¹ 𝔸 M⁻¹ sandwich needs for the energy degeneracy. This identity is the whole
# reason the state is j and not u.
let Kμ = Matrix(gs_stiffness(s)), M = Matrix(mass_matrix(s)),
    e = maximum(abs, box.Λ' * Kμ - M) / maximum(abs, M)

    check("Λᵀ K^μ = M, so ∂H/∂ĵ = M ψ̂", e < 1e-10,
        @sprintf("max |ΛᵀK^μ − M| / max|M| = %.3e", e))
end

let e = maximum(abs, box.MΛ - box.MΛ') / maximum(abs, box.MΛ)
    check("M Λ is symmetric, so H is a quadratic form", e < 1e-14,
        @sprintf("max |MΛ − (MΛ)ᵀ| / max|MΛ| = %.3e", e))
end

# H = ½ ĵᵀMΛĵ has to BE the poloidal magnetic energy ½∫|∇ψ|²dμ and not a surrogate that happens
# to be conserved. ½ψ̂ᵀK^μψ̂ is that energy by the definition of `gs_stiffness`.
let ψ̂ = box.Λ * ĵ₀, H = hamiltonian(flow, ĵ₀), Hψ = dot(ψ̂, gs_stiffness(s), ψ̂) / 2,
    e = abs(H - Hψ) / abs(Hψ)

    check("H IS ½∫|∇ψ|² dμ, the poloidal magnetic energy", e < 1e-12,
        @sprintf("½ĵᵀMΛĵ = %.14e   ½ψ̂ᵀK^μψ̂ = %.14e   relative %.3e", H, Hψ, e))
end

# S = ½ ĵᵀWĵ has to be ∫ u²/2(Cr²+D) dμ, the manuscript's own entropy in its own variable. Both
# sides use the same quadrature, so this is exact rather than approximate.
let u = gs_current(box, ĵ₀), wμ = quadrature_weights(s) .* gs_density.(quadrature_nodes(s)),
    Sq = dot(wμ, u .^ 2 ./ (2 .* (C .* [p[1] for p in quadrature_nodes(s)] .^ 2 .+ D))),
    S = entropy(flow, ĵ₀), e = abs(S - Sq) / abs(Sq)

    check("S IS ∫ u²/2(Cr²+D) dμ", e < 1e-13,
        @sprintf("½ĵᵀWĵ = %.14e   quadrature %.14e   relative %.3e", S, Sq, e))
end

# The mass Casimir is the μ-integral of u, which is the dx-integral of j — the same cancellation
# the state variable is chosen for.
let x = quadrature_nodes(s), wμ = quadrature_weights(s) .* gs_density.(x),
    a = integrate(box, ĵ₀), b = dot(wμ, gs_current(box, ĵ₀)), e = abs(a - b) / abs(b)

    check("∫ j dx = ∫ u dμ", e < 1e-13,
        @sprintf("b·ĵ = %.14e   ∫u dμ = %.14e   relative %.3e", a, b, e))
end

# =================================================================================================
header("3. CONTROLS on Δ*: the u formulation and the plain Laplacian both break the degeneracy")

# The degeneracy 𝔾 ∂H/∂ĵ = 0 is what makes H conserved by the BRACKET rather than by the
# integrator, and it holds here to round-off. `degeneracy_residual(flow, ·)` takes the gradient
# from the flow, not from the bracket, which is the form that can fail.
for (name, ĵ) in (("the initial state", ĵ₀), (
    "a random state", randn(MersenneTwister(20260907),
        nbasis(box))))
    e = degeneracy_residual(flow, ĵ)
    check("𝔾 ∂H/∂ĵ = 0 at $(name)", e < 1e-12, @sprintf("normalised residual %.3e", e))
end

# CONTROL 1 -- the u formulation, i.e. the manuscript's own state variable with the manuscript's
# own μ-pairing. Then Λ_u = (K^μ)⁻¹M^μ, so ∂H/∂û = M^μψ̂ and the sandwich's plain M⁻¹ no longer
# recovers the bracket's generating field. This is the measurement that says the choice of state
# variable is FORCED and not presentational.
let x = quadrature_nodes(s),
    Mμ = Matrix(weighted_matrix(s, gs_density.(x), (0, 0), (0, 0))),
    Kμ = Matrix(gs_stiffness(s)), Λu = Symmetric(Kμ) \ Mμ, A = Mμ * Λu,
    Gu = CollisionBracket(s, Λu; mobility = (p, y) -> herrnegger_mobility(p[1]),
        mobility_derivative = 0, density = gs_density),
    fu = MetriplecticFlow(s, Gu, QuadraticHamiltonian(Matrix((A .+ A') ./ 2)),
        QuadraticHamiltonian(box.W)), e = degeneracy_residual(fu, ĵ₀)

    check("CONTROL: the u formulation loses the degeneracy at O(1)", e > 1e-3,
        @sprintf("normalised residual %.3e against %.3e for the j formulation — %.0fx",
            e, degeneracy_residual(flow, ĵ₀), e / degeneracy_residual(flow, ĵ₀)))
end

# CONTROL 2 -- the measure dropped from Δ* alone: the elliptic solve becomes the ordinary
# Dirichlet Laplacian while the energy stays the poloidal magnetic energy. Λᵀ K^μ ≠ M, and the
# degeneracy goes with it.
let K = Matrix(weighted_matrix(s, ones(length(quadrature_weights(s))), (1, 0), (1, 0))) .+
        Matrix(weighted_matrix(s, ones(length(quadrature_weights(s))), (0, 1), (0, 1))),
    M = Matrix(mass_matrix(s)), Λp = Symmetric(K) \ M, Kμ = Matrix(gs_stiffness(s)),
    A = Λp' * Kμ * Λp,
    Gp = CollisionBracket(s, Λp; mobility = (p, y) -> herrnegger_mobility(p[1]),
        mobility_derivative = 0, density = gs_density),
    fp = MetriplecticFlow(s, Gp, QuadraticHamiltonian(Matrix((A .+ A') ./ 2)),
        QuadraticHamiltonian(box.W)), e = degeneracy_residual(fp, ĵ₀),
    mismatch = maximum(abs, Λp' * Kμ - M) / maximum(abs, M)

    check(
        "CONTROL: a plain Laplacian in place of Δ* loses the degeneracy at O(1)", e > 1e-3,
        @sprintf("normalised residual %.3e   and ΛᵀK^μ misses M by %.3f relative",
            e, mismatch))
end

# =================================================================================================
header("4. CONTROLS on the bracket's measure: brute force against the collapsed form")

# `metric_operator` evaluates the manuscript's double integral through 14 global moments. Here
# the SAME operator is formed as a literal O(N_q²) double sum, once per measure combination, so
# that the measure in the inner integral and the measure in the outer one can be wrong
# separately. The grid is deliberately tiny: this is quadratic in the number of quadrature nodes.
function brute_operator(space, Λ, ĵ, inner_μ::Bool, outer_μ::Bool)
    x = quadrature_nodes(space)
    w = quadrature_weights(space)
    ρ = gs_density.(x)
    Mmob = [herrnegger_mobility(p[1]) for p in x]
    wi = inner_μ ? w .* ρ : w
    wo = outer_μ ? w .* ρ : w
    P = (basis_values(space, (1, 0)), basis_values(space, (0, 1)))
    φ̂ = Λ * ĵ
    # β = (∇ψ)^⊥, the perp that carries the degeneracy — fixed here as it is in the package.
    β = (-field(space, φ̂, (0, 1)), field(space, φ̂, (1, 0)))
    N, Q = nbasis(space), length(w)
    # 𝔸_KL = ∫∫ κ (∇Φ_K(x) − ∇Φ_K(x'))ᵀ Q₂(β(x) − β(x')) (∇Φ_L(x) − ∇Φ_L(x')) dμ' dμ, halved,
    # with Q₂(z) = z^⊥ ⊗ z^⊥ in two dimensions and κ = M(x)M(x').
    A = zeros(N, N)
    G = [Matrix(P[k]) for k in 1:2]
    for q in 1:Q, qq in 1:Q

        c = Mmob[q] * Mmob[qq] * wo[q] * wi[qq] / 2
        iszero(c) && continue
        # `β` is ALREADY the perp of ∇ψ, so `d = β(x) − β(x')` is the perp of the argument of
        # Q₂ and Q₂(z) = z^⊥ ⊗ z^⊥ = d ⊗ d. Do NOT perp `d` again: a second perp rotates the
        # rank-one factor by π/2 and the assembled bracket stops matching `CollisionBracket`
        # altogether, at an 0.85 relative disagreement rather than a small one.
        d = (β[1][q] - β[1][qq], β[2][q] - β[2][qq])
        # (d ⊗ d) contracted with the two difference gradients: one scalar per basis index.
        v = [d[1] * (G[1][K, q] - G[1][K, qq]) + d[2] * (G[2][K, q] - G[2][K, qq])
             for K in 1:N]
        A .+= c .* (v * v')
    end
    return A
end

let b = GradShafranovBox((4, 5), 2), ĵ = gs_state(b, spec),
    G = CollisionBracket(b.space, b.Λ; mobility = (p, y) -> herrnegger_mobility(p[1]),
        mobility_derivative = 0, density = gs_density), A = Matrix(metric_operator(G, ĵ)),
    ref = brute_operator(b.space, b.Λ, ĵ, true, true)

    println(@sprintf("      brute force on 4×5 cells: %d quadrature nodes, %d basis functions",
        length(quadrature_weights(b.space)), nbasis(b.space)))
    let e = maximum(abs, A - ref) / maximum(abs, ref)
        check(
            "the collapsed operator reproduces the μ-weighted double integral", e < 1e-11,
            @sprintf("max |𝔸 − 𝔸_brute| / max|𝔸_brute| = %.3e", e))
    end

    # The two halves of the measure, wrong one at a time -- `dμ(x')` in the inner moment
    # quadrature and `dμ(x)` in the outer div_μ assembly. Both are O(1), and they come out
    # IDENTICAL, which is worth stating rather than hiding: the integrand is symmetric under
    # x ↔ x' -- the difference gradient is antisymmetric and appears twice -- so exchanging
    # which of the two measures is wrong renames the integration variables and nothing more.
    # A stray factor of r cannot be localised to one of the two integrals by this test; it can
    # only be detected, which is what the test is for.
    for (name, inner, outer) in (("the inner moment quadrature", false, true),
        ("the outer div_μ assembly", true, false),
        ("both", false, false))
        wrong = brute_operator(b.space, b.Λ, ĵ, inner, outer)
        e = maximum(abs, wrong - ref) / maximum(abs, ref)
        check("CONTROL: dropping 1/r from $(name) is O(1)", e > 0.1,
            @sprintf("max |𝔸_wrong − 𝔸| / max|𝔸| = %.3f", e))
    end

    # `metric_apply` evaluates the same operator through the 𝔽_s moment route rather than
    # through the ℝ,𝕊,𝕋 factorisation `metric_operator` uses, so this is a second independent
    # path through the collapse -- with the measure -- and not a tautology.
    let v = randn(MersenneTwister(3), nbasis(b.space)), Ga = metric_apply(G, ĵ, v),
        Gm = metric_matrix(G, ĵ) * v, e = maximum(abs, Ga - Gm) / maximum(abs, Gm)

        check("metric_apply's 𝔽_s route agrees with the assembled operator", e < 1e-10,
            @sprintf("max |𝔾v − 𝔸-route| / max|𝔾v| = %.3e", e))
    end

    # And the collapsed form with `density = 1`, which is the same error reached through the
    # package's own knob rather than through the brute-force reference.
    let G1 = CollisionBracket(b.space, b.Λ; mobility = (p, y) -> herrnegger_mobility(p[1]),
            mobility_derivative = 0, density = 1.0),
        e = maximum(abs, Matrix(metric_operator(G1, ĵ)) - ref) / maximum(abs, ref)

        check("CONTROL: `density = 1` in the bracket is O(1)", e > 0.1,
            @sprintf("max |𝔸 − 𝔸_μ| / max|𝔸_μ| = %.3f", e))
    end
end

# =================================================================================================
header("5. CONTROL on the entropy's measure: eq:entropy-2D's dx is a typo for dμ")

# The manuscript writes S(u) = ∫ s(x,u) dx at eq:entropy-2D and S(u) = ∫ s(r,u) dμ at
# eq:GradShafranov-S-H, and asserts δS/δu = ∂_y s in both. Under the μ-pairing those disagree by
# a factor of r. Which one eq:gs-ref belongs to is settled here by measurement.
#
# In the j variable the two readings are two entropy weights: σ = r/(Cr²+D) for dμ and
# σ_dx = r²/(Cr²+D) for dx. Each has its own relaxed state -- the lowest eigenvector of
# (W, MΛ) -- and the question is which of the two satisfies u/(Cr²+D) = λψ.
let x = quadrature_nodes(s),
    Wdx = Matrix(weighted_matrix(s,
        [p[1]^2 / (C * p[1]^2 + D) for p in x], (0, 0), (0, 0)))

    function relaxed(W)
        F = eigen(Symmetric(box.MΛ), Symmetric((W .+ W') ./ 2))
        (λ = inv(maximum(real, F.values)), ĵ = F.vectors[:, argmax(real(F.values))])
    end

    a, b = relaxed(box.W), relaxed(Wdx)
    (λa, _, rela) = gs_fit(box, a.ĵ)
    (λb, _, relb) = gs_fit(box, b.ĵ)

    check("the dμ reading's relaxed state satisfies eq:gs-ref", rela < 5e-3,
        @sprintf("‖σj − λψ‖/‖σj‖ = %.4e   fitted λ = %.10f   λ_h = %.10f",
            rela, λa, λh))
    check("CONTROL: the dx reading's relaxed state does NOT — an O(1) miss", relb > 0.1,
        @sprintf("‖σj − λψ‖/‖σj‖ = %.4f, i.e. %.0f× the dμ reading's; its own λ = %.8f",
            relb, relb / rela, λb))
    check("CONTROL: and its eigenvalue is O(1) away from the reference",
        abs(b.λ - GS_LAMBDA_CONTINUUM) / GS_LAMBDA_CONTINUUM > 0.1,
        @sprintf("dx reading λ = %.8f   dμ reading %.10f   continuum %.10f",
            b.λ, a.λ, GS_LAMBDA_CONTINUUM))

    # The relaxed state of the dμ reading satisfies eq:gs-ref only up to the L² projection of
    # σj: the eigenvector satisfies Π(σj) = λψ exactly, and `gs_fit` measures the pointwise
    # residual. So the number above must FALL with the mesh, which is what says the remaining
    # disagreement is a projection error and not a formulation error.
    for cells in ((6, 7), (10, 12), (14, 16), (18, 21))
        bb = GradShafranovBox(cells, 2)
        F = eigen(Symmetric(bb.MΛ), Symmetric(bb.W))
        v = F.vectors[:, argmax(real(F.values))]
        println(@sprintf("      %2d×%2d cells   ‖σj − λψ‖/‖σj‖ = %.4e   λ_h = %.10f",
            cells..., gs_fit(bb, v)[3], gs_eigenvalue(bb)))
    end
end

# =================================================================================================
header("6. the bracket is symmetric, positive semi-definite and dissipative")

let rng = MersenneTwister(20260907)
    for (name, ĵ) in (("the initial state", ĵ₀), (
        "a random state", randn(rng, nbasis(box))))
        check("𝔾 is symmetric at $(name)", issymmetric(flow.metric, ĵ), "")
        check("𝔾 is positive semi-definite at $(name)",
            ispositive_semidefinite(flow.metric, ĵ), "")
    end
end

check("(S,S) > 0 at the initial state", entropy_production(flow, ĵ₀) > 0,
    @sprintf("(S,S) = %.6e", entropy_production(flow, ĵ₀)))

# The mobility M = Cr²+D is what makes κ = M(x)M(x') a POSITIVE kernel, and eq:M-condition
# requires it. A sign-changing mobility is a bracket that is not positive semi-definite, which is
# the failure `ispositive_semidefinite` exists to name.
let Gneg = CollisionBracket(s, box.Λ; mobility = (p, y) -> p[1] - 4.0,
        mobility_derivative = 0, density = gs_density),
    λs = eigvals(Symmetric(Matrix(metric_matrix(Gneg, ĵ₀))))

    check("CONTROL: a sign-changing mobility breaks positive semi-definiteness",
        !ispositive_semidefinite(Gneg, ĵ₀),
        @sprintf("λ_min = %.4e   λ_max = %.4e   ratio %.3e",
            minimum(λs), maximum(λs), minimum(λs) / maximum(λs)))
end

# =================================================================================================
header("7. the Poincaré floor S/H ≥ λ_h holds for every admissible state")

# `gs_rayleigh` is the Rayleigh quotient of the constrained minimisation, so its infimum over the
# whole space is λ_h -- not merely along the flow. That makes it a precondition to assert at
# t = 0, and the control is that no state can get under it.
let rng = MersenneTwister(551), worst = Inf, arg = 0
    for k in 1:400
        q = gs_rayleigh(box, randn(rng, nbasis(box)))
        q < worst && (worst = q; arg = k)
    end
    check("no random state gets below λ_h", worst ≥ λh * (1 - 1e-12),
        @sprintf("smallest S/H over 400 draws = %.10f at draw %d   λ_h = %.10f   excess %+.3e",
            worst, arg, λh, worst - λh))
end

let q = gs_rayleigh(box, ĵ₀)
    check("the initial state is above the floor", q > λh,
        @sprintf("S₀/H₀ = %.10f   λ_h = %.10f   excess %+.4e", q, λh, q - λh))
end

# The two λ_h routes: (MΛ, W), which is the invariants the flow carries, and (K^μ, B), which is
# eq:Grad-Shafranov-equation tested against the basis. They are DIFFERENT discretisations of the
# same continuum eigenvalue -- the first projects σj, the second does not -- so they agree at the
# order of the space and not exactly, and the gap must fall with the mesh.
for cells in ((6, 7), (10, 12), (14, 16), (18, 21))
    b = GradShafranovBox(cells, 2)
    λa, λb = gs_eigenvalue(b), gs_eigenvalue(b.space)
    println(@sprintf("      %2d×%2d cells   (MΛ,W) %.12f   (K^μ,B) %.12f   relative %.3e",
        cells..., λa, λb, abs(λa - λb) / λb))
end

let b = GradShafranovBox((18, 21), 2),
    e = abs(gs_eigenvalue(b) - gs_eigenvalue(b.space)) / gs_eigenvalue(b.space)

    check("the two λ_h routes agree at the order of the space", e < 1e-7,
        @sprintf("relative gap %.3e at 18×21 cells", e))
end

let b = GradShafranovBox((18, 21), 2),
    e = abs(gs_eigenvalue(b) - GS_LAMBDA_CONTINUUM) / GS_LAMBDA_CONTINUUM

    check("and both are within 1e-5 of the continuum eigenvalue", e < 1e-5,
        @sprintf("λ_h = %.12f   continuum %.10f   relative %+.3e",
            gs_eigenvalue(b), GS_LAMBDA_CONTINUUM, e))
end

# =================================================================================================
header("8. the residual is a cubic polynomial and the Jacobian is analytic")

# M = Cr²+D does not depend on the state, so 𝔾 is quadratic in ĵ and ∂S/∂ĵ is linear: the vector
# field is an exact cubic. Third differences of a cubic in a fixed direction are constant, which
# is a statement about the whole discretisation and not about one evaluation.
let rng = MersenneTwister(7), v = randn(rng, nbasis(box)), d = randn(rng, nbasis(box)),
    g(t) = vectorfield(flow, v .+ t .* d),
    third(t, h) = (g(t + 2h) .- 3 .* g(t + h) .+ 3 .* g(t) .- g(t - h)) ./ h^3,
    a = third(0.0, 0.25), c = third(0.7, 0.25), e = maximum(abs, a - c) / maximum(abs, a)

    check("the vector field is an exact cubic in the degrees of freedom", e < 1e-8,
        @sprintf("third differences at t = 0 and t = 0.7 agree to %.3e", e))
end

# =================================================================================================
header("9. a short run conserves H and dissipates S")

# The long run is `run_c1.jl`; five steps here are enough to say the structure survives the time
# discretisation, and they are what would catch a Jacobian error before an hour is spent.
let ĵ = copy(ĵ₀), H₀ = hamiltonian(flow, ĵ₀), S = Float64[entropy(flow, ĵ₀)],
    integ = Integrator(flow, ImplicitMidpoint(), spec.Δt; û₀ = ĵ₀), drift = 0.0

    for _ in 1:5
        integrate_step!(ĵ, integ)
        push!(S, entropy(flow, ĵ))
        drift = max(drift, abs(hamiltonian(flow, ĵ) - H₀))
    end
    tol = default_f_abstol(Float64, nbasis(box), ĵ₀)
    check("H is conserved at the Newton residual tolerance over 5 steps",
        drift / abs(H₀) < 1e-11,
        @sprintf("max |ΔH| = %.3e   |ΔH|/H₀ = %.3e   H₀ = %.12e   f_abstol/H₀ = %.2e",
            drift, drift / abs(H₀), H₀, tol / abs(H₀)))
    check("S decreases at every step", all(diff(S) .< 0),
        @sprintf("S: %.10e → %.10e   worst increment %+.3e",
            S[1], S[end], maximum(diff(S))))
    check("S stays above the floor λ_h H₀", minimum(S) > λh * abs(H₀),
        @sprintf("min S − λ_h H₀ = %+.4e   λ_h H₀ = %.10e", minimum(S) - λh * H₀, λh * H₀))
end

# The scatter data is what the figure and the reference check are drawn from, so its two columns
# have to be ψ and u/(Cr²+D) and not the state.
let d = Diagnostics(box, spec), (_, y) = scatter_data(d, ĵ₀),
    e = maximum(abs, y - gs_ordinate(box, ĵ₀))

    check("scatter_data returns (ψ, u/(Cr²+D)) and not (ψ, j)",
        e == 0 &&
            maximum(abs,
                y - field(s, ĵ₀, (0, 0))) > 1e-3,
        @sprintf("ordinate matches gs_ordinate exactly; it differs from j by %.4e",
            maximum(abs, y - field(s, ĵ₀, (0, 0)))))
end

# =================================================================================================
header("10. CONTROL: the state space, and why eq:gs-ref needs the Dirichlet one")

# `eq:gs-ref` is the μ = 0, c = 0 member of the equilibrium family δS/δj = λψ + c·x + μ, whose
# extra multipliers belong to the bracket's mass and momentum Casimirs. Only a space that does
# not contain 1, x₁, x₂ forces them to vanish. This is §5.4's finding on §5.5's problem.
#
# IT HAS TO BE MEASURED ALONG THE FLOW, and an eigenvalue problem cannot stand in for it. The
# pencil (MΛ, W) carries no mass constraint, so its lowest eigenvector is the μ = 0, c = 0
# member in BOTH spaces — measured, the two spaces' eigenvalues agree to 5e-9 and both
# eigenvectors fit eq:gs-ref to 4e-4, so an eigen-based control here would report NOTHING. What
# separates the spaces is that the FLOW in the free space conserves the mass and the momenta and
# therefore cannot reach that member from an initial state whose mass is nonzero. A short
# relaxation in each space is the control, exactly as §5.4's B3 needed one.
function multiplier_fit(b, ĵ)
    sp = b.space
    wμ = quadrature_weights(sp) .* gs_density.(quadrature_nodes(sp))
    y = gs_ordinate(b, ĵ)
    ψ = field(sp, b.Λ * ĵ, (0, 0))
    x = quadrature_nodes(sp)
    # The four members of the equilibrium family: ψ, and the derivatives 1, r, z of the mass and
    # momentum Casimirs. `extra` is the largest of the three multipliers that eq:gs-ref sets to
    # zero, as a fraction of the ordinate's own size.
    B = [ψ ones(length(ψ)) [p[1] for p in x] [p[2] for p in x]]
    θ = (B' * (wμ .* B)) \ (B' * (wμ .* y))
    scale = sqrt(max(dot(wμ, y .^ 2) / sum(wμ), 0.0))
    return (λfit = θ[1], extra = maximum(abs, θ[2:4]) / scale)
end

let rows = NamedTuple[], steps = 240
    for st in (:dirichlet, :free)
        b = GradShafranovBox((10, 12), 2; state = st)
        f = gs_flow(b)
        ĵ = gs_state(b, spec)
        m₀ = integrate(b, ĵ)
        integ = Integrator(f, ImplicitMidpoint(), spec.Δt; û₀ = copy(ĵ))
        for _ in 1:steps
            integrate_step!(ĵ, integ)
        end
        push!(rows,
            (; state = st, N = nbasis(b), λ = gs_eigenvalue(b), rel = gs_fit(b, ĵ)[3],
                multiplier_fit(b, ĵ)..., mass = (integrate(b, ĵ) - m₀) / m₀))
    end
    println(@sprintf("      %d steps at Δt = %.4g, i.e. to t = %.4g, on 10×12 cells",
        steps, spec.Δt, steps * spec.Δt))
    for r in rows
        println(@sprintf("      %-10s N = %4d   ‖σj − λψ‖/‖σj‖ = %.4e   |extra|/‖σj‖ = %.4e   Δmass %+.3e",
            ":" * String(r.state), r.N, r.rel, r.extra, r.mass))
    end
    d, f = rows[1], rows[2]
    check("the mass is conserved in the free space and drifts in the Dirichlet one",
        abs(f.mass) < 1e-10 && abs(d.mass) > 0.1,
        @sprintf(":free %+.3e   :dirichlet %+.3e — the constants are in one space and not the other",
            f.mass, d.mass))
    check("CONTROL: the free space's relaxed state misses eq:gs-ref by orders more",
        f.rel / d.rel > 50,
        @sprintf(":free %.4e against :dirichlet %.4e — %.0f×", f.rel, d.rel, f.rel / d.rel))
    check("and it carries the extra multipliers the Dirichlet space cannot",
        f.extra > 50 * d.extra,
        @sprintf("max|(μ, c₁, c₂)|/‖σj‖: :free %.4e   :dirichlet %.4e — %.0f×",
            f.extra, d.extra, f.extra / d.extra))
    check(
        "both spaces have the SAME eigenvalue, so λ_h is not what separates them",
        abs(f.λ - d.λ) / d.λ < 1e-8,
        @sprintf("λ_h: :free %.12f   :dirichlet %.12f   relative %+.3e",
            f.λ, d.λ, (f.λ - d.λ) / d.λ))
end

summary("verify_gradshafranov.jl")
