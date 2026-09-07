#!/usr/bin/env julia
#
# Everything the Section 5.4 runs rest on, settled before any of them is believed.
#
#     julia --project=scripts scripts/verify_euler.jl
#
# Six sections, in the order a defect in one would invalidate the next:
#
#   1. the initial conditions, and specifically that Section 5 states w^2 AND NOT w.  That is
#      the single most likely transcription error in the whole of Section 5 -- it is a factor of
#      ten in the width of the narrow direction -- so it is checked against the e-folding
#      distance, with a control that a `w` reading must fail;
#   2. the homogeneous-Dirichlet space: the domain it integrates, its first eigenvalue, and the
#      fact that CONSTANTS ARE NOT IN IT, which is what makes both of Section 5.4's closed-form
#      references exact -- and then, in 2b, the two separate things that decide whether B3's
#      initial condition is admissible at all: the mesh, and the floor;
#   3. the collision-like bracket on that space: symmetry, positive semi-definiteness, the
#      degeneracy (F,H) = 0, and the two independent evaluations of the same operator -- plus
#      two controls that MUST fail;
#   4. the two entropies, their analytic gradients and Hessians against central differences, and
#      both functionals against an independent quadrature;
#   5. the closed-form references S_eta = lambda_11 H_0 and lambda = (M+S)/2H_0;
#   6. Crank-Nicolson on a short run: H to machine precision, S monotone, and the step-halving
#      difference that says the time step is a cost choice rather than an accuracy one.
#
# Everything here runs on spaces of at most a few hundred degrees of freedom, so the whole script
# is seconds. The claims that need resolution to be true at all -- the relaxed states themselves
# -- are the drivers' business.

using MetriplecticRelaxation
using MetriplecticRelaxation: SECTION5_RUNS, SECTION5_ORDER, EulerSpec, EulerSquare,
                              GibbsEntropy,
                              euler_state, euler_flow, euler_entropy_floor, eigenmode_fit,
                              gibbs_lambda, gibbs_fit, gibbs_residual, interior_weights,
                              state_extrema, dirichlet_eigenvalue,
                              perturbation_b2, Diagnostics, Trace, record!,
                              entropy_monotone, energy_error, l2inner, l2norm, integrate,
                              DIRICHLET_EIGENVALUE, SQUARE_LENGTH, B3_FLOOR
using PoissonBrackets
using PoissonBrackets: CollisionBracket, MetriplecticFlow, QuadraticHamiltonian,
                       Integrator, ImplicitMidpoint, integrate_step!,
                       nbasis, project, evaluate, field,
                       mass_matrix, stiffness_matrix, domainvolume,
                       quadrature_weights,
                       hamiltonian, gradient, hessian, vectorfield,
                       entropy_gradient, entropy_production,
                       issymmetric, ispositive_semidefinite, degeneracy_residual,
                       metric_matrix, metric_apply
using LinearAlgebra
using Printf
using Random
using SparseArrays

include(joinpath(@__DIR__, "check.jl"))
using .Checks: header, check, check_refined, summary

Random.seed!(0x5eb1a704)

println("verify_euler.jl  --  the Section 5.4 problem, space, bracket and references")

# =============================================================================================
header("1. Section 5 states w², not w")

for name in SECTION5_ORDER
    spec = SECTION5_RUNS[name]
    g = spec.gaussian
    w² = (g.w[1]^2, g.w[2]^2)
    check(@sprintf("%-3s w² = (0.01, 0.07) as printed", name),
        isapprox(w²[1], 0.01; rtol = 1e-14) && isapprox(w²[2], 0.07; rtol = 1e-14),
        @sprintf("w² = (%.16f, %.16f)   w = (%.10f, %.10f)", w²[1], w²[2], g.w[1], g.w[2]))
end

# The e-folding distance is the whole content of the distinction: ω₀ falls by exactly one factor
# of e at |x - x₀| = w, i.e. at √(w²) = 0.1 along the narrow axis. A `w = 0.01` reading would put
# it at 0.01, ten times closer, and would make the vortex ten times narrower than the manuscript
# says.
let g = SECTION5_RUNS["b1"].gaussian, x₀ = g.x₀
    peak = g(x₀[1], x₀[2])
    check("ω₀(x₀) = 1/N exactly", peak == 1.0, @sprintf("%.16f", peak))
    for (k, w²) in ((1, 0.01), (2, 0.07))
        d = ntuple(j -> j == k ? sqrt(w²) : 0.0, 2)
        e = g(x₀[1] + d[1], x₀[2] + d[2]) / peak
        check(@sprintf("axis %d: ω₀ falls by e at √(w²) = %.10f", k, sqrt(w²)),
            isapprox(e, exp(-1); rtol = 1e-14), @sprintf("ratio = %.16f   1/e = %.16f", e,
                exp(-1)))
        # The control: reading the printed number as `w` itself must be wrong by order one.
        d² = ntuple(j -> j == k ? w² : 0.0, 2)
        ew = g(x₀[1] + d²[1], x₀[2] + d²[2]) / peak
        check(@sprintf("axis %d: the `w = w²` reading is WRONG by order one", k),
            abs(ew - exp(-1)) > 0.1, @sprintf("ratio at %.4f = %.10f   (vs 1/e = %.6f)",
                w², ew, exp(-1)))
    end
end

# The three amplitudes, which is the other easy transcription error: B2 gives N and B3 gives 1/N.
for (name, amplitude) in (("b1", 1.0), ("b2", 0.01), ("b3", 10.0))
    g = SECTION5_RUNS[name].gaussian
    check(@sprintf("%-3s peak amplitude = %.4g", name, amplitude),
        isapprox(g(g.x₀[1], g.x₀[2]), amplitude; rtol = 1e-14),
        @sprintf("ω_G(x₀) = %.10f", g(g.x₀[1], g.x₀[2])))
end

# B2's added mode is an eigenfunction of -Δ on the square, which is what makes B2 a PERTURBED
# EQUILIBRIUM rather than an arbitrary second initial condition.
let λ = 52π^2, ε = 1e-5
    f(x₁, x₂) = perturbation_b2(x₁, x₂)
    x₁, x₂ = 0.31, 0.57
    lap = (f(x₁ + ε, x₂) - 2f(x₁, x₂) + f(x₁ - ε, x₂)) / ε^2 +
          (f(x₁, x₂ + ε) - 2f(x₁, x₂) + f(x₁, x₂ - ε)) / ε^2
    check("B2's mode is an eigenfunction of -Δ with λ = 52π²",
        isapprox(-lap / f(x₁, x₂), λ; rtol = 1e-6),
        @sprintf("-Δf/f = %.6f   52π² = %.6f", -lap / f(x₁, x₂), λ))
    check("B2's mode vanishes on ∂Ω",
        maximum(abs, [f(0.0, 0.37), f(1.0, 0.37), f(0.42, 0.0), f(0.42, 1.0)]) < 1e-15,
        @sprintf("max |f| on the four edges = %.3e",
            maximum(abs, [f(0.0, 0.37), f(1.0, 0.37), f(0.42, 0.0), f(0.42, 1.0)])))
end

# =============================================================================================
header("2. the homogeneous-Dirichlet space on [0,1]²")

const sq8 = EulerSquare(8, 2)
const sq16 = EulerSquare(16, 2)

check("|Ω| = 1 and the quadrature weights sum to it",
    isapprox(domainvolume(sq8.space), 1.0; atol = 1e-14) &&
        isapprox(sum(quadrature_weights(sq8.space)), 1.0; atol = 1e-13),
    @sprintf("|Ω| = %.16f   Σw = %.16f", domainvolume(sq8.space),
        sum(quadrature_weights(sq8.space))))

# ∫ sin(πx₁) sin(πx₂) = (2/π)². A function that vanishes on ∂Ω is representable here, so the
# projection reproduces its integral to the order of the space.
let exact = (2 / π)^2
    for sq in (sq8, sq16)
        v = integrate(sq, project(sq.space, x -> sin(π * x[1]) * sin(π * x[2])))
        check(@sprintf("∫sin(πx₁)sin(πx₂) = (2/π)² at n = %d", ncells(sq.space)[1]),
            abs(v - exact) / exact < 1e-3,
            @sprintf("%.12f vs %.12f   rel %.3e", v, exact, abs(v - exact) / exact))
    end
end

# The projection is idempotent on the space, which is what makes `euler_state` a projection
# rather than an interpolation with a name.
let û = project(sq16.space, x -> sin(π * x[1]) * sin(2π * x[2]))
    e = norm(project(sq16.space, x -> evaluate(sq16.space, û, (x[1], x[2]))) .- û) / norm(û)
    check("the L² projection is idempotent", e < 1e-11, @sprintf("rel %.3e", e))
end

# The stiffness matrix is NONSINGULAR, unlike the periodic case: there is no constant in the
# kernel because there is no constant in the space. That is what removes `spline.jl`'s bordered
# solve, and it is the same fact that makes §5.4's references exact.
let K = Matrix(stiffness_matrix(sq16.space)), M = Matrix(sq16.M)
    λ = eigvals(Symmetric(K), Symmetric(M))
    check("K is positive definite (no constant in the kernel)", minimum(λ) > 1.0,
        @sprintf("λ_min = %.8f   λ_max = %.4e", minimum(λ), maximum(λ)))
    # The constant function is not in the space: its best approximation is far from constant.
    c = project(sq16.space, x -> 1.0)
    v = [evaluate(sq16.space, c, (x, 0.5)) for x in (0.5, 0.1, 0.02, 0.0)]
    check("the constant function is NOT in the space",
        abs(v[end]) < 1e-14 && v[3] < 0.9,
        @sprintf("proj(1) at x₁ = 0.5, 0.1, 0.02, 0: %.4f %.4f %.4f %.4f", v...))
end

# λ₁,₁ = 2π², and it converges. The run can only ever reach the space's own eigenvalue, so both
# numbers get reported wherever a relaxed state is measured against one.
let e8 = dirichlet_eigenvalue(sq8), e16 = dirichlet_eigenvalue(sq16)
    check_refined("the discrete first eigenvalue → 2π² = 19.7392088022",
        e8 - DIRICHLET_EIGENVALUE, e16 - DIRICHLET_EIGENVALUE; atol = 1e-12, minrate = 8.0)
    check(@sprintf("λ_h = %.10f at n = 16", e16),
        abs(e16 - DIRICHLET_EIGENVALUE) / DIRICHLET_EIGENVALUE < 1e-5,
        @sprintf("rel %.3e   (2π² = %.10f)",
            abs(e16 - DIRICHLET_EIGENVALUE) /
            DIRICHLET_EIGENVALUE, DIRICHLET_EIGENVALUE))
end

# Λ inverts -Δ with the Dirichlet condition: on the first eigenmode it is division by λ_h, and
# φ vanishes on ∂Ω identically because every basis function does.
let sq = sq16, λh = dirichlet_eigenvalue(sq)
    e = project(sq.space, x -> sin(π * x[1]) * sin(π * x[2]))
    r = norm(sq.Λ * e .- e ./ λh) / norm(e ./ λh)
    check("Λ scales the first eigenmode by 1/λ_h", r < 1e-6, @sprintf("rel %.3e", r))
    φ̂ = sq.Λ * project(sq.space, x -> sin(π * x[1]) * sin(2π * x[2]) + 0.3sin(3π * x[1]))
    b = maximum(abs, [evaluate(sq.space, φ̂, p)
                      for p in ((0.0, 0.4), (1.0, 0.4), (0.4, 0.0), (0.4, 1.0))])
    check("φ = Λω vanishes on ∂Ω exactly", b < 1e-15, @sprintf("max |φ| on ∂Ω = %.3e", b))
    A = sq.MΛ
    check("MΛ is symmetric, as the energy's QuadraticHamiltonian needs",
        norm(A - A') / norm(A) < 1e-14, @sprintf("rel %.3e", norm(A - A') / norm(A)))
end

# The mesh and B3's floor are both fixed by one requirement: `y log y` needs ω > 0 at every
# quadrature node, and so does its mobility M = y. Two separate things stand in the way of that
# — resolution and amplitude — and this section separates them.
header("2b. the mesh, B3's initial condition, and the floor")

const b3 = SECTION5_RUNS["b3"]
const b3printed = EulerSpec("b3-printed", b3.section, b3.entropy, b3.state, b3.gaussian,
    nothing, b3.Δt, b3.T)

# The undershoot is asserted to a tolerance and not merely printed. It is quoted in `B3_FLOOR`'s
# and `EulerSquare`'s docstrings, and a sign-only check let the 20-cell row drift 9 % away from
# the value in the prose before anyone noticed. `rtol = 1e-3` is loose enough for the platform
# spread across the CI matrix and tight enough to catch that.
const B3_UNDERSHOOT = Dict(12 => -1.374746e-2, 16 => -5.830994e-4, 20 => -4.249944e-7)

# (i) The mesh threshold, and it is about the NARROW direction rather than about the boundary:
# w₁ = 0.1, and below 26 cells the L² projection oscillates below zero next to the peak. The
# three coarser rows are the control that 26 is a threshold and not a preference.
for n in (12, 16, 20, 26)
    sqn = EulerSquare(n, 2)
    (lo, hi) = state_extrema(sqn, euler_state(sqn, b3printed))
    check(
        @sprintf("n = %2d: the printed ω₀ is %s, as the 26-cell threshold says", n,
            n >= 26 ? "positive" : "NOT positive"),
        (lo > 0) == (n >= 26),
        @sprintf("min = %+.4e   max = %.4f", lo, hi))
    # At 26 cells min ω₀ is a cancellation residual, so only its sign is reproducible; the three
    # coarser rows are real undershoots and carry the numbers the docstrings quote.
    haskey(B3_UNDERSHOOT, n) && check(
        @sprintf("n = %2d: and it is the undershoot the docstrings quote, %+.3e", n,
            B3_UNDERSHOOT[n]),
        isapprox(lo, B3_UNDERSHOOT[n]; rtol = 1e-3),
        @sprintf("min = %+.6e   expected %+.6e   rel = %.2e", lo, B3_UNDERSHOOT[n],
            abs(lo - B3_UNDERSHOOT[n]) / abs(B3_UNDERSHOOT[n])))
end

# (ii) Positive is not admissible. Even at 26 cells the printed condition sits at the Gaussian's
# own far-corner value, 1.4e-11 of its peak, and `s = y log y` has no admissible set around a
# state a Galerkin scheme cannot hold away from zero.
for n in (26, 32)
    sqn = EulerSquare(n, 2)
    (lo, _) = state_extrema(sqn, euler_state(sqn, b3printed))
    check(@sprintf("n = %2d: the printed ω₀ has no MARGIN", n), 0 < lo < 1e-9,
        @sprintf("min ω₀ = %+.4e against a peak of 10 — the Gaussian's own value there",
            lo))
end

# (iii) The floor is what restores one, and this is the number that justifies its size.
for n in (12, 26)
    sqn = EulerSquare(n, 2)
    (lo, hi) = state_extrema(sqn, euler_state(sqn, b3))
    check(@sprintf("n = %2d: the floored ω₀ has a margin", n), lo > 1e-3,
        @sprintf("min ω₀ = %+.4e — %.1f %% of the floor %.3g; max = %.4f",
            lo, 100lo / B3_FLOOR, B3_FLOOR, hi))
end

# The free space's Λ still solves the Dirichlet problem: φ vanishes on ∂Ω while ω need not.
let sqf = EulerSquare(16, 2; state = :free)
    ω̂ = project(sqf.space, x -> 1.0)
    φ̂ = sqf.Λ * ω̂
    b = maximum(abs, [evaluate(sqf.space, φ̂, p)
                      for p in ((0.0, 0.4), (1.0, 0.4), (0.4, 0.0), (0.4, 1.0))])
    check("free space: φ = Λω still vanishes on ∂Ω while ω = 1 does not", b < 1e-14,
        @sprintf("max |φ| on ∂Ω = %.3e   with ω ≡ 1 in the space (ω(0,0.4) = %.6f)",
            b, evaluate(sqf.space, ω̂, (0.0, 0.4))))
    # -Δφ = 1 on the unit square has φ(½,½) = 0.0736713…, the classic torsion constant.
    check("free space: -Δφ = 1 gives the known centre value 0.07367135",
        abs(evaluate(sqf.space, φ̂, (0.5, 0.5)) - 0.0736713532) < 1e-6,
        @sprintf("φ(½,½) = %.10f   exact = 0.0736713532",
            evaluate(sqf.space, φ̂, (0.5, 0.5))))
    A = sqf.MΛ
    check("free space: MΛ is symmetric and positive semi-definite",
        norm(A - A') / norm(A) < 1e-14 && minimum(eigvals(Symmetric(A))) > -1e-14,
        @sprintf("asymmetry %.3e   λ_min = %+.3e", norm(A - A') / norm(A),
            minimum(eigvals(Symmetric(A)))))
    check("free space: the constant IS in it, so the mass is a discrete Casimir",
        abs(evaluate(sqf.space, project(sqf.space, x -> 1.0), (0.0, 0.5)) - 1) < 1e-12,
        @sprintf("proj(1) at the boundary = %.14f",
            evaluate(sqf.space, project(sqf.space, x -> 1.0), (0.0, 0.5))))
end

# =============================================================================================
header("3. the collision-like bracket on the Dirichlet space")

const sq = EulerSquare(12, 2)
const N = nbasis(sq.space)

for name in SECTION5_ORDER
    spec = SECTION5_RUNS[name]
    # Each run's own space, because that is part of the run — see `SECTION5_RUNS`.
    sqn = EulerSquare(12, 2; state = spec.state)
    f = euler_flow(sqn, spec)
    # A state that is positive at every node, so that B3's mobility M = ω is admissible, and
    # not an eigenmode, so that the vector field is not accidentally zero. The `+0.4` is what
    # keeps it away from zero at the boundary, and it survives only in the free space.
    ω̂ = project(sqn.space,
        x -> sin(π * x[1]) * sin(π * x[2]) * (1.5 + 0.4sin(2π * x[1]) - 0.3cos(π * x[2])) +
             (spec.state === :free ? 0.4 : 0.0))
    (lo, _) = state_extrema(sqn, ω̂)
    check(@sprintf("%-3s the test state is admissible (ω > 0) on :%s", name, spec.state),
        lo > 0, @sprintf("min ω = %.4e", lo))
    check(@sprintf("%-3s 𝔾 is symmetric", name), issymmetric(f.metric, ω̂), "")
    check(@sprintf("%-3s 𝔾 is positive semi-definite", name),
        ispositive_semidefinite(f.metric, ω̂), "")
    r = degeneracy_residual(f, ω̂)
    check(@sprintf("%-3s (F, H) = 0 — 𝔾 ∂H/∂û = 0 on this space", name), r < 1e-11,
        @sprintf("normalised residual %.3e", r))
    # The two evaluations of the same operator: `metric_apply` goes through the 𝔽_s moments,
    # `metric_matrix` through the R/S/T factorisation. Agreement is a cross-check.
    c = Vector(sqn.M * ω̂)
    a, m = metric_apply(f.metric, ω̂, c), metric_matrix(f.metric, ω̂) * c
    check(@sprintf("%-3s metric_apply agrees with metric_matrix", name),
        norm(a .- m) / norm(m) < 1e-10, @sprintf("rel %.3e", norm(a .- m) / norm(m)))
    v = vectorfield(f, ω̂)
    cosH = abs(dot(gradient(f, ω̂), v)) / (norm(gradient(f, ω̂)) * norm(v))
    check(@sprintf("%-3s the vector field is orthogonal to ∂H/∂û", name), cosH < 1e-11,
        @sprintf("cos = %.3e", cosH))
    p = entropy_production(f, ω̂)
    check(@sprintf("%-3s (S,S) > 0 and equals -dS/dt", name),
        p > 0 && isapprox(p, -dot(entropy_gradient(f, ω̂), v); rtol = 1e-12),
        @sprintf("(S,S) = %.10e   -Ṡ = %.10e", p, -dot(entropy_gradient(f, ω̂), v)))
end

# The two controls. Neither is a hypothetical: the first is the trap `MetriplecticFlow`'s own
# docstring warns about, and the second is `eq:M-condition`'s requirement M > 0 being violated.
let ω̂ = project(sq.space, x -> sin(π * x[1]) * sin(π * x[2]) * (1.5 + 0.4sin(2π * x[1])))
    # A bracket generated by an UNRELATED field is a perfectly good metric bracket that is
    # degenerate on the wrong energy. The two-argument degeneracy check cannot see it, because
    # it takes the gradient from the bracket itself.
    #
    # A prescribed field rather than another elliptic solve, and that is not laziness: with
    # Λ' = (K+M)⁻¹M the wrong stream function is 0.95 φ, and degeneracy only sees direction, so
    # the residual comes out at 3e-6 and the control is nearly vacuous.
    ĥ = project(sq.space, x -> sin(2π * x[1]) * sin(3π * x[2]))
    G = CollisionBracket(sq.space, ĥ)
    fw = MetriplecticFlow(sq.space, G, QuadraticHamiltonian(sq.MΛ),
        QuadraticHamiltonian(Matrix(sq.M)))
    rb, rf = degeneracy_residual(G, ω̂), degeneracy_residual(fw, ω̂)
    # The content of this control is the RATIO, not an absolute threshold. `degeneracy_residual`
    # normalises by max|𝔾| max|g|, which bounds a single product rather than the N-term sum with
    # its cancellation, so a complete misalignment still normalises to 1e-3 — thirteen orders
    # above the clean 1e-16, and that gap is what the check asserts.
    check("CONTROL: a mismatched generator breaks the FLOW's degeneracy", rf > 1e8 * rb,
        @sprintf("bracket residual %.3e (clean), flow residual %.3e, ratio %.2e",
            rb, rf, rf / rb))
    # A sign-indefinite mobility violates eq:M-condition's M > 0 and must lose positivity.
    # x₁ - 0.3, not x₁ - 0.5: the latter integrates to zero over the square, so m₀ = ∫M dμ = 0
    # and the bracket's own recentring divides by it. That breaks the degeneracy too, which
    # would make the next check fail for a reason that has nothing to do with the sign.
    Gi = CollisionBracket(sq.space, sq.Λ; mobility = (x, u) -> x[1] - 0.3,
        mobility_derivative = 0)
    check("CONTROL: an indefinite mobility is NOT positive semi-definite",
        !ispositive_semidefinite(Gi, ω̂), "")
    # ...and is still symmetric and still degenerate, which is why positivity is a separate
    # check rather than a corollary.
    check("CONTROL: it is nevertheless symmetric and degenerate",
        issymmetric(Gi, ω̂) && degeneracy_residual(Gi, ω̂) < 1e-11,
        @sprintf("degeneracy %.3e", degeneracy_residual(Gi, ω̂)))
end

# =============================================================================================
header("4. the two entropies")

"An independent quadrature: the tensor midpoint rule on an n-by-n grid, which shares no code
path with the space's own Gauss-Legendre tables."
function midpoint_integral(g, s, û, n)
    h = SQUARE_LENGTH / n
    total = 0.0
    for j in 1:n, i in 1:n

        x = ((i - 0.5) * h, (j - 0.5) * h)
        total += g(evaluate(s, û, x))
    end
    return total * h^2
end

# The free space, and bounded away from zero, for both halves of this section. `y log y` on a
# state that vanishes on ∂Ω has an unbounded derivative there: its midpoint quadrature then
# converges at first order rather than second, and the comparison would be measuring the
# reference rule instead of the entropy.
const sqe = EulerSquare(16, 2; state = :free)

let s = sqe.space,
    ω̂ = project(s, x -> sin(π * x[1]) * sin(π * x[2]) * (1.5 + 0.4sin(2π * x[1])) + 0.4)

    (lo, _) = state_extrema(sqe, ω̂)
    check("the entropy test state is admissible (ω > 0)", lo > 0,
        @sprintf("min ω = %.4e", lo))

    Squad = hamiltonian(QuadraticHamiltonian(Matrix(sqe.M)), s, ω̂)
    check_refined("½∫ω² matches an independent midpoint quadrature",
        Squad - midpoint_integral(y -> y^2 / 2, s, ω̂, 120),
        Squad - midpoint_integral(y -> y^2 / 2, s, ω̂, 240); atol = 1e-9, minrate = 3.5)

    Sg = hamiltonian(GibbsEntropy(), s, ω̂)
    check_refined("∫ω log ω matches an independent midpoint quadrature",
        Sg - midpoint_integral(y -> y * log(y), s, ω̂, 120),
        Sg - midpoint_integral(y -> y * log(y), s, ω̂, 240); atol = 1e-9, minrate = 3.5)

    # The analytic gradient and Hessian against central differences. Not an alternative
    # implementation for the flow to use -- a differenced Jacobian has a floor four orders above
    # the conservation the runs measure -- but the right way to test the analytic one.
    for (label, H) in (("½∫ω²", QuadraticHamiltonian(Matrix(sqe.M))),
        ("∫ω log ω", GibbsEntropy()))
        g = gradient(H, s, ω̂)
        h = hessian(H, s, ω̂)
        worst_g, worst_h = 0.0, 0.0
        for _ in 1:4
            v = randn(length(ω̂))
            v ./= norm(v)
            ε = 1e-6 * norm(ω̂)
            dS = (hamiltonian(H, s, ω̂ .+ ε .* v) - hamiltonian(H, s, ω̂ .- ε .* v)) / 2ε
            worst_g = max(worst_g, abs(dS - dot(g, v)) / abs(dot(g, v)))
            dg = (gradient(H, s, ω̂ .+ ε .* v) .- gradient(H, s, ω̂ .- ε .* v)) ./ 2ε
            worst_h = max(worst_h, norm(dg .- h * v) / norm(h * v))
        end
        check(@sprintf("%-10s analytic gradient = central difference", label),
            worst_g < 1e-6, @sprintf("worst rel %.3e", worst_g))
        check(@sprintf("%-10s analytic Hessian = central difference", label),
            worst_h < 1e-5, @sprintf("worst rel %.3e", worst_h))
    end

    # eq:M-condition, M ∂²_y s = 1, ties the mobility to the entropy through its SECOND
    # derivative, and the Gibbs entropy is where that has content: ∂²_y s = 1/y, so the
    # Hessian's weight is 1/ω and the mobility that inverts it is M = ω. On a constant state
    # ω ≡ c the condition is readable straight off the assembled Hessian — c ∇²S must be the
    # plain mass matrix — and a Hessian weighted by ω rather than by 1/ω misses it by c².
    # The quadratic entropy needs no such row: its Hessian IS the mass matrix by construction,
    # with M = 1, so comparing the two would assert nothing.
    let Mm = Matrix(mass_matrix(s)), c = 2.5, ĉ = project(s, x -> c),
        e = maximum(abs, c .* hessian(GibbsEntropy(), s, ĉ) .- Mm) / maximum(abs, Mm)

        check("eq:M-condition: M ∂²_y s = 1, with M = ω for s = ω log ω", e < 1e-10,
            @sprintf("at ω ≡ %.1f: max |c ∇²S − 𝕄| / max|𝕄| = %.3e   (M = 1 for ω²/2 is 𝕄 itself)",
                c, e))
    end
end

# =============================================================================================
header("5. the closed-form references")

let sq = sq16, λh = dirichlet_eigenvalue(sq)
    K = Matrix(stiffness_matrix(sq.space))
    F = eigen(Symmetric(K), Symmetric(Matrix(sq.M)))
    ê = F.vectors[:, 1] ./ norm(F.vectors[:, 1])

    # S_η = λ H is an identity at ω = λφ, and it is the identity `euler_entropy_floor` states.
    Hη = dot(ê, sq.MΛ, ê) / 2
    Sη = l2inner(sq, ê, ê) / 2
    check("S = λ_h H at the eigenmode, which is S_η",
        isapprox(Sη, euler_entropy_floor(Hη; λ = λh); rtol = 1e-10),
        @sprintf("S = %.12e   λ_h H = %.12e   rel %.3e", Sη,
            euler_entropy_floor(Hη; λ = λh),
            abs(Sη - euler_entropy_floor(Hη; λ = λh)) / Sη))

    # `eigenmode_fit` returns λ_h and a vanishing residual on the eigenmode exactly.
    (λ, _, rel) = eigenmode_fit(sq, ê)
    check("eigenmode_fit recovers λ_h with no residual",
        isapprox(λ, λh; rtol = 1e-10) && rel < 1e-10,
        @sprintf("λ = %.10f   λ_h = %.10f   rel residual %.3e", λ, λh, rel))

    # The Poincaré floor S ≥ λ_h H holds for EVERY state of the space, not only along the flow,
    # and the eigenmode attains it. That is why the drivers assert it at t = 0 as well.
    worst = Inf
    for _ in 1:200
        v = randn(nbasis(sq.space))
        worst = min(worst, l2inner(sq, v, v) / dot(v, sq.MΛ, v))
    end
    check("S/H ≥ λ_h for every state, attained at the eigenmode",
        worst > λh * (1 - 1e-9) && isapprox(Sη / Hη, λh; rtol = 1e-10),
        @sprintf("min S/H over 200 random states = %.8f   λ_h = %.8f", worst, λh))
    check("...and λ_h ≥ 2π², so S ≥ 2π² H is implied", λh >= DIRICHLET_EIGENVALUE,
        @sprintf("λ_h - 2π² = %+.3e", λh - DIRICHLET_EIGENVALUE))
end

# `gibbs_lambda` is an identity rather than a fit, and its two inputs are what have to be right:
# ∫ωφ = 2H and ∫ω = M(ω). Given those, λ = (M+S)/2H follows from log ω = λφ - 1 in one line.
let sq = sq16,
    ω̂ = project(sq.space, x -> sin(π * x[1]) * sin(π * x[2]) * (1.5 + 0.4sin(2π * x[1])))

    φ̂ = sq.Λ * ω̂
    twoH = dot(ω̂, sq.MΛ, ω̂)
    check("∫ωφ = 2H, which is gibbs_lambda's denominator",
        isapprox(l2inner(sq, ω̂, φ̂), twoH; rtol = 1e-12),
        @sprintf("(ω,φ) = %.12e   2H = %.12e", l2inner(sq, ω̂, φ̂), twoH))
    check_refined("∫ω = integrate(sq, ω̂) against an independent quadrature",
        integrate(sq, ω̂) - midpoint_integral(identity, sq.space, ω̂, 120),
        integrate(sq, ω̂) - midpoint_integral(identity, sq.space, ω̂, 240);
        atol = 1e-12, minrate = 3.5)

    # The identity itself, on a state MANUFACTURED to satisfy log ω = λφ + μ - 1 exactly at the
    # quadrature nodes -- which is the algebra `gibbs_lambda` rests on, isolated from whether any
    # run reaches such a state.
    w = quadrature_weights(sq.space)
    φ = field(sq.space, φ̂, (0, 0))
    for (λ, μ) in ((3.0, 0.0), (7.5, 0.0), (3.0, 0.4))
        u = exp.(λ .* φ .+ μ .- 1)
        S = dot(w, u .* log.(u))
        Mm = dot(w, u)
        Hm = dot(w, u .* φ) / 2
        check(@sprintf("λ = (M + (1-μ)S… ) closes at λ = %.1f, μ = %.1f", λ, μ),
            isapprox(S, 2λ * Hm + (μ - 1) * Mm; rtol = 1e-12),
            @sprintf("S = %.12e   2λH + (μ-1)M = %.12e", S, 2λ * Hm + (μ - 1) * Mm))
        # And the manuscript's own formula is the μ = 0 case -- which is why the space matters.
        check(
            @sprintf("gibbs_lambda is exact at μ = 0 and wrong by %.0f%% at μ = %.1f",
                100abs(μ), μ),
            isapprox(gibbs_lambda(Mm, S, Hm), λ; rtol = 1e-10) == (μ == 0),
            @sprintf("(M+S)/2H = %.10f   λ = %.10f", gibbs_lambda(Mm, S, Hm), λ))
    end
end

# `gibbs_fit` must recover both multipliers from data that really is of that form. Manufactured
# in the FREE space, because e^{λφ+μ-1} does not vanish on ∂Ω and is not representable in the
# Dirichlet one -- which is B3's whole story.
let sqf = EulerSquare(16, 2; state = :free)
    ψ̂ = sqf.Λ * project(sqf.space, x -> 1.0)                # a smooth positive Dirichlet field
    for (λ, μ) in ((3.0, 0.0), (7.5, 0.4))
        # The fit reads φ from Λω̂ internally, so `ω̂` has to be a state whose OWN stream
        # function is the one in the exponent: a few fixed-point sweeps get there.
        ω̂ = project(sqf.space, x -> exp(λ * evaluate(sqf.space, ψ̂, (x[1], x[2])) + μ - 1))
        for _ in 1:60
            φ̂ = sqf.Λ * ω̂
            ω̂ = project(sqf.space,
                x -> exp(λ * evaluate(sqf.space, φ̂, (x[1], x[2])) + μ - 1))
        end
        (λf, μf, rf) = gibbs_fit(sqf, ω̂)
        check(@sprintf("gibbs_fit recovers (λ, μ) = (%.1f, %.1f)", λ, μ),
            abs(λf - λ) / λ < 1e-3 && abs(μf - μ) < 1e-3,
            @sprintf("fit λ = %.8f  μ = %+.8f   relative fit residual %.3e", λf, μf, rf))
        check(
            @sprintf("...and gibbs_residual sees it as the family member it is (μ = %.1f)",
                μ),
            gibbs_residual(sqf, ω̂, λf, μf) < 1e-4,
            @sprintf("‖ω - e^{λφ+μ-1}‖/‖ω‖ = %.3e   (against μ = 0: %.3e)",
                gibbs_residual(sqf, ω̂, λf, μf), gibbs_residual(sqf, ω̂, λf)))
    end
end

# `interior_weights` is only used to say that B3's residual is not a boundary artefact, so what
# has to be right about it is the area it drops.
let sq = sq16, h = 1 / 16
    w0 = interior_weights(sq)
    for m in (h, 2h, 4h)
        wm = interior_weights(sq; margin = m)
        lost = (sum(w0) - sum(wm)) / sum(w0)
        expected = 1 - (1 - 2m)^2
        check(@sprintf("interior_weights(margin = %.4f) drops the right area", m),
            abs(lost - expected) < 0.05,
            @sprintf("dropped %.4f of |Ω|, boundary strip is %.4f", lost, expected))
    end
    check("e^{λφ-1} is e^{-1} on ∂Ω, hence NOT in the Dirichlet space",
        isapprox(exp(-1), 0.36787944117144233; rtol = 1e-15),
        @sprintf("e^{-1} = %.16f   while every ω_h ∈ V_D is 0 there", exp(-1)))
end

# =============================================================================================
header("6. Crank-Nicolson on a short run")

# ImplicitMidpoint IS Crank-Nicolson for this field. What is asserted here is the pair of
# semi-discrete conservation laws surviving the time discretisation, at a resolution and a step
# count that take seconds; the relaxation itself is the drivers' business.
for name in SECTION5_ORDER
    spec = SECTION5_RUNS[name]
    sqc = EulerSquare(12, 2; state = spec.state)
    d = Diagnostics(sqc, spec)
    f = euler_flow(sqc, spec)
    ω̂ = euler_state(sqc, spec)
    # B3's step size is bounded by admissibility rather than by accuracy — see section 2b and
    # the sweep in `run_b3.jl` — so it gets its own, well inside the boundary measured there.
    Δt = name == "b3" ? 0.01 : 0.5
    integ = Integrator(f, ImplicitMidpoint(), Δt; û₀ = ω̂)
    tr = Trace(ω̂)
    record!(tr, d, 0.0, ω̂)
    for k in 1:6
        integrate_step!(ω̂, integ)
        record!(tr, d, k * Δt, ω̂)
    end
    e = maximum(energy_error(tr))
    check(@sprintf("%-3s H conserved to machine precision over 6 steps", name), e < 1e-13,
        @sprintf("max |ΔH|/|H₀| = %.3e   H₀ = %.12e   Δt = %.4g", e, tr.H[1], Δt))
    (ok, worst) = entropy_monotone(tr)
    check(@sprintf("%-3s S monotone over 6 steps", name), ok,
        @sprintf("worst increment %+.3e   S: %.10e → %.10e", worst, tr.S[1], tr.S[end]))
    # Only B3 has an admissible set to stay inside: `y log y` and its mobility M = y are
    # undefined at ω ≤ 0. For the quadratic entropy the sign is unconstrained — B2's own ω₀ is
    # deliberately sign-indefinite — so the extrema are REPORTED there rather than asserted. A
    # row reading `[PASS] b2 ω stays admissible   min ω = -9.6e-01` is a verdict no run earned.
    (lo, _) = state_extrema(sqc, ω̂)
    if name == "b3"
        check(@sprintf("%-3s ω stays admissible", name), lo > 0,
            @sprintf("min ω = %+.4e", lo))
    else
        check(@sprintf("%-3s ω sign is unconstrained  [REPORTED]", name), true,
            @sprintf("min ω = %+.4e", lo))
    end
end

# The step-halving difference: Δt is a COST choice here and not an accuracy one, and this is the
# measurement that says so. Six steps at Δt against twelve at Δt/2, compared at the same time.
let name = "b1", spec = SECTION5_RUNS[name], sqc = EulerSquare(12, 2)
    function endpoint(Δt, n)
        f = euler_flow(sqc, spec)
        ω̂ = euler_state(sqc, spec)
        integ = Integrator(f, ImplicitMidpoint(), Δt; û₀ = ω̂)
        for _ in 1:n
            integrate_step!(ω̂, integ)
        end
        return ω̂
    end
    a, b = endpoint(1.0, 6), endpoint(0.5, 12)
    e = l2norm(sqc, a .- b) / l2norm(sqc, b)
    check("Δt = 1 and Δt = 0.5 land on the same state at t = 6", e < 1e-3,
        @sprintf("relative L² difference %.3e", e))
end

summary("verify_euler.jl")
