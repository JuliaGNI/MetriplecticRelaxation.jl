using MetriplecticRelaxation
using MetriplecticRelaxation: SpectralTorus, torus_field, poisson_periodic,
                              projector_bracket_field, l2inner, l2norm, euler_minimiser,
                              euler_entropy_minimum, entropy
using GeometricBrackets: project, field
using LinearAlgebra
using Random
using SparseArrays
using Test

# The runs here excite fields with random degrees of freedom, which is deliberate: a
# conservation or degeneracy identity tested on a single smooth mode reports round-off and hides
# a structural defect entirely. The seed is fixed so that a failure is reproducible.
Random.seed!(0x5c1e9a3b)

# Every resolution below is deliberately far coarser than a Section 4 run, and every final time
# far shorter. These tests assert the STRUCTURAL claims -- conservation, monotonicity,
# degeneracy, the closed forms -- which hold on any mesh. The claims that need resolution to be
# true at all, the rates and the contour averages, are the run drivers' business and are
# measured in `scripts/`, not here.
@testset "$(rpad("Projector Rate Tests", 80))" begin
    # The linearised spectrum of the projector flow is exactly {1 - 1/λ} over the eigenvalues
    # of -Δ. That is what makes the manuscript's "≈ 1/2" and "≈ 1" exact rather than fitted.
    g = SpectralTorus(32)
    H₀ = 0.35
    ω★ = torus_field(g, euler_minimiser(H₀, 0.6, 0.0, 0.0))
    φ★ = poisson_periodic(g, ω★)

    function jac(v; ε = 1e-6)
        f(w) = projector_bracket_field(g, w, poisson_periodic(g, w))
        (f(ω★ .+ ε .* v) .- f(ω★ .- ε .* v)) ./ 2ε
    end

    @testset "$(rpad("The relaxed state IS a fixed point", 76))" begin
        @test φ★≈ω★ rtol=1e-13
        @test l2norm(g, projector_bracket_field(g, ω★, φ★)) / l2norm(g, ω★) < 1e-13
        @test l2inner(g, ω★, ω★) / 2≈euler_entropy_minimum(H₀) rtol=1e-13
    end

    @testset "$(rpad("An eigenmode of eigenvalue lambda DECAYS at 1 - 1/lambda", 76))" begin
        for (f, λ) in (((a, b) -> cos(a + b), 2.0), ((a, b) -> sin(a - b), 2.0),
            ((a, b) -> cos(2a), 4.0), ((a, b) -> cos(2a + b), 5.0),
            ((a, b) -> cos(3a), 9.0))
            v = torus_field(g, f)
            Jv = jac(v)
            @test -l2inner(g, Jv, v) / l2inner(g, v, v)≈1 - 1 / λ atol=1e-6
            # And J v really is parallel to v, so "the rate" is well defined.
            @test l2norm(g, Jv .+ (1 - 1 / λ) .* v) / l2norm(g, Jv) < 1e-5
        end
    end

    @testset "$(rpad("The whole lambda = 1 eigenspace IS neutral", 76))" begin
        for f in ((a, b) -> sin(a), (a, b) -> cos(b), (a, b) -> sin(b))
            v = torus_field(g, f)
            v .-= (l2inner(g, v, φ★) / l2inner(g, φ★, φ★)) .* φ★
            l2norm(g, v) < 1e-10 && continue
            @test abs(l2inner(g, jac(v), v) / l2inner(g, v, v)) < 1e-6
        end
        # Parallel to φ★ too: scaling moves along the family to another H₀, and every member
        # of it is a fixed point.
        @test abs(l2inner(g, jac(φ★), φ★) / l2inner(g, φ★, φ★)) < 1e-6
    end
end
