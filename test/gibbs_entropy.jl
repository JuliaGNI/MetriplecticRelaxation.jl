using MetriplecticRelaxation
using MetriplecticRelaxation: entropy
using MetriplecticRelaxation: EulerSquare, GibbsEntropy, state_extrema
using GeometricBrackets: project, gradient, hamiltonian, hessian, QuadraticHamiltonian,
                         field, quadrature_weights, basis_values
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
@testset "$(rpad("Gibbs Entropy Tests", 80))" begin
    # The free space, because that is where `y log y` is usable at all: a state bounded away
    # from zero needs a nonzero boundary trace.
    sq = EulerSquare(10, 2; state = :free)
    s = sq.space
    ω̂ = project(s, x -> sin(π * x[1]) * sin(π * x[2]) * (1.5 + 0.4sin(2π * x[1])) + 0.4)

    @testset "$(rpad("The analytic gradient and Hessian ARE the derivatives", 76))" begin
        @test state_extrema(sq, ω̂)[1] > 0
        for H in (QuadraticHamiltonian(Matrix(sq.M)), GibbsEntropy())
            g = gradient(H, s, ω̂)
            h = hessian(H, s, ω̂)
            for _ in 1:3
                v = randn(length(ω̂))
                v ./= norm(v)
                ε = 1e-6 * norm(ω̂)
                dS = (hamiltonian(H, s, ω̂ .+ ε .* v) -
                      hamiltonian(H, s, ω̂ .- ε .* v)) / 2ε
                @test dS≈dot(g, v) rtol=1e-6
                dg = (gradient(H, s, ω̂ .+ ε .* v) .- gradient(H, s, ω̂ .- ε .* v)) ./ 2ε
                @test norm(dg .- h * v) / norm(h * v) < 1e-5
            end
        end
    end

    @testset "$(rpad("The Hessian of y log y IS 1/M, which is eq:M-condition", 76))" begin
        # M ∂²_y s = 1 with s = y log y gives M = y, so the Hessian's weight is 1/ω exactly.
        u = field(s, ω̂, (0, 0))
        Φ = basis_values(s, (0, 0))
        @test hessian(GibbsEntropy(), s, ω̂) ≈
              Matrix(Φ * Diagonal(quadrature_weights(s) ./ u) * Φ') rtol=1e-14
        # And the premise that makes the weight 1/u finite in the first place: ω > 0 at every
        # quadrature node. `u .* (1 ./ u) ≈ 1` stood here, which is true of any nonzero float
        # and so asserted nothing.
        @test all(>(0), u)
    end

    @testset "$(rpad("y log y REFUSES a state it does not admit", 76))" begin
        # No guard: `log` raises, which is the correct outcome for a state the entropy is not
        # defined on. A guard would turn it into a silently wrong number.
        @test_throws DomainError hamiltonian(GibbsEntropy(), s, .-ω̂)
    end
end
