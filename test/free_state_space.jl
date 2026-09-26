using MetriplecticRelaxation
using MetriplecticRelaxation: EulerSquare, dirichlet_eigenvalue, DIRICHLET_EIGENVALUE
using GeometricBrackets: project, evaluate
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
@testset "$(rpad("Free State Space Tests", 80))" begin
    # B3's space: ω unconstrained on ∂Ω, φ still the homogeneous-Dirichlet solve, reached
    # through the recombination matrix. That split is the finding §5.4 turned up.
    sqf = EulerSquare(12, 2; state = :free)

    @testset "$(rpad("The CONSTANT is in the free space and phi still vanishes", 76))" begin
        c = project(sqf.space, x -> 1.0)
        @test evaluate(sqf.space, c, (0.0, 0.5))≈1.0 atol=1e-12
        φ̂ = sqf.Λ * c
        for p in ((0.0, 0.4), (1.0, 0.4), (0.4, 0.0), (0.4, 1.0))
            @test abs(evaluate(sqf.space, φ̂, p)) < 1e-14
        end
        # -Δφ = 1 on the unit square: φ(½,½) is the classical torsion constant 0.07367135…,
        # reached to 2e-6 at these 12 cells and to 7e-7 at the 16 of `verify_euler.jl`.
        @test evaluate(sqf.space, φ̂, (0.5, 0.5))≈0.0736713532 atol=1e-5
    end

    @testset "$(rpad("M Lambda IS symmetric positive semi-definite there too", 76))" begin
        @test norm(sqf.MΛ - sqf.MΛ') / norm(sqf.MΛ) < 1e-14
        @test minimum(eigvals(Symmetric(sqf.MΛ))) > -1e-14
        @test dirichlet_eigenvalue(sqf)≈DIRICHLET_EIGENVALUE rtol=1e-3
        @test dirichlet_eigenvalue(sqf) >= DIRICHLET_EIGENVALUE
    end

    @testset "$(rpad("A bad state space name IS rejected", 76))" begin
        @test_throws ArgumentError EulerSquare(4, 2; state = :neumann)
    end
end
