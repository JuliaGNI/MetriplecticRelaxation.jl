using MetriplecticRelaxation
using MetriplecticRelaxation: l2inner
using MetriplecticRelaxation: EulerSquare, euler_entropy_floor, eigenmode_fit,
                              dirichlet_eigenvalue, gibbs_lambda, interior_weights,
                              DIRICHLET_EIGENVALUE
using GeometricBrackets: nbasis, project, stiffness_matrix, field, quadrature_weights
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
@testset "$(rpad("Section 5 Reference Tests", 80))" begin
    sq = EulerSquare(12, 2)
    λh = dirichlet_eigenvalue(sq)
    F = eigen(Symmetric(Matrix(stiffness_matrix(sq.space))), Symmetric(Matrix(sq.M)))
    ê = F.vectors[:, 1] ./ norm(F.vectors[:, 1])

    @testset "$(rpad("S = lambda H IS an identity at the eigenmode", 76))" begin
        Hη = dot(ê, sq.MΛ, ê) / 2
        Sη = l2inner(sq, ê, ê) / 2
        @test Sη≈euler_entropy_floor(Hη; λ = λh) rtol=1e-10
        (λ, _, rel) = eigenmode_fit(sq, ê)
        @test λ≈λh rtol=1e-10
        @test rel < 1e-10
    end

    @testset "$(rpad("S/H >= lambda_h for EVERY state of the space", 76))" begin
        # The Poincaré inequality, which is why the drivers assert S ≥ S_η at t = 0 as a
        # precondition and not only at the end.
        worst = Inf
        for _ in 1:100
            v = randn(nbasis(sq.space))
            worst = min(worst, l2inner(sq, v, v) / dot(v, sq.MΛ, v))
        end
        @test worst > λh * (1 - 1e-9)
        @test λh >= DIRICHLET_EIGENVALUE
    end

    @testset "$(rpad("gibbs_lambda IS exact at mu = 0 and WRONG otherwise", 76))" begin
        # The manuscript's λ = (M+S)/2H₀ follows from log ω = λφ - 1 by substitution, so it
        # holds exactly when the mass multiplier μ vanishes -- which is what the
        # homogeneous-Dirichlet space forces -- and not otherwise.
        w = quadrature_weights(sq.space)
        ω̂ = project(sq.space, x -> sin(π * x[1]) * sin(π * x[2]) *
                                   (1.5 + 0.4sin(2π * x[1])))
        φ = field(sq.space, sq.Λ * ω̂, (0, 0))
        for (λ, μ) in ((3.0, 0.0), (7.5, 0.0), (3.0, 0.4))
            u = exp.(λ .* φ .+ μ .- 1)
            S, M, H = dot(w, u .* log.(u)), dot(w, u), dot(w, u .* φ) / 2
            @test S≈2λ * H + (μ - 1) * M rtol=1e-12
            @test isapprox(gibbs_lambda(M, S, H), λ; rtol = 1e-10) == (μ == 0)
        end
    end

    @testset "$(rpad("interior_weights DROPS exactly the boundary strip", 76))" begin
        # e^{λφ-1} is e^{-1} on ∂Ω while every ω_h of the space is zero there, so B3's
        # reference cannot hold in the last cell. This is the exclusion, and it is a strip.
        w0 = interior_weights(sq)
        @test sum(w0)≈1.0 atol=1e-13
        for m in (1 / 12, 2 / 12, 4 / 12)
            lost = (sum(w0) - sum(interior_weights(sq; margin = m))) / sum(w0)
            @test abs(lost - (1 - (1 - 2m)^2)) < 0.05
        end
    end
end
