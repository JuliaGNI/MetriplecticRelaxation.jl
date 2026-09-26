using MetriplecticRelaxation
using MetriplecticRelaxation: integrate
using MetriplecticRelaxation: EulerSquare, dirichlet_eigenvalue, DIRICHLET_EIGENVALUE,
                              euler_axis, euler_grid
using GeometricBrackets: nbasis, project, evaluate, domainvolume, stiffness_matrix,
                         quadrature_weights
using SimpleSplines: Dirichlet
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
@testset "$(rpad("Euler Square Tests", 80))" begin
    sq = EulerSquare(10, 2)
    N = nbasis(sq.space)

    @testset "$(rpad("The space IS the unit square, and integrates it", 76))" begin
        @test domainvolume(sq.space)≈1.0 atol=1e-14
        @test sum(quadrature_weights(sq.space))≈1.0 atol=1e-13
        # ∫ sin(πx₁)sin(πx₂) = (2/π)², for a function that vanishes on ∂Ω and is therefore
        # representable here.
        v = integrate(sq, project(sq.space, x -> sin(π * x[1]) * sin(π * x[2])))
        @test v≈(2 / π)^2 rtol=1e-3
        û = project(sq.space, x -> sin(π * x[1]) * sin(2π * x[2]))
        @test project(sq.space, x -> evaluate(sq.space, û, (x[1], x[2])))≈û rtol=1e-11
    end

    @testset "$(rpad("Constants are NOT in the space, so K is nonsingular", 76))" begin
        # Both halves of one fact, and it is the fact both §5.4 references rest on.
        K = Matrix(stiffness_matrix(sq.space))
        @test minimum(eigvals(Symmetric(K), Symmetric(Matrix(sq.M)))) > 1.0
        c = project(sq.space, x -> 1.0)
        @test abs(evaluate(sq.space, c, (0.0, 0.5))) < 1e-14
        @test evaluate(sq.space, c, (0.02, 0.5)) < 0.9
    end

    @testset "$(rpad("Lambda INVERTS the Dirichlet Laplacian", 76))" begin
        λh = dirichlet_eigenvalue(sq)
        @test λh≈DIRICHLET_EIGENVALUE rtol=1e-4
        @test λh >= DIRICHLET_EIGENVALUE
        e = project(sq.space, x -> sin(π * x[1]) * sin(π * x[2]))
        @test norm(sq.Λ * e .- e ./ λh) / norm(e) < 1e-5
        # φ = Λω vanishes on ∂Ω identically, because every basis function does.
        φ̂ = sq.Λ * project(sq.space, x -> sin(π * x[1]) * sin(2π * x[2]))
        for p in ((0.0, 0.4), (1.0, 0.4), (0.4, 0.0), (0.4, 1.0))
            @test abs(evaluate(sq.space, φ̂, p)) < 1e-15
        end
        @test norm(sq.MΛ - sq.MΛ') / norm(sq.MΛ) < 1e-14
    end

    @testset "$(rpad("euler_grid SAMPLES x_i along the FIRST index, endpoints included", 76))" begin
        xs = euler_axis(21)
        @test length(xs) == 21
        @test xs[1] == 0.0
        @test xs[end] == 1.0

        # Every function in V_D vanishes on ∂Ω, for RANDOM degrees of freedom and not merely
        # for the projection of something that vanished there already.
        Z = euler_grid(sq, randn(N), 21)
        @test maximum(abs, vcat(Z[1, :], Z[end, :], Z[:, 1], Z[:, end])) < 1e-13

        # An asymmetric polynomial the degree-3 space contains exactly: reproduced to round-off,
        # and wrong by order one read the other way round. Degree 3 rather than the runs' own
        # degree 2 because the boundary-vanishing bi-degree-(2,2) polynomials are
        # span{x(1-x)y(1-y)}, which is symmetric and would pass a transposed implementation.
        f(x, y) = x * (1 - x) * y^2 * (1 - y)
        sq3 = EulerSquare(5, 3)
        W = euler_grid(sq3, project(sq3.space, x -> f(x[1], x[2])), 21)
        @test maximum(abs(W[i, j] - f(xs[i], xs[j])) for i in 1:21, j in 1:21) < 1e-14
        @test maximum(abs(W[i, j] - f(xs[j], xs[i])) for i in 1:21, j in 1:21) > 1e-3
    end
end
