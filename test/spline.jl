using MetriplecticRelaxation
using MetriplecticRelaxation: SplineTorus, spline_state, fixed_double_operator, spline_flow,
                              integrate, mean_value, SECTION4_RUNS, SECTION4_ORDER, entropy,
                              EllipticEnergy
using GeometricBrackets: nbasis, project, evaluate, vectorfield, gradient, entropy_gradient,
                         issymmetric, ispositive_semidefinite, degeneracy_residual,
                         domainvolume, hamiltonian, hessian, field
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
@testset "$(rpad("Spline Solver Tests", 80))" begin
    t = SplineTorus(16, 3)
    N = nbasis(t.space)

    @testset "$(rpad("The space INTEGRATES the domain and projects idempotently", 76))" begin
        @test integrate(t, ones(N))≈4π^2 atol=1e-12
        @test domainvolume(t.space)≈4π^2 atol=1e-12
        @test project(t.space, x -> 2.5)≈fill(2.5, N) rtol=1e-12
        û = project(t.space, x -> cos(x[1]) + 0.4sin(2x[2]))
        @test project(t.space, x -> evaluate(t.space, û, (x[1], x[2])))≈û rtol=1e-11
    end

    @testset "$(rpad("The bordered Poisson solve IS mean-free and consistent", 76))" begin
        ω̂ = project(t.space, x -> cos(x[1]) + 0.4sin(2x[2]))
        φ̂ = t.Λ * ω̂
        @test abs(mean_value(t, φ̂)) < 1e-12
        # The Lagrange multiplier vanishes exactly when the right-hand side is mean-free.
        sol = t.Λ.F \ vcat(t.M * ω̂, 0.0)
        @test abs(sol[N + 1]) < 1e-12
        # -Δ has cos x₁ as an eigenfunction with eigenvalue 1, so φ = ω there.
        e = project(t.space, x -> cos(x[1]))
        @test norm(t.Λ * e .- e) / norm(e) < 1e-3
    end

    # `*` takes a `StridedVector`, to keep the method unambiguous. Anything else reaches the
    # operator through Base's generic `*`, which goes via `mul!` — and the `mul!` method in
    # `spline.jl` is what makes that one solve instead of one solve per entry. Both halves are
    # asserted: the value, and that the method reached is ours. Without the second, a regression
    # is silent, because the `getindex` path returns the same numbers and only costs N² solves.
    @testset "$(rpad("A NON-strided vector solves once, not once per entry", 76))" begin
        r = range(0.25, 0.75; length = N)
        @test !(typeof(r) <: StridedVector)
        @test t.Λ * r ≈ t.Λ * collect(r)
        @test which(mul!, (Vector{Float64}, typeof(t.Λ), typeof(r))).module ===
              MetriplecticRelaxation
    end

    @testset "$(rpad("M*Lambda IS symmetric, as EllipticEnergy assumes", 76))" begin
        MΛ = t.M * Matrix(t.Λ)
        @test norm(MΛ - MΛ') / norm(MΛ) < 1e-11
        Hen = EllipticEnergy(t.Λ, t.M)
        û = randn(N)
        û .-= mean_value(t, û)
        @test hamiltonian(Hen, t.space, û) > 0
        @test_throws ArgumentError hessian(Hen, t.space, û)
    end

    @testset "$(rpad("Every bracket IS degenerate on the flow's own Hamiltonian", 76))" begin
        for name in SECTION4_ORDER
            f = spline_flow(t, SECTION4_RUNS[name])
            ω̂, _ = spline_state(t, SECTION4_RUNS[name])
            ŵ = ω̂ .+ 0.1 .* randn(N)
            ŵ .-= mean_value(t, ŵ)
            @test issymmetric(f.metric, ŵ)
            @test ispositive_semidefinite(f.metric, ŵ)
            @test degeneracy_residual(f, ŵ) < 1e-11
            v = vectorfield(f, ŵ)
            @test abs(dot(gradient(f, ŵ), v)) / (norm(gradient(f, ŵ)) * norm(v)) < 1e-11
            @test dot(entropy_gradient(f, ŵ), v) < 0
        end
    end

    @testset "$(rpad("A1's assembled operator IS the generic vector field", 76))" begin
        spec = SECTION4_RUNS["a1"]
        f = spline_flow(t, spec)
        A = fixed_double_operator(t, spec)
        ŵ = randn(N)
        ŵ .-= mean_value(t, ŵ)
        @test -(t.Mfac \ (A * ŵ))≈vectorfield(f, ŵ) rtol=1e-11
        @test norm(A - A') / norm(A) < 1e-12
        @test norm(A * ones(N)) / norm(A) < 1e-12    # constants are in the kernel
        @test nnz(A) < N^2 / 2
    end
end
