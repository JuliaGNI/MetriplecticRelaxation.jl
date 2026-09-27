using MetriplecticRelaxation
using MetriplecticRelaxation: SpectralTorus, torus_field, ∂₁, ∂₂, laplacian,
                              poisson_periodic, canonical_bracket, hamiltonian_field,
                              double_bracket_field, parallel_diffusion,
                              projector_bracket_field, l2inner, l2norm, mean_value
using GeometricBrackets: spectral_grid, ∂x, ∂y
using LinearAlgebra
using Random
using SparseArrays
using Test

Random.seed!(0x5c1e9a3b)

# Every resolution below is deliberately far coarser than a Section 4 run, and every final time
# far shorter. These tests assert the STRUCTURAL claims -- conservation, monotonicity,
# degeneracy, the closed forms -- which hold on any mesh. The claims that need resolution to be
# true at all, the rates and the contour averages, are the run drivers' business and are
# measured in `scripts/`, not here.
@testset "$(rpad("Spectral Solver Tests", 80))" begin
    g = SpectralTorus(24)
    U = torus_field(g, (a, b) -> 0.7cos(a) + 0.4sin(2b) + 0.3cos(a - b) + 1.1)
    H = torus_field(g, (a, b) -> cos(a) + 0.6sin(b) + 0.25cos(a + b))

    @testset "$(rpad("The FFT derivatives ARE GeometricBrackets' matrix ones", 76))" begin
        gp = spectral_grid(24)
        @test ∂₁(g, U)≈∂x(gp, U) atol=1e-12
        @test ∂₂(g, U)≈∂y(gp, U) atol=1e-12
        @test canonical_bracket(g, U, H)≈∂x(gp, U) .* ∂y(gp, H) .- ∂y(gp, U) .* ∂x(gp, H) atol=1e-12
    end

    @testset "$(rpad("The Poisson solve INVERTS the Laplacian, mean-free", 76))" begin
        ω = torus_field(g, (a, b) -> cos(a) + 0.5sin(2b))
        φ = poisson_periodic(g, ω)
        @test φ≈torus_field(g, (a, b) -> cos(a) + 0.5sin(2b) / 4) atol=1e-13
        @test abs(mean_value(g, φ)) < 1e-14
        @test .-laplacian(g, φ)≈ω atol=1e-12
    end

    @testset "$(rpad("The hoisted parallel diffusion IS the nested bracket", 76))" begin
        X = hamiltonian_field(g, H)
        @test parallel_diffusion(g, U, X)≈double_bracket_field(g, U, H) atol=1e-13
        @test maximum(abs, ∂₁(g, X[1]) .+ ∂₂(g, X[2])) < 1e-12   # ∇·X_h = 0
    end

    @testset "$(rpad("Both vector fields CONSERVE H and dissipate S", 76))" begin
        ω = U .- mean_value(g, U)
        # Double bracket, prescribed h.
        f = double_bracket_field(g, ω, H)
        hz = H .- mean_value(g, H)
        @test abs(l2inner(g, hz, f)) / (l2norm(g, hz) * l2norm(g, f)) < 1e-13
        @test l2inner(g, ω, f) < 0
        # Projector bracket.
        φ = poisson_periodic(g, ω)
        fp = projector_bracket_field(g, ω, φ)
        @test abs(l2inner(g, φ, fp)) / (l2norm(g, φ) * l2norm(g, fp)) < 1e-13
        @test l2inner(g, ω, fp) < 0
    end

    @testset "$(rpad("The factor 2 in the projector field IS what eq:SS-projector needs", 76))" begin
        # The manuscript's prose prints H/‖φ‖²; its own eq:SS-projector requires 2H/‖φ‖².
        ω = U .- mean_value(g, U)
        φ = poisson_periodic(g, ω)
        H₀ = l2inner(g, φ, ω) / 2
        S = l2inner(g, ω, ω) / 2
        SS = 2S - 4H₀^2 / l2inner(g, φ, φ)
        @test -l2inner(g, ω, projector_bracket_field(g, ω, φ))≈SS rtol=1e-12
        # With the printed factor it is neither energy-conserving nor consistent with (S,S).
        wrong = .-(ω .- (H₀ / l2inner(g, φ, φ)) .* φ)
        @test abs(l2inner(g, φ, wrong)) / (l2norm(g, φ) * l2norm(g, wrong)) > 1e-2
        @test !isapprox(-l2inner(g, ω, wrong), SS; rtol = 1e-2)
    end
end
