using MetriplecticRelaxation
using MetriplecticRelaxation: SpectralTorus, Diagnostics, Trace, torus_field,
                              spectral_state, l2inner, l2norm, SECTION4_RUNS,
                              euler_minimiser, euler_entropy_minimum, energy,
                              best_fit_euler, fit_rate, record!, entropy_monotone
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
@testset "$(rpad("Diagnostics Tests", 80))" begin
    g = SpectralTorus(24)

    @testset "$(rpad("best_fit_euler IS the minimiser over the three phases", 76))" begin
        spec = SECTION4_RUNS["a3"]
        d = Diagnostics(g, spec)
        ω, _ = spectral_state(g, spec)
        H₀ = energy(d, ω)
        fit, res, _ = best_fit_euler(d, ω, H₀)
        # The fit is a member of the family, so its entropy is S_η exactly.
        @test l2inner(g, fit, fit) / 2≈euler_entropy_minimum(H₀) rtol=1e-12
        # A coarse scan cannot beat it.
        best = Inf
        for i in 0:23, j in 0:23, k in 0:23
            w = torus_field(g, euler_minimiser(H₀, 2π * i / 24, 2π * j / 24, 2π * k / 24))
            best = min(best, l2inner(g, ω .- w, ω .- w))
        end
        @test res <= sqrt(best) + 1e-12
        # A state already in the family is fitted to round-off.
        w = torus_field(g, euler_minimiser(H₀, 0.7, 1.3, 2.9))
        @test best_fit_euler(d, w, H₀)[2] / l2norm(g, w) < 1e-13
    end

    @testset "$(rpad("fit_rate RECOVERS a known exponential and flags one that is not", 76))" begin
        t = collect(0.0:0.01:10.0)
        for λ in (0.5, 1.0, 2.0)
            (r, r², _) = fit_rate(t, 3.7 .* exp.(-λ .* t))
            @test r≈λ atol=1e-10
            @test r² > 0.9999999
        end
        @test fit_rate(t, 1 ./ (1 .+ t) .^ 2)[2] < 0.999
        # A round-off tail must not drag the slope.
        @test fit_rate(t, max.(3.7 .* exp.(-t), 1e-17))[1]≈1.0 atol=1e-8
    end

    @testset "$(rpad("entropy_monotone CATCHES a rising entropy", 76))" begin
        spec = SECTION4_RUNS["a3"]
        d = Diagnostics(g, spec)
        ω, _ = spectral_state(g, spec)
        tr = Trace(ω)
        for i in 1:8
            record!(tr, d, 0.1i, (1 - 0.01i) .* ω)
        end
        @test entropy_monotone(tr)[1]
        tr2 = Trace(ω)
        for i in 1:8
            record!(tr2, d, 0.1i, (1 + 0.01i) .* ω)
        end
        @test !entropy_monotone(tr2)[1]
    end
end
