using MetriplecticRelaxation
using MetriplecticRelaxation: agm, contour_length, contour_length_quadrature,
                              contour_average, contour_deviation, relaxation_time,
                              islands_h, CENTRAL_ISLANDS, periodise, SECTION4_RUNS
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
@testset "$(rpad("Torus Geometry Tests", 80))" begin
    @testset "$(rpad("The closed form for l_h IS the arclength integral", 76))" begin
        for h in (0.9, 0.5, 0.1, 0.01, 1e-4)
            @test contour_length(h)≈contour_length_quadrature(h; n = 400) rtol=1e-12
        end
        # The control: a wrong closed form must disagree.
        @test !isapprox(π / (sqrt(0.5) * agm(1.0, 0.5)),
            contour_length_quadrature(0.5; n = 400); rtol = 1e-3)
    end

    @testset "$(rpad("tau_h has the two limits the manuscript states", 76))" begin
        # At the island centre the harmonic approximation h ≈ 1 - s² - t² rotates at angular
        # frequency 2, so the period is π and τ = 1/4.
        @test contour_length(1.0)≈π atol=1e-14
        @test relaxation_time(1.0)≈0.25 atol=1e-14
        # At the separatrix it diverges.
        @test relaxation_time(1e-8) > 1e6
        @test relaxation_time(0.9) < relaxation_time(0.1) < relaxation_time(0.01)
    end

    @testset "$(rpad("A field that is a function of h IS its own contour average", 76))" begin
        for h in (0.8, 0.3, 0.05)
            f(x₁, x₂) = 3islands_h(x₁, x₂)^2 - islands_h(x₁, x₂) + 1
            @test contour_average(f, h, CENTRAL_ISLANDS[1]; n = 400)≈3h^2 - h + 1 atol=1e-11
            @test contour_deviation(f, h, CENTRAL_ISLANDS[1]; n = 400) < 1e-11
        end
    end

    @testset "$(rpad("A4's Gaussian IS discontinuous where A1's and A2's are not", 76))" begin
        # The value the printed formula still has at the far edge of the domain.
        # Measured maxima over the whole boundary: A4 reaches 0.153, which is 8.5 % of its own
        # peak; A2/A3 reach 5.2e-5; A1 reaches 1.2e-25.
        @test SECTION4_RUNS["a4"].gaussian(π, 2π) > 0.15
        @test SECTION4_RUNS["a2"].gaussian(π, 0.0) < 1e-4
        @test SECTION4_RUNS["a1"].gaussian(π, 2π) < 1e-20
        # Periodising changes A4 and leaves the others alone.
        @test periodise(SECTION4_RUNS["a4"].gaussian)(π, 0.0) > 0.15
        @test periodise(SECTION4_RUNS["a1"].gaussian)(π, 2π) < 1e-20
    end
end
