using MetriplecticRelaxation
using MetriplecticRelaxation: entropy
using MetriplecticRelaxation: SECTION5_RUNS, SECTION5_ORDER, EulerSpec, EulerSquare,
                              B3_FLOOR, euler_state, state_extrema, gaussian_w2,
                              perturbation_b2, DIRICHLET_EIGENVALUE
using GeometricBrackets: project
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
# =============================================================================================
# Section 5.4: reduced Euler on [0,1]² with homogeneous Dirichlet conditions.
#
# The mesh here is 10 cells against the runs' 26, and every final time is a handful of steps.
# What is asserted is again the STRUCTURE -- the transcription of the initial conditions, the
# space, the bracket's three defining properties, the analytic entropy derivatives and the
# algebra behind the two closed-form references -- all of which hold on any mesh.
#
# One deliberate exception to "coarser is fine": B3's own initial state is NOT admissible on a
# coarse mesh, because the L² projection of its narrow Gaussian undershoots below zero and
# `y log y` is undefined there. So the Gibbs tests use an admissible state of their own, and the
# threshold mesh itself is a scripts/ measurement rather than a test -- 26 cells is 10 s per
# implicit step.

@testset "$(rpad("Section 5 Problem Tests", 80))" begin
    @testset "$(rpad("Section 5 states w SQUARED, not w", 76))" begin
        for name in SECTION5_ORDER
            g = SECTION5_RUNS[name].gaussian
            @test g.w[1]^2≈0.01 rtol=1e-14
            @test g.w[2]^2≈0.07 rtol=1e-14
        end
        g = SECTION5_RUNS["b1"].gaussian
        @test g(0.5, 0.5) == 1.0
        # The e-folding distance is √(w²) = 0.1 and 0.2646, not 0.01 and 0.07.
        @test g(0.5 + 0.1, 0.5) / g(0.5, 0.5)≈exp(-1) rtol=1e-14
        @test g(0.5, 0.5 + sqrt(0.07)) / g(0.5, 0.5)≈exp(-1) rtol=1e-14
        # The control: reading the printed number as w itself is wrong by order one.
        @test !isapprox(g(0.5 + 0.01, 0.5) / g(0.5, 0.5), exp(-1); rtol = 0.1)
        @test !isapprox(g(0.5, 0.5 + 0.07) / g(0.5, 0.5), exp(-1); rtol = 0.1)
        # B2 states N and B3 states 1/N, which is the other easy transcription error.
        @test SECTION5_RUNS["b2"].gaussian(0.5, 0.5)≈0.01 rtol=1e-14
        @test SECTION5_RUNS["b3"].gaussian(0.5, 0.5)≈10.0 rtol=1e-14
        @test all(isapprox.(gaussian_w2((0.25, 0.75), (0.04, 0.09), 2.0).w, (0.2, 0.3);
            rtol = 1e-14))
    end

    @testset "$(rpad("B3 carries a POSITIVE FLOOR, and it is not decoration", 76))" begin
        # `s = y log y` is undefined at y ≤ 0 and eq:M-condition gives it the mobility M = y,
        # which must be positive. The printed Gaussian decays to 1.4e-11 of its peak at the far
        # corner — strictly positive as a function, indistinguishable from zero to a Galerkin
        # scheme. The floor is what gives the admissible set an interior.
        b3 = SECTION5_RUNS["b3"]
        @test b3.entropy === :gibbs
        @test b3.background !== nothing
        @test b3.background(0.13, 0.87) == B3_FLOOR
        @test B3_FLOOR / b3.gaussian(0.5, 0.5)≈0.01 rtol=1e-14
        # The printed condition at the far corner, which is what the floor has to dominate.
        @test b3.gaussian(0.0, 0.0) / b3.gaussian(0.5, 0.5) < 1e-10
        # And the floored condition is admissible on a mesh where the printed one is not.
        sq = EulerSquare(12, 2)
        printed = EulerSpec("printed", b3.section, b3.entropy, b3.state, b3.gaussian,
            nothing, b3.Δt, b3.T)
        @test state_extrema(sq, euler_state(sq, printed))[1] < 0
        @test state_extrema(sq, euler_state(sq, b3))[1] > 1e-3
    end

    @testset "$(rpad("Every run uses the homogeneous-Dirichlet state space", 76))" begin
        # It is the absence of the constants that forces the equilibrium's mass multiplier to
        # zero, and both of §5.4's closed-form references are the μ = 0 case.
        for name in SECTION5_ORDER
            @test SECTION5_RUNS[name].state === :dirichlet
        end
    end

    @testset "$(rpad("B2's added mode IS an unstable equilibrium", 76))" begin
        # An eigenfunction of -Δ satisfies ω = λφ, hence ∇(δS/δu) ∥ ∇(δH/δu): a stationary
        # state of the relaxation. 52π² is not the lowest eigenvalue, so it is unstable.
        ε = 1e-5
        f = perturbation_b2
        x₁, x₂ = 0.31, 0.57
        lap = (f(x₁ + ε, x₂) - 2f(x₁, x₂) + f(x₁ - ε, x₂)) / ε^2 +
              (f(x₁, x₂ + ε) - 2f(x₁, x₂) + f(x₁, x₂ - ε)) / ε^2
        @test -lap / f(x₁, x₂)≈52π^2 rtol=1e-6
        @test 52π^2 > DIRICHLET_EIGENVALUE
        for p in ((0.0, 0.37), (1.0, 0.37), (0.42, 0.0), (0.42, 1.0))
            @test abs(f(p...)) < 1e-15
        end
    end
end
