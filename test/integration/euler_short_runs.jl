using MetriplecticRelaxation
using MetriplecticRelaxation: SpectralTorus, Diagnostics, Trace, spectral_state,
                              SECTION4_RUNS, energy, entropy, record!
using MetriplecticRelaxation: SECTION5_RUNS, EulerSquare, GibbsEntropy, euler_state,
                              euler_flow, euler_entropy_floor, dirichlet_eigenvalue,
                              state_extrema, entropy_plateau
using GeometricBrackets: project, hamiltonian, Integrator, ImplicitMidpoint, integrate_step!
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
@testset "$(rpad("Section 5 Short Run Tests", 80))" begin
    @testset "$(rpad("Crank-Nicolson CONSERVES H and dissipates S", 76))" begin
        # ImplicitMidpoint IS Crank-Nicolson for this field. Four steps is enough: H is
        # conserved because the CONVERGED midpoint increment lies in the range of 𝔾(ū), which
        # holds from the first step, and the relaxation itself is scripts/ business.
        for name in ("b1", "b2")
            spec = SECTION5_RUNS[name]
            sq = EulerSquare(10, 2)
            d = Diagnostics(sq, spec)
            f = euler_flow(sq, spec)
            ω̂ = euler_state(sq, spec)
            H₀, S₀ = energy(d, ω̂), entropy(d, ω̂)
            integ = Integrator(f, ImplicitMidpoint(), 0.5; û₀ = ω̂)
            for _ in 1:4
                integrate_step!(ω̂, integ)
            end
            @test abs(energy(d, ω̂) - H₀) / abs(H₀) < 1e-13
            @test entropy(d, ω̂) <= S₀
            # The Poincaré floor holds at every state, so it holds here.
            @test entropy(d, ω̂) >= euler_entropy_floor(H₀; λ = dirichlet_eigenvalue(sq))
        end
    end

    @testset "$(rpad("The Gibbs flow CONSERVES H on an admissible state", 76))" begin
        # B3's own initial state needs 26 cells to be admissible at all, which is 14 s per
        # step; the structural claim does not, so it is made on a state that is positive here.
        spec = SECTION5_RUNS["b3"]
        sq = EulerSquare(10, 2; state = spec.state)
        d = Diagnostics(sq, spec)
        f = euler_flow(sq, spec)
        ω̂ = project(sq.space,
            x -> 6sin(π * x[1]) * sin(π * x[2]) * (1.5 + 0.4sin(2π * x[1])) + 0.4)
        @test state_extrema(sq, ω̂)[1] > 0
        H₀, S₀ = energy(d, ω̂), entropy(d, ω̂)
        @test S₀≈hamiltonian(GibbsEntropy(), sq.space, ω̂) rtol=1e-14
        integ = Integrator(f, ImplicitMidpoint(), 1e-3; û₀ = ω̂)
        for _ in 1:4
            integrate_step!(ω̂, integ)
        end
        @test abs(energy(d, ω̂) - H₀) / abs(H₀) < 1e-12
        @test entropy(d, ω̂) <= S₀
    end

    @testset "$(rpad("entropy_plateau SEPARATES a plateau from an immediate decay", 76))" begin
        # The diagnostic B2's claim rests on, and the control that it is not vacuous. The
        # entropy of `c·ω` is `c²S₀`, so a factor sequence is a chosen entropy trace.
        g = SpectralTorus(16)
        spec = SECTION4_RUNS["a3"]
        d = Diagnostics(g, spec)
        ω, _ = spectral_state(g, spec)
        function trace_of(factors)
            tr = Trace(ω)
            for (i, c) in enumerate(factors)
                record!(tr, d, 0.5i, c .* ω)
            end
            return tr
        end
        flat = trace_of([i <= 10 ? 1.0 : 1.0 - 0.05 * (i - 10) for i in 1:20])
        early = trace_of([1.0 - 0.02i for i in 1:20])
        @test entropy_plateau(flat)[2] > 8
        @test entropy_plateau(flat)[3] < 1e-12
        @test entropy_plateau(early)[2] <= 2
        # A trace that never falls has no break to find.
        @test entropy_plateau(trace_of(ones(5)))[2] == 0
    end
end
