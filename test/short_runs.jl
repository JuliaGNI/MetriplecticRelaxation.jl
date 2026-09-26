using MetriplecticRelaxation
using MetriplecticRelaxation: SpectralTorus, SplineTorus, Diagnostics, Trace,
                              spectral_state, spectral_rhs, spectral_step, spline_state,
                              spline_rhs, spline_step, integrate, SECTION4_RUNS,
                              SECTION4_ORDER, energy, entropy, potential_norm², record!,
                              cone_residual
using GeometricBrackets: project, field
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
@testset "$(rpad("Short Run Tests", 80))" begin
    # Twenty steps at a coarse resolution. Far too short to relax anything, which is the point:
    # what is asserted here is that the STRUCTURE survives time stepping, and that holds from
    # the first step. The relaxation claims themselves are the run drivers' business.
    # 48 rather than 24 points, and A2 is the only reason. Its energy conservation rests on an
    # integration by parts of a CUBIC nonlinearity -- ∫φ ∇·(X_φ ⊗ X_φ ∇ω), with X_φ built from
    # φ = φ(ω) -- and the product of band-limited fields is not band-limited, so the trapezoidal
    # rule that makes the by-parts exact for A1 aliases here. The error is spatial and converges
    # spectrally: measured, it is 4.7e-6 at 24 points and 7.5e-12 at 48, and INDEPENDENT of Δt,
    # which is what identifies it as aliasing rather than as time-integration error.
    #
    # A1, A3 and A4 conserve H to round-off at 24 points already: A1 because its H is linear in
    # ω with a fixed generating field, so every Runge-Kutta stage lies in the kernel of the
    # exactly-degenerate operator, and A3/A4 because the projector bracket needs no quadrature
    # at all -- it is algebra in φ̂.
    @testset "$(rpad("Both discretisations CONSERVE H over a short run", 76))" begin
        for name in SECTION4_ORDER
            spec = SECTION4_RUNS[name]
            g = SpectralTorus(48)
            dg = Diagnostics(g, spec)
            ω, _ = spectral_state(g, spec)
            rhs = spectral_rhs(g, spec)
            H₀, S₀ = energy(dg, ω), entropy(dg, ω)
            for _ in 1:20
                ω = spectral_step(rhs, ω, spec.Δt)
            end
            @test abs(energy(dg, ω) - H₀) / abs(H₀) < 1e-10
            @test entropy(dg, ω) <= S₀
            @test abs(integrate(g, ω)) < 1e-12

            t = SplineTorus(16, 3)
            dt = Diagnostics(t, spec)
            ω̂, _ = spline_state(t, spec)
            rhs = spline_rhs(t, spec)
            Ĥ₀, Ŝ₀ = energy(dt, ω̂), entropy(dt, ω̂)
            for _ in 1:20
                ω̂ = spline_step(rhs, ω̂, spec.Δt)
            end
            @test abs(energy(dt, ω̂) - Ĥ₀) / abs(Ĥ₀) < 1e-10
            @test entropy(dt, ω̂) <= Ŝ₀
            @test abs(integrate(t, ω̂)) < 1e-12
        end
    end

    @testset "$(rpad("The cone bounds HOLD for the initial state of every run", 76))" begin
        # eq:theoretical-limits constrain any state of energy H₀, so they must hold at t = 0.
        g = SpectralTorus(24)
        for name in ("a2", "a3", "a4")
            spec = SECTION4_RUNS[name]
            d = Diagnostics(g, spec)
            ω, _ = spectral_state(g, spec)
            H₀ = energy(d, ω)
            @test entropy(d, ω) >= H₀                          # eq:tb2
            @test potential_norm²(d, ω) <= 2H₀ + 1e-12          # eq:tb3
            tr = Trace(ω)
            record!(tr, d, 0.0, ω)
            @test all(<(1e-9), cone_residual(tr, H₀))
        end
    end
end
