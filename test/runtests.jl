using MetriplecticRelaxation
using MetriplecticRelaxation: SpectralTorus, SplineTorus, Diagnostics, Trace,
                              torus_field, spectral_state, spectral_rhs, spectral_step,
                              spline_state, spline_rhs, spline_step,
                              fixed_double_operator, spline_flow,
                              ∂₁, ∂₂, laplacian, poisson_periodic, canonical_bracket,
                              hamiltonian_field, double_bracket_field, parallel_diffusion,
                              projector_bracket_field, integrate, l2inner, l2norm,
                              mean_value, agm, contour_length, contour_length_quadrature,
                              contour_average, contour_deviation, relaxation_time,
                              islands_h, CENTRAL_ISLANDS, periodise,
                              SECTION4_RUNS, SECTION4_ORDER,
                              euler_minimiser, euler_entropy_minimum,
                              energy, entropy,
                              potential_norm², best_fit_euler, fit_rate, record!,
                              entropy_monotone, cone_residual,
                              EllipticEnergy
using PoissonBrackets: nbasis, project, evaluate, spectral_grid, ∂x, ∂y, vectorfield,
                       gradient, entropy_gradient, issymmetric, ispositive_semidefinite,
                       degeneracy_residual, domainvolume, hamiltonian, hessian
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

@testset "$(rpad("Spectral Solver Tests", 80))" begin
    g = SpectralTorus(24)
    U = torus_field(g, (a, b) -> 0.7cos(a) + 0.4sin(2b) + 0.3cos(a - b) + 1.1)
    H = torus_field(g, (a, b) -> cos(a) + 0.6sin(b) + 0.25cos(a + b))

    @testset "$(rpad("The FFT derivatives ARE PoissonBrackets' matrix ones", 76))" begin
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

@testset "$(rpad("Projector Rate Tests", 80))" begin
    # The linearised spectrum of the projector flow is exactly {1 - 1/λ} over the eigenvalues
    # of -Δ. That is what makes the manuscript's "≈ 1/2" and "≈ 1" exact rather than fitted.
    g = SpectralTorus(32)
    H₀ = 0.35
    ω★ = torus_field(g, euler_minimiser(H₀, 0.6, 0.0, 0.0))
    φ★ = poisson_periodic(g, ω★)

    function jac(v; ε = 1e-6)
        f(w) = projector_bracket_field(g, w, poisson_periodic(g, w))
        (f(ω★ .+ ε .* v) .- f(ω★ .- ε .* v)) ./ 2ε
    end

    @testset "$(rpad("The relaxed state IS a fixed point", 76))" begin
        @test φ★≈ω★ rtol=1e-13
        @test l2norm(g, projector_bracket_field(g, ω★, φ★)) / l2norm(g, ω★) < 1e-13
        @test l2inner(g, ω★, ω★) / 2≈euler_entropy_minimum(H₀) rtol=1e-13
    end

    @testset "$(rpad("An eigenmode of eigenvalue lambda DECAYS at 1 - 1/lambda", 76))" begin
        for (f, λ) in (((a, b) -> cos(a + b), 2.0), ((a, b) -> sin(a - b), 2.0),
            ((a, b) -> cos(2a), 4.0), ((a, b) -> cos(2a + b), 5.0),
            ((a, b) -> cos(3a), 9.0))
            v = torus_field(g, f)
            Jv = jac(v)
            @test -l2inner(g, Jv, v) / l2inner(g, v, v)≈1 - 1 / λ atol=1e-6
            # And J v really is parallel to v, so "the rate" is well defined.
            @test l2norm(g, Jv .+ (1 - 1 / λ) .* v) / l2norm(g, Jv) < 1e-5
        end
    end

    @testset "$(rpad("The whole lambda = 1 eigenspace IS neutral", 76))" begin
        for f in ((a, b) -> sin(a), (a, b) -> cos(b), (a, b) -> sin(b))
            v = torus_field(g, f)
            v .-= (l2inner(g, v, φ★) / l2inner(g, φ★, φ★)) .* φ★
            l2norm(g, v) < 1e-10 && continue
            @test abs(l2inner(g, jac(v), v) / l2inner(g, v, v)) < 1e-6
        end
        # Parallel to φ★ too: scaling moves along the family to another H₀, and every member
        # of it is a fixed point.
        @test abs(l2inner(g, jac(φ★), φ★) / l2inner(g, φ★, φ★)) < 1e-6
    end
end
