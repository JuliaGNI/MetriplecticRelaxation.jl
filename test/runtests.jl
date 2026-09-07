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
using MetriplecticRelaxation: SECTION5_RUNS, SECTION5_ORDER, EulerSpec, EulerSquare,
                              GibbsEntropy, B3_FLOOR,
                              euler_state, euler_flow, euler_entropy_floor, eigenmode_fit,
                              dirichlet_eigenvalue, gibbs_lambda,
                              interior_weights, state_extrema, gaussian_w2, perturbation_b2,
                              entropy_plateau, DIRICHLET_EIGENVALUE
using MetriplecticRelaxation: SECTION55_RUNS, GradShafranovBox, gs_state, gs_flow, gs_fit,
                              gs_rayleigh, gs_current, gs_density, gs_stiffness,
                              gs_eigenvalue, separable_eigenvalue, TakedaGrid,
                              takeda_iterate, takeda_eigenvalue, herrnegger_mobility,
                              herrnegger_profile, GS_RADIAL, GS_AXIAL, HERRNEGGER_C,
                              HERRNEGGER_D, GS_LAMBDA_RECTANGLE, GS_LAMBDA_CONTINUUM,
                              disk_map, DiskTriangulation, disk_area, disk_eigenvalue,
                              GS_LAMBDA_DISK, GS_LAMBDA_DISK_CONTINUUM
using PoissonBrackets: nbasis, project, evaluate, spectral_grid, ∂x, ∂y, vectorfield,
                       gradient, entropy_gradient, issymmetric, ispositive_semidefinite,
                       degeneracy_residual, domainvolume, hamiltonian, hessian,
                       CollisionBracket, MetriplecticFlow, QuadraticHamiltonian,
                       Integrator, ImplicitMidpoint, integrate_step!, entropy_production,
                       metric_matrix, metric_apply, stiffness_matrix, field,
                       quadrature_weights, quadrature_nodes, basis_values, weighted_matrix,
                       mass_matrix
using SimpleSplines: UniformMesh, Dirichlet
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
end

@testset "$(rpad("Free State Space Tests", 80))" begin
    # B3's space: ω unconstrained on ∂Ω, φ still the homogeneous-Dirichlet solve, reached
    # through the recombination matrix. That split is the finding §5.4 turned up.
    sqf = EulerSquare(12, 2; state = :free)

    @testset "$(rpad("The CONSTANT is in the free space and phi still vanishes", 76))" begin
        c = project(sqf.space, x -> 1.0)
        @test evaluate(sqf.space, c, (0.0, 0.5))≈1.0 atol=1e-12
        φ̂ = sqf.Λ * c
        for p in ((0.0, 0.4), (1.0, 0.4), (0.4, 0.0), (0.4, 1.0))
            @test abs(evaluate(sqf.space, φ̂, p)) < 1e-14
        end
        # -Δφ = 1 on the unit square: φ(½,½) is the classical torsion constant 0.07367135…,
        # reached to 2e-6 at these 12 cells and to 7e-7 at the 16 of `verify_euler.jl`.
        @test evaluate(sqf.space, φ̂, (0.5, 0.5))≈0.0736713532 atol=1e-5
    end

    @testset "$(rpad("M Lambda IS symmetric positive semi-definite there too", 76))" begin
        @test norm(sqf.MΛ - sqf.MΛ') / norm(sqf.MΛ) < 1e-14
        @test minimum(eigvals(Symmetric(sqf.MΛ))) > -1e-14
        @test dirichlet_eigenvalue(sqf)≈DIRICHLET_EIGENVALUE rtol=1e-3
        @test dirichlet_eigenvalue(sqf) >= DIRICHLET_EIGENVALUE
    end

    @testset "$(rpad("A bad state space name IS rejected", 76))" begin
        @test_throws ArgumentError EulerSquare(4, 2; state = :neumann)
    end
end

@testset "$(rpad("Collision Bracket Tests", 80))" begin
    @testset "$(rpad("Every run's bracket IS symmetric, PSD and DEGENERATE", 76))" begin
        for name in SECTION5_ORDER
            spec = SECTION5_RUNS[name]
            # Each run's own space, because that is part of the run.
            sq = EulerSquare(10, 2; state = spec.state)
            # Positive at every quadrature node, so B3's mobility M = ω is admissible, and not
            # an eigenmode, so the vector field is not accidentally zero. The constant offset
            # keeps it away from zero at the boundary and survives only in the free space.
            ω̂ = project(sq.space,
                x -> sin(π * x[1]) * sin(π * x[2]) *
                     (1.5 + 0.4sin(2π * x[1]) - 0.3cos(π * x[2])) +
                     (spec.state === :free ? 0.4 : 0.0))
            @test state_extrema(sq, ω̂)[1] > 0
            f = euler_flow(sq, spec)
            @test issymmetric(f.metric, ω̂)
            @test ispositive_semidefinite(f.metric, ω̂)
            # The flow's form, not the bracket's: the bracket supplies its own generator and is
            # clean by construction even when it is degenerate on the wrong energy.
            @test degeneracy_residual(f, ω̂) < 1e-11
            v = vectorfield(f, ω̂)
            @test abs(dot(gradient(f, ω̂), v)) / (norm(gradient(f, ω̂)) * norm(v)) < 1e-11
            @test dot(entropy_gradient(f, ω̂), v) < 0
            @test entropy_production(f, ω̂) > 0
            @test entropy_production(f, ω̂)≈-dot(entropy_gradient(f, ω̂), v) rtol=1e-12
            # Two independent evaluations of the same operator: the 𝔽_s moments against the
            # R/S/T factorisation.
            c = Vector(sq.M * ω̂)
            @test metric_apply(f.metric, ω̂, c)≈metric_matrix(f.metric, ω̂) * c rtol=1e-10
        end
    end

    sq = EulerSquare(10, 2)
    ω̂ = project(sq.space,
        x -> sin(π * x[1]) * sin(π * x[2]) * (1.5 + 0.4sin(2π * x[1]) - 0.3cos(π * x[2])))

    @testset "$(rpad("A MISMATCHED generator breaks the flow's degeneracy", 76))" begin
        # The trap `MetriplecticFlow` documents: a bracket generated by a different elliptic
        # problem is a perfectly good metric bracket, degenerate on the wrong energy.
        # A prescribed field unrelated to φ, not another elliptic solve: (K+M)⁻¹M gives a
        # stream function 0.95 φ, and degeneracy sees only direction.
        ĥ = project(sq.space, x -> sin(2π * x[1]) * sin(3π * x[2]))
        G = CollisionBracket(sq.space, ĥ)
        fw = MetriplecticFlow(sq.space, G, QuadraticHamiltonian(sq.MΛ),
            QuadraticHamiltonian(Matrix(sq.M)))
        @test degeneracy_residual(G, ω̂) < 1e-11
        # The ratio, not an absolute floor: `degeneracy_residual` normalises by max|𝔾| max|g|,
        # which bounds one product rather than the N-term sum, so a complete misalignment
        # normalises to 1e-3 -- thirteen orders above the clean value, and that is the gap.
        @test degeneracy_residual(fw, ω̂) > 1e8 * degeneracy_residual(G, ω̂)
    end

    @testset "$(rpad("An INDEFINITE mobility loses positivity and nothing else", 76))" begin
        # eq:M-condition requires M > 0. A mobility that changes sign keeps the symmetry and
        # the degeneracy, which is why positivity is a separate assertion.
        #
        # x₁ - 0.3, not x₁ - 0.5: the latter integrates to zero over the square, so the
        # bracket's moment m₀ = ∫M dμ vanishes and its recentring divides by it -- which
        # destroys the degeneracy too, for a reason that has nothing to do with the sign.
        Gi = CollisionBracket(sq.space, sq.Λ; mobility = (x, u) -> x[1] - 0.3,
            mobility_derivative = 0)
        @test !ispositive_semidefinite(Gi, ω̂)
        @test issymmetric(Gi, ω̂)
        @test degeneracy_residual(Gi, ω̂) < 1e-11
    end
end

@testset "$(rpad("Gibbs Entropy Tests", 80))" begin
    # The free space, because that is where `y log y` is usable at all: a state bounded away
    # from zero needs a nonzero boundary trace.
    sq = EulerSquare(10, 2; state = :free)
    s = sq.space
    ω̂ = project(s, x -> sin(π * x[1]) * sin(π * x[2]) * (1.5 + 0.4sin(2π * x[1])) + 0.4)

    @testset "$(rpad("The analytic gradient and Hessian ARE the derivatives", 76))" begin
        @test state_extrema(sq, ω̂)[1] > 0
        for H in (QuadraticHamiltonian(Matrix(sq.M)), GibbsEntropy())
            g = gradient(H, s, ω̂)
            h = hessian(H, s, ω̂)
            for _ in 1:3
                v = randn(length(ω̂))
                v ./= norm(v)
                ε = 1e-6 * norm(ω̂)
                dS = (hamiltonian(H, s, ω̂ .+ ε .* v) -
                      hamiltonian(H, s, ω̂ .- ε .* v)) / 2ε
                @test dS≈dot(g, v) rtol=1e-6
                dg = (gradient(H, s, ω̂ .+ ε .* v) .- gradient(H, s, ω̂ .- ε .* v)) ./ 2ε
                @test norm(dg .- h * v) / norm(h * v) < 1e-5
            end
        end
    end

    @testset "$(rpad("The Hessian of y log y IS 1/M, which is eq:M-condition", 76))" begin
        # M ∂²_y s = 1 with s = y log y gives M = y, so the Hessian's weight is 1/ω exactly.
        u = field(s, ω̂, (0, 0))
        Φ = basis_values(s, (0, 0))
        @test hessian(GibbsEntropy(), s, ω̂) ≈
              Matrix(Φ * Diagonal(quadrature_weights(s) ./ u) * Φ') rtol=1e-14
        # And the premise that makes the weight 1/u finite in the first place: ω > 0 at every
        # quadrature node. `u .* (1 ./ u) ≈ 1` stood here, which is true of any nonzero float
        # and so asserted nothing.
        @test all(>(0), u)
    end

    @testset "$(rpad("y log y REFUSES a state it does not admit", 76))" begin
        # No guard: `log` raises, which is the correct outcome for a state the entropy is not
        # defined on. A guard would turn it into a silently wrong number.
        @test_throws DomainError hamiltonian(GibbsEntropy(), s, .-ω̂)
    end
end

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

@testset "$(rpad("Grad-Shafranov Problem Tests", 80))" begin
    @testset "$(rpad("Section 5.5 states w SQUARED too, and its own values", 76))" begin
        g = SECTION55_RUNS["c1"].gaussian
        @test g.w[1]^2≈0.5 rtol=1e-14
        @test g.w[2]^2≈3.2 rtol=1e-14
        @test g(4.0, 0.0) == 1.0
        @test g(4.0 + sqrt(0.5), 0.0) / g(4.0, 0.0)≈exp(-1) rtol=1e-14
        @test g(4.0, sqrt(3.2)) / g(4.0, 0.0)≈exp(-1) rtol=1e-14
        # The control: the printed numbers read as w rather than w² are wrong by order one.
        @test !isapprox(g(4.0 + 0.5, 0.0) / g(4.0, 0.0), exp(-1); rtol = 0.1)
        @test !isapprox(g(4.0, 3.2) / g(4.0, 0.0), exp(-1); rtol = 0.1)
    end

    @testset "$(rpad("f(r,·) INVERTS the entropy derivative, and the current form does not", 76))" begin
        # The manuscript defines f as ∂_y s(r,·)⁻¹, and for s = y²/2(Cr²+D) that inverse is
        # (Cr²+D)y — hence eq:Grad-Shafranov-equation is a LINEAR eigenvalue problem here.
        ∂ys(r, y) = y / (HERRNEGGER_C * r^2 + HERRNEGGER_D)
        for r in range(GS_RADIAL...; length = 5), y in (-2.0, -0.3, 0.7, 3.0)

            @test ∂ys(r, herrnegger_profile(r, y))≈y rtol=1e-14
            @test herrnegger_mobility(r) * (1 / (HERRNEGGER_C * r^2 + HERRNEGGER_D))≈1 rtol=1e-14
        end
        # The paper also prints the CURRENT of the same equilibrium, λ(Cr + D/r)ψ. Reading that
        # as f is off by a factor of r, which on [1,7] is not a perturbation.
        @test !isapprox(∂ys(5.0, (HERRNEGGER_C * 5.0 + HERRNEGGER_D / 5.0) * 1.0), 1.0;
            rtol = 0.1)
    end

    @testset "$(rpad("The measure is dmu = dr dz / r, in Delta-star and in the profile", 76))" begin
        s = MetriplecticRelaxation.TensorSplineSpace(
            (UniformMesh(8, GS_RADIAL), UniformMesh(9, GS_AXIAL)), 2, Dirichlet())
        # ∫ dμ over Ω is (b−a in z) log(7), which the weighted quadrature has to reproduce. The
        # tolerance is the quadrature's own error on a rational integrand: 3.2e-6 at 8 cells,
        # 7.1e-8 at 16 and 1.1e-10 at 48, so this is not a slack bound but a measured one.
        wμ = quadrature_weights(s) .* gs_density.(quadrature_nodes(s))
        @test sum(wμ)≈(GS_AXIAL[2] - GS_AXIAL[1]) * log(7.0) rtol=1e-5
        # gs_stiffness is symmetric positive definite and is NOT the plain stiffness matrix.
        Kμ = Matrix(gs_stiffness(s))
        @test maximum(abs, Kμ - Kμ') / maximum(abs, Kμ) < 1e-14
        @test minimum(eigvals(Symmetric(Kμ))) > 0
        @test maximum(abs, Kμ - Matrix(stiffness_matrix(s))) / maximum(abs, Kμ) > 0.1
    end

    @testset "$(rpad("The eigenvalue is 0.0302346260, by three independent routes", 76))" begin
        # Separation of variables on the rectangle, a 2D Galerkin solve that uses none, and the
        # classical finite-volume iteration. Coarse here; `verify_takeda.jl` refines.
        sep = separable_eigenvalue(; cells = 32, degree = 3).λ
        @test sep≈GS_LAMBDA_CONTINUUM rtol=1e-8
        s = MetriplecticRelaxation.TensorSplineSpace(
            (UniformMesh(12, GS_RADIAL), UniformMesh(14, GS_AXIAL)), 3, Dirichlet())
        @test gs_eigenvalue(s)≈GS_LAMBDA_CONTINUUM rtol=1e-6
        g = TakedaGrid(24, 24)
        it = takeda_iterate(g)
        @test it.converged
        @test it.λ≈takeda_eigenvalue(g) rtol=1e-10
        @test it.λ≈GS_LAMBDA_CONTINUUM rtol=1e-2
        @test it.λ < GS_LAMBDA_CONTINUUM          # the finite-volume value comes from below
        # The printed λ = 0.030302 is a DISCRETE number, 0.22 % above the continuum one.
        @test GS_LAMBDA_RECTANGLE > GS_LAMBDA_CONTINUUM
        @test abs(GS_LAMBDA_RECTANGLE - GS_LAMBDA_CONTINUUM) / GS_LAMBDA_CONTINUUM≈2.2e-3 rtol=0.1
    end

    @testset "$(rpad("The state j = u/r is what makes the degeneracy EXACT", 76))" begin
        box = GradShafranovBox((6, 7), 2)
        f = gs_flow(box)
        ĵ = gs_state(box, SECTION55_RUNS["c1"])
        s = box.space
        # Λᵀ K^μ = M is the identity the whole formulation rests on, and with it
        # ∂H/∂ĵ = M ψ̂ and 𝔾 ∂H/∂ĵ = 0.
        @test maximum(abs, box.Λ' * Matrix(gs_stiffness(s)) - Matrix(mass_matrix(s))) /
              maximum(abs, Matrix(mass_matrix(s))) < 1e-10
        @test degeneracy_residual(f, ĵ) < 1e-12
        @test degeneracy_residual(f, randn(nbasis(box))) < 1e-12
        @test issymmetric(f.metric, ĵ)
        @test ispositive_semidefinite(f.metric, ĵ)
        @test entropy_production(f, ĵ) > 0
        # H is the poloidal magnetic energy and S the manuscript's own entropy, both checked
        # against their integrals rather than against the forms that compute them.
        ψ̂ = box.Λ * ĵ
        @test hamiltonian(f, ĵ)≈dot(ψ̂, gs_stiffness(s), ψ̂) / 2 rtol=1e-12
        wμ = quadrature_weights(s) .* gs_density.(quadrature_nodes(s))
        u = gs_current(box, ĵ)
        r = [p[1] for p in quadrature_nodes(s)]
        @test entropy(f, ĵ)≈dot(wμ, u .^ 2 ./ (2 .* herrnegger_mobility.(r))) rtol=1e-12
        # The mass Casimir is ∫u dμ = ∫j dx, the cancellation the state variable is chosen for.
        @test integrate(box, ĵ)≈dot(wμ, u) rtol=1e-12
    end

    @testset "$(rpad("CONTROL: the u formulation LOSES the degeneracy", 76))" begin
        # Λ_u = (K^μ)⁻¹M^μ, so ∂H/∂û = M^μψ̂ and the bracket's plain M⁻¹ sandwich no longer
        # recovers the generating field. This is why the state is j and not u.
        box = GradShafranovBox((6, 7), 2)
        s = box.space
        ĵ = gs_state(box, SECTION55_RUNS["c1"])
        Mμ = Matrix(weighted_matrix(s, gs_density.(quadrature_nodes(s)), (0, 0), (0, 0)))
        Λu = Symmetric(Matrix(gs_stiffness(s))) \ Mμ
        A = Mμ * Λu
        fu = MetriplecticFlow(s,
            CollisionBracket(s, Λu; mobility = (p, y) -> herrnegger_mobility(p[1]),
                mobility_derivative = 0, density = gs_density),
            QuadraticHamiltonian(Matrix((A .+ A') ./ 2)), QuadraticHamiltonian(box.W))
        @test degeneracy_residual(fu, ĵ) > 1e-4
    end

    @testset "$(rpad("The relaxed state SATISFIES u/(Cr^2+D) = lambda psi", 76))" begin
        box = GradShafranovBox((10, 12), 2)
        λh = gs_eigenvalue(box)
        @test λh≈GS_LAMBDA_CONTINUUM rtol=1e-4
        @test λh≈gs_eigenvalue(box.space) rtol=1e-6
        # The discrete equilibrium is the lowest eigenvector of (MΛ, W); its scatter residual is
        # the L² projection error of σj and falls with the mesh rather than being zero.
        F = eigen(Symmetric(box.MΛ), Symmetric(box.W))
        v = F.vectors[:, argmax(real(F.values))]
        (λ, _, rel) = gs_fit(box, v)
        @test λ≈λh rtol=1e-6
        @test rel < 1e-3
        @test gs_rayleigh(box, v)≈λh rtol=1e-12
        # The Poincaré floor holds off the flow, not only along it.
        for _ in 1:20
            @test gs_rayleigh(box, randn(nbasis(box))) >= λh
        end
        @test gs_rayleigh(box, gs_state(box, SECTION55_RUNS["c1"])) > λh
    end

    @testset "$(rpad("CONTROL: the dx reading of eq:entropy-2D MISSES eq:gs-ref", 76))" begin
        # eq:entropy-2D writes ∫s dx and eq:GradShafranov-S-H writes ∫s dμ. In the j variable
        # those are two entropy weights, σ = r/(Cr²+D) and σ_dx = r²/(Cr²+D), and only the first
        # has u/(Cr²+D) = λψ as its equilibrium.
        box = GradShafranovBox((10, 12), 2)
        s = box.space
        Wdx = Matrix(weighted_matrix(s,
            [p[1]^2 / herrnegger_mobility(p[1]) for p in quadrature_nodes(s)],
            (0, 0), (0, 0)))
        F = eigen(Symmetric(box.MΛ), Symmetric((Wdx .+ Wdx') ./ 2))
        v = F.vectors[:, argmax(real(F.values))]
        @test gs_fit(box, v)[3] > 0.1
        @test !isapprox(inv(maximum(real(F.values))), GS_LAMBDA_CONTINUUM; rtol = 0.1)
    end

    @testset "$(rpad("The residual is a CUBIC, so the Jacobian is analytic", 76))" begin
        # M = Cr²+D does not depend on the state, so 𝔾 is quadratic and ∂S/∂ĵ linear. Third
        # differences of a cubic in a fixed direction are constant.
        box = GradShafranovBox((6, 7), 2)
        f = gs_flow(box)
        v, d = randn(nbasis(box)), randn(nbasis(box))
        g(t) = vectorfield(f, v .+ t .* d)
        third(t, h) = (g(t + 2h) .- 3 .* g(t + h) .+ 3 .* g(t) .- g(t - h)) ./ h^3
        a, c = third(0.0, 0.25), third(0.9, 0.25)
        @test maximum(abs, a - c) / maximum(abs, a) < 1e-8
    end

    @testset "$(rpad("A short run CONSERVES H and dissipates S monotonically", 76))" begin
        spec = SECTION55_RUNS["c1"]
        box = GradShafranovBox((8, 9), 2)
        f = gs_flow(box)
        ĵ = gs_state(box, spec)
        λh = gs_eigenvalue(box)
        H₀, S = hamiltonian(f, ĵ), entropy(f, ĵ)
        integ = Integrator(f, ImplicitMidpoint(), spec.Δt; û₀ = copy(ĵ))
        for _ in 1:6
            integrate_step!(ĵ, integ)
            Snew = entropy(f, ĵ)
            @test Snew < S                       # exact for a quadratic entropy at any Δt
            @test Snew > λh * H₀                 # the Poincaré floor
            S = Snew
            @test abs(hamiltonian(f, ĵ) - H₀) / abs(H₀) < 1e-12
        end
    end

    @testset "$(rpad("C2's map is transcribed right, and DEGENERATES at s = 0", 76))" begin
        # The geometry, from three consequences of the printed constants — none of the four
        # numbers below appears in the paper.
        t = DiskTriangulation(32, 64)
        @test minimum(t.r)≈8.0 atol=1e-12
        @test maximum(t.r)≈16.0 atol=1e-12
        @test t.r[1]≈11.412924654785932 rtol=1e-12
        # The true area is 114.777, the Jacobian integrated over the parameter disk; an
        # inscribed polygon approaches it from BELOW at second order, so the tolerance names the
        # mesh. 114.605 at 32×64, 114.734 at 64×128, 114.766 at 128×256.
        @test disk_area(t)≈114.77699 rtol=2e-3
        @test disk_area(t) < 114.77699
        @test disk_area(DiskTriangulation(128, 256))≈114.77699 rtol=1e-4
        # The obstruction that defers C2's relaxation: s = 0 is one point, so a tensor-product
        # space on the parameter square is not even single-valued there.
        for θ in range(0, 2π; length = 9)
            @test all(disk_map(0.0, θ) .≈ disk_map(0.0, 0.0))
        end
        # Its reference eigenvalue is NOT deferred: P₁ on the physical triangulation reproduces
        # the printed 0.002599, from above, as a conforming Galerkin eigenvalue must.
        @test disk_eigenvalue(32, 64).λ≈GS_LAMBDA_DISK rtol=3e-3
        @test disk_eigenvalue(32, 64).λ > GS_LAMBDA_DISK_CONTINUUM
        @test GS_LAMBDA_DISK > GS_LAMBDA_DISK_CONTINUUM
    end
end
