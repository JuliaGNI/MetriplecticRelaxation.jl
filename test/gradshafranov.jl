using MetriplecticRelaxation
using MetriplecticRelaxation: integrate, energy, entropy
using MetriplecticRelaxation: SECTION55_RUNS, GradShafranovBox, gs_state, gs_flow, gs_fit,
                              gs_rayleigh, gs_current, gs_density, gs_stiffness,
                              gs_eigenvalue, separable_eigenvalue, TakedaGrid,
                              takeda_iterate, takeda_eigenvalue, herrnegger_mobility,
                              herrnegger_profile, GS_RADIAL, GS_AXIAL, HERRNEGGER_C,
                              HERRNEGGER_D, GS_LAMBDA_RECTANGLE, GS_LAMBDA_CONTINUUM,
                              disk_map, DiskTriangulation, disk_area, disk_eigenvalue,
                              GS_LAMBDA_DISK, GS_LAMBDA_DISK_CONTINUUM, gs_axes, gs_grid,
                              gs_ordinate_grid, gs_entropy_weight
using GeometricBrackets: nbasis, project, ∂y, vectorfield, issymmetric,
                         ispositive_semidefinite, degeneracy_residual, hamiltonian,
                         CollisionBracket, MetriplecticFlow, QuadraticHamiltonian,
                         Integrator, ImplicitMidpoint, integrate_step!, entropy_production,
                         stiffness_matrix, field, quadrature_weights, quadrature_nodes,
                         weighted_matrix, mass_matrix
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

    @testset "$(rpad("gs_grid SAMPLES r_i along the FIRST index, and sigma goes with it", 76))" begin
        (rs, zs) = gs_axes(21, 31)
        @test length(rs) == 21 && length(zs) == 31
        @test (rs[1], rs[end]) == GS_RADIAL
        @test (zs[1], zs[end]) == GS_AXIAL

        # Every function in V_D vanishes on ∂Ω, for RANDOM degrees of freedom and not merely for
        # the projection of something that vanished there already.
        box = GradShafranovBox((8, 9), 2)
        Z = gs_grid(box, randn(nbasis(box.space)), 21, 31)
        @test maximum(abs, vcat(Z[1, :], Z[end, :], Z[:, 1], Z[:, end])) < 1e-13

        # An asymmetric polynomial the degree-3 space contains exactly: reproduced to round-off,
        # and wrong by order one read the other way round. Degree 3 rather than the run's own
        # degree 2 because the boundary-vanishing bi-degree-(2,2) polynomials are a
        # one-dimensional span per axis, which would pass a transposed implementation. Equal node
        # counts, because at the 21x31 of the rows above a transposed read is a dimension
        # mismatch instead.
        f(r, z) = (r - 1) * (7 - r) * (z + 9.5)^2 * (9.5 - z)
        box3 = GradShafranovBox((5, 7), 3)
        W = gs_grid(box3, project(box3.space, x -> f(x[1], x[2])), 21, 21)
        (qs, ws) = gs_axes(21, 21)
        amp = maximum(abs(f(r, z)) for r in qs, z in ws)
        @test maximum(abs(W[i, j] - f(qs[i], ws[j])) for i in 1:21, j in 1:21) / amp < 1e-13
        @test maximum(abs(W[j, i] - f(qs[i], ws[j])) for i in 1:21, j in 1:21) / amp > 1e-3

        # And sigma(r) rides the RADIAL index. Getting this wrong is not a shape error — both
        # grids are 21x31 — so it needs its own assertion. Interior only: Z is zero on ∂Ω.
        ĵ = randn(nbasis(box.space))
        J = gs_grid(box, ĵ, 21, 31)
        Y = gs_ordinate_grid(box, ĵ, 21, 31)
        @test maximum(abs(Y[i, j] / J[i, j] - gs_entropy_weight(rs[i]))
        for i in 2:20, j in 2:30) < 1e-12
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
        # The obstruction C2's relaxation needs a `PolarSplineSpace` for: s = 0 is one point,
        # so a tensor-product space on the parameter square is not even single-valued there.
        for θ in range(0, 2π; length = 9)
            @test all(disk_map(0.0, θ) .≈ disk_map(0.0, 0.0))
        end
        # The reference eigenvalue needs no such space: P₁ on the physical triangulation
        # reproduces the printed 0.002599, as a conforming Galerkin eigenvalue must.
        @test disk_eigenvalue(32, 64).λ≈GS_LAMBDA_DISK rtol=3e-3
        @test disk_eigenvalue(32, 64).λ > GS_LAMBDA_DISK_CONTINUUM
        @test GS_LAMBDA_DISK > GS_LAMBDA_DISK_CONTINUUM
    end

    @testset "$(rpad("C2's flow runs on the map: H, S and the MASS Casimir", 76))" for state in (:dirichlet,
        :free)
        # The mesh is coarse on purpose — the run itself is 12×24, and what is under test here is
        # the mapped-domain path rather than the run's numbers.
        #
        # Both state spaces are run, because the mass Casimir is the one thing they disagree
        # about and a test of one of them says nothing about the other. Everything else — the
        # degeneracy, symmetry, semidefiniteness, energy, entropy, the Poincaré floor — must
        # hold in both, and does.
        spec = SECTION55_RUNS["c2"]
        disk = GradShafranovDisk((6, 12), 3; state = state)
        f = gs_flow(disk)
        ĵ = gs_state(disk, spec)
        λh = gs_eigenvalue(disk.space)

        # λ_h does not depend on the state space: the pencil (K^μ, B) carries no mass
        # constraint, so both spaces give the same lowest eigenvalue. That is exactly why an
        # eigenvalue check cannot see the difference this testset is about.
        @test λh ≈ gs_eigenvalue(GradShafranovDisk((6, 12), 3;
            state = state === :dirichlet ? :free : :dirichlet).space)
        # The degeneracy is asked against the flow's own ∂H/∂ĵ = 𝕄Λĵ, which carries the PHYSICAL
        # mass. A bracket built on the parameter frame misses this by twelve orders of magnitude;
        # `verify_gradshafranov_disk.jl` §4 measures both controls.
        @test degeneracy_residual(f, ĵ) < 1e-12
        @test degeneracy_residual(f, randn(nbasis(disk))) < 1e-12
        @test issymmetric(f.metric, ĵ)
        @test ispositive_semidefinite(f.metric, ĵ)
        @test entropy_production(f, ĵ) > 0

        H₀, S = hamiltonian(f, ĵ), entropy(f, ĵ)
        m₀ = integrate(disk, ĵ)
        integ = Integrator(f, ImplicitMidpoint(), spec.Δt; û₀ = copy(ĵ))
        for _ in 1:6
            integrate_step!(ĵ, integ)
            Snew = entropy(f, ĵ)
            @test Snew < S                       # exact for a quadratic entropy at any Δt
            @test Snew > λh * H₀                 # the Poincaré floor
            S = Snew
            @test abs(hamiltonian(f, ĵ) - H₀) / abs(H₀) < 1e-12
        end

        # The mass is the diagnostic that separates the two spaces, and the assertion inverts
        # with them. `:free` keeps the whole polar space, which contains the constant exactly by
        # the partition of unity the pole triangle preserves, so the bracket's mass Casimir
        # survives and ∫j dx is conserved. `:dirichlet` removes the constant with the rim row,
        # so it must **drift** — and that drift is the positive control that the constant has
        # actually left the space, without which "the run reaches eq:gs-ref" would be consistent
        # with the rim condition having done nothing.
        drift = abs(integrate(disk, ĵ) - m₀) / abs(m₀)
        if state === :free
            @test drift < 1e-12
        else
            @test drift > 1e-3
        end
    end
end
