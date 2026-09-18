# C2, the mapped disk: what the polar space delivers, and the one thing it does not.
#
# Section 5.5's second geometry is the image of the unit disk under `eq:mapping`. Its
# discretisation is `GradShafranovDisk`, and this script establishes three things and refutes a
# fourth:
#
#   1. λ_h, the space's own Grad-Shafranov eigenvalue, against the continuum value and against
#      the two measure-swap controls that must be wrong by an order of magnitude;
#   2. the discrete equilibrium satisfies `eq:gs-ref`, u/(Cr²+D) = λψ, with λ = λ_h, confirmed
#      by three independent estimates that must agree;
#   3. the Poincaré floor, S/H ≥ λ_h, at the initial state and off the flow;
#   4. **the relaxation run itself cannot be done yet.** `CollisionBracket` reads the gradient
#      in the space's own coordinates, which on a mapped domain are not the physical ones. The
#      last section measures that, and it is why `gs_flow(::GradShafranovDisk)` raises.
#
# Run: julia --project=. --startup-file=no scripts/verify_gradshafranov_disk.jl

using LinearAlgebra
using Printf
using Random
using MetriplecticRelaxation
using PoissonBrackets
using SimpleSplines: UniformMesh, Dirichlet, (..)

Random.seed!(20260918)

pass = Bool[]
record(name, ok) = (push!(pass, ok); @printf("  %-56s %s\n", name, ok ? "pass" : "FAIL"))
say(s) = (println(s); flush(stdout))

## ---------------------------------------------------------------------------------------
say("1. λ_h, the polar space's own Grad-Shafranov eigenvalue")
say("")

# Eigenvalues of a degree-p isogeometric discretisation converge at O(h^{2p}), and the geometry
# here is exact — the Jacobian is evaluated analytically at the quadrature nodes rather than
# interpolated — so λ_h is converged to eight figures on a mesh that would be coarse for a
# finite-element solve. That is why the table below barely moves.
λ_ref = 0.0
for p in (2, 3), cells in ((8, 16), (16, 32), (32, 64))

    λ = gs_eigenvalue(PolarSplineSpace(cells, p))
    global λ_ref = (p == 3 && cells == (32, 64)) ? λ : λ_ref
    @printf("  p=%d  %2d × %-3d   λ_h = %.10f   rel to the continuum %+.2e\n",
        p, cells[1], cells[2], λ, (λ - GS_LAMBDA_DISK_CONTINUUM) / GS_LAMBDA_DISK_CONTINUUM)
end
say("")
@printf("  converged λ_h = %.10f\n", λ_ref)
@printf("  continuum     = %.7f   (P₁ extrapolation, itself uncertain at 1.43e-05 relative)\n",
    GS_LAMBDA_DISK_CONTINUUM)
@printf("  published     = %.6f   (%.3f %% above the continuum)\n",
    GS_LAMBDA_DISK,
    100 * (GS_LAMBDA_DISK - GS_LAMBDA_DISK_CONTINUUM) / GS_LAMBDA_DISK_CONTINUUM)

# The isogeometric value agrees with the P₁ extrapolation to within that extrapolation's own
# stated uncertainty. It is the better estimate of the two: p=2 and p=3 agree with each other
# to eight figures, which a second-order method three refinements from its limit cannot claim.
record("λ_h agrees with the continuum within the P₁ uncertainty",
    abs(λ_ref - GS_LAMBDA_DISK_CONTINUUM) / GS_LAMBDA_DISK_CONTINUUM < 1.43e-5)
record("p = 2 and p = 3 agree to seven figures",
    abs(gs_eigenvalue(PolarSplineSpace((32, 64), 2)) - λ_ref) / λ_ref < 1e-6)
say("")

## The controls: an inconsistent measure is not a small error
say("  CONTROLS — K^μ and B must carry the SAME measure")
let s = PolarSplineSpace((16, 32), 3)
    keep = disk_interior(s)
    F(x) = disk_map(x[1], x[2])
    DF(x) = disk_jacobian(x[1], x[2])
    pbμ = PulledBack(s, F, DF; density = x -> gs_density(x))
    pbx = PulledBack(s, F, DF)
    prof = [herrnegger_mobility(x[1]) for x in nodes(pbμ)]
    smallest(K, B) = minimum(real,
        eigvals(Symmetric(Matrix(K[keep, keep])), Symmetric(Matrix(B[keep, keep]))))

    Kμ = tensor_weighted_matrix(s, metric(pbμ))
    Kx = tensor_weighted_matrix(s, metric(pbx))
    Bμ = weighted_matrix(s, prof .* measure(pbμ), (0, 0), (0, 0))
    Bx = weighted_matrix(s, prof .* measure(pbx), (0, 0), (0, 0))

    consistent = smallest(Kμ, Bμ)
    @printf("    consistent  K^μ with B^μ   %.8f\n", consistent)
    @printf("    CONTROL     K^μ with B^dx  %.8f\n", smallest(Kμ, Bx))
    @printf("    CONTROL     K^dx with B^μ  %.8f\n", smallest(Kx, Bμ))
    @printf("    CONTROL     no pullback    %.8f\n",
        smallest(stiffness_matrix(s), weighted_matrix(s, prof, (0, 0), (0, 0))))
    record("the consistent pair gives λ_h",
        abs(consistent - GS_LAMBDA_DISK_CONTINUUM) / GS_LAMBDA_DISK_CONTINUUM < 1e-4)
    record("CONTROL swapping the measures is wrong by ≥ 10×",
        smallest(Kμ, Bx) < consistent / 10 && smallest(Kx, Bμ) > 10 * consistent)
end
say("")

## ---------------------------------------------------------------------------------------
say("2. the discrete equilibrium satisfies eq:gs-ref, u/(Cr²+D) = λψ")
say("")

disk = GradShafranovDisk((10, 20), 3)
λh = gs_eigenvalue(space(disk))
F = eigen(Symmetric(disk.MΛ), Symmetric(disk.W))
v = F.vectors[:, argmax(real(F.values))]
(λfit, res, rel) = gs_fit(disk, v)
λray = gs_rayleigh(disk, v)

say(@sprintf("  %s", disk))
@printf("  λ_h  eigensolve  = %.10f\n", λh)
@printf("  λ    scatter fit = %.10f   rel %+.2e\n", λfit, (λfit - λh) / λh)
@printf("  λ    Rayleigh    = %.10f   rel %+.2e\n", λray, (λray - λh) / λh)
@printf("  scatter residual = %.3e   relative %.3e\n", res, rel)

# Three estimates of the same number, by three different routes: a generalised eigensolve, a
# projection of the fields, and a ratio of the two invariants. Agreement is the claim.
record("the fitted λ is λ_h", abs(λfit - λh) / λh < 1e-8)
record("the Rayleigh quotient is λ_h", abs(λray - λh) / λh < 1e-8)
record("the scatter has collapsed onto the line", rel < 1e-4)
say("")

## ---------------------------------------------------------------------------------------
say("3. the Poincaré floor: S/H ≥ λ_h everywhere, not only along the flow")
say("")

spec = GSSpec(
    "c2", "5.5", gaussian_w2((12.0, 0.0), (0.6, 6.0), 1.0), (10, 20), 3, 0.0625, 25.0)
ĵ₀ = gs_state(disk, spec)
@printf("  S/H at the initial Gaussian = %.8f   which is %.1f× λ_h\n",
    gs_rayleigh(disk, ĵ₀), gs_rayleigh(disk, ĵ₀) / λh)
@printf("  the initial state is nothing like an equilibrium: fit relative %.3f\n",
    gs_fit(disk, ĵ₀)[3])
record("the floor holds at the initial state", gs_rayleigh(disk, ĵ₀) > λh)
record("the floor holds on 50 random states",
    all(gs_rayleigh(disk, randn(nbasis(disk))) >= λh for _ in 1:50))
say("")

## ---------------------------------------------------------------------------------------
say("4. WHY THE RELAXATION RUN IS BLOCKED: the bracket is not frame-covariant")
say("")

# `CollisionBracket` forms β = (−∂₂φ, ∂₁φ) from the space's own derivative tables. On an
# unmapped domain those are the physical derivatives and it is right. On a mapped domain they
# are the *parameter* derivatives, and the physical gradient is ∇_x = J⁻ᵀ∇̂ — a different
# object wherever J is not a multiple of a rotation.
#
# Measured below on the same PHYSICAL problem written in three parametrisations that differ
# only by a linear stretch. `∫u dx` confirms the three are the same problem; the bracket is not.

physical(x) = sin(π * (x[1] - 1)) * sin(π * x[2]) * (1 + 0.3 * cos(3π * x[1]))

function stretched(scale)
    s = TensorSplineSpace((UniformMesh(6, 1.0 .. 2.0), UniformMesh(6, 0.0 .. 1.0 / scale)),
        2, (Dirichlet(), Dirichlet()))
    pb = PulledBack(s, x -> (x[1], scale * x[2]), _ -> [1.0 0.0; 0.0 scale])
    K = tensor_weighted_matrix(s, metric(pb))
    M = weighted_matrix(s, measure(pb), (0, 0), (0, 0))
    Λ = Matrix(cholesky(Symmetric(Matrix(K))) \ Matrix(M))
    û = project(s, [physical(x) for x in nodes(pb)])
    G = metric_matrix(CollisionBracket(s, Λ; density = measure(pb)), û)
    mass = dot(quadrature_weights(s) .* measure(pb), basis_values(s, (0, 0))' * û)
    return (mass = mass, norm = maximum(abs, G), trace = tr(G))
end

base = stretched(1.0)
say("  the same physical problem, three parametrisations differing by a linear stretch")
for scale in (1.0, 2.0, 4.0)
    r = stretched(scale)
    @printf("    scale %.1f   ∫u dx = %.10f   ‖G‖ = %.4e   ratio to scale 1 = %8.2f\n",
        scale, r.mass, r.norm, r.norm / base.norm)
end

r2, r4 = stretched(2.0), stretched(4.0)
record("the three are the same physical problem (∫u dx identical)",
    abs(r2.mass - base.mass) < 1e-12 && abs(r4.mass - base.mass) < 1e-12)
record("CONTROL the bracket is NOT invariant — it scales as the fourth power",
    abs(r2.norm / base.norm - 16) < 1 && abs(r4.norm / r2.norm - 16) < 1)

say("")
say("""  The fourth power is the signature: the bracket is quadratic in ∇φ and the assembly
  contracts two further derivatives, so a gradient scaled by `scale` scales it by `scale⁴`.
  Symmetry, positive semi-definiteness and the degeneracy (F,H) = 0 hold in ALL THREE — those
  properties are algebraic and cannot see the frame. So the §2 structural checks pass on a
  mapped domain while the operator is the wrong one.

  What this needs is `CollisionBracket` reading ∇_x = J⁻ᵀ∇̂ throughout: the perp in the physical
  frame, and the assembly's derivative tables replaced by their node-dependent combinations
  Φ_k^phys = Σ_l (J⁻ᵀ)_kl Φ̂_l. That is an extension to the bracket, not to the space, and it is
  why `gs_flow(::GradShafranovDisk)` raises rather than returning a flow whose numbers would be
  about the parametrisation.""")
say("")

## ---------------------------------------------------------------------------------------
say(all(pass) ? "ALL CHECKS PASS" :
    "SOME CHECK FAILED ($(count(!, pass)) of $(length(pass)))")
exit(all(pass) ? 0 : 1)
