#!/usr/bin/env julia
#
# Verify `src/takeda.jl`: the Grad-Shafranov problem of Section 5.5 and the classical iterative
# solver the manuscript's reference eigenvalues come from.
#
#     julia --project=scripts scripts/verify_takeda.jl
#
# THE CLAIM: both eigenvalues the manuscript prints -- lambda = 0.030302 for C1's rectangle and
# lambda = 0.002599 for C2's mapped disk -- are reproducible from its own description,
# INDEPENDENTLY of any relaxation run.
#
# They are reproducible, and neither is exact.  C1's continuum eigenvalue is 0.0302346260 --
# nine digits, three independent routes -- and the printed number sits 0.22 % above it; C2's is
# 0.0025970 and its printed number sits 0.077 % above.  A conforming Galerkin eigenvalue
# converges from above, so the sign is right in both cases and the size is what a low-order
# element on the stated grid gives.  Sections 4 and 7 measure that rather than asserting it.
#
# The controls are the point of the script.  Section 5.5's measure dmu = dr dz / r has to appear
# in the Delta-star operator AND in the profile matrix, and a stray factor of r in either is not
# a small error: on the rectangle it moves lambda by between four and five times its own value.
# So does mistaking the manuscript's printed CURRENT profile (Cr + D/r) for its state profile
# (Cr^2 + D).  Ten controls run in all and none of them can pass.
#
# One of them behaves differently on the two geometries, and that is reported as a finding
# rather than tuned away: dropping the measure from BOTH matrices largely CANCELS, because the
# weight sits on both sides of the Rayleigh quotient.  What survives is set by how far r varies
# -- a factor of 7 on the rectangle, where the error is 13 %, and a factor of 2 on the disk,
# where it is 2 %.  So C1 is the geometry that discriminates between the readings, and a check
# run only on C2 would be much weaker than it looks.

using MetriplecticRelaxation
using MetriplecticRelaxation: TensorSplineSpace, GS_RADIAL, GS_AXIAL,
                              GS_LAMBDA_RECTANGLE, GS_LAMBDA_CONTINUUM,
                              HERRNEGGER_C, HERRNEGGER_D,
                              gs_density, herrnegger_mobility, herrnegger_profile,
                              gs_stiffness, gs_eigenvalue, separable_eigenvalue,
                              TakedaGrid, takeda_iterate, takeda_eigenvalue, takeda_field,
                              takeda_order, takeda_extrapolate,
                              GS_LAMBDA_DISK, GS_LAMBDA_DISK_CONTINUUM, disk_map,
                              DiskTriangulation, disk_area, disk_eigenvalue
using PoissonBrackets: quadrature_nodes, quadrature_weights, weighted_matrix, project
using SimpleSplines: UniformMesh, Dirichlet
using LinearAlgebra
using SparseArrays
using Printf

include(joinpath(@__DIR__, "check.jl"))
using .Checks: header, check, summary

const C, D = HERRNEGGER_C, HERRNEGGER_D

"A tensor-product homogeneous-Dirichlet space on C1's rectangle."
function rectangle(n, p; nz = n)
    TensorSplineSpace((UniformMesh(n, GS_RADIAL), UniformMesh(nz, GS_AXIAL)), p, Dirichlet())
end

# The eigenvalue of `-Δ*ψ = λ f(r,ψ)` with the measure put in `lhs` places and `rhs` places
# independently, and with the profile `w`. `(true, true, herrnegger_mobility)` is the
# manuscript's reading; the other combinations are the controls of Section 5.
function variant_eigenvalue(s, lhs::Bool, rhs::Bool, w = herrnegger_mobility)
    x = quadrature_nodes(s)
    a = lhs ? gs_density.(x) : ones(length(x))
    b = [w(p[1]) * (rhs ? gs_density(p) : one(eltype(a))) for p in x]
    K = weighted_matrix(s, a, (1, 0), (1, 0)) .+ weighted_matrix(s, a, (0, 1), (0, 1))
    B = weighted_matrix(s, b, (0, 0), (0, 0))
    minimum(real, eigvals(Symmetric(Matrix(K)), Symmetric(Matrix(B))))
end

println("verify_takeda.jl  --  the Grad-Shafranov problem and its classical solver (§5.5)")
@printf("    Ω = [%.1f, %.1f] × [%.1f, %.1f]   C = %.1f   D = %.1f\n",
    GS_RADIAL..., GS_AXIAL..., C, D)

# =================================================================================================
header("1. the profile and the operator are the manuscript's")

# f(r,·) is DEFINED as the inverse of ∂_y s(r,·). Checking the inversion rather than the closed
# form is what makes the check independent of having done the algebra the same way twice.
∂ys(r, y) = y / (C * r^2 + D)

let rs = range(GS_RADIAL...; length = 7), ys = range(-3.0, 3.0; length = 9),
    e = maximum(abs(∂ys(r, herrnegger_profile(r, y)) - y) for r in rs, y in ys)

    check("f(r,·) inverts ∂_y s(r,·) — the manuscript's definition of f", e < 1e-15,
        @sprintf("max |∂_y s(r, f(r,y)) − y| = %.3e over 7×9 points", e))
end

# eq:M-condition, M ∂_y²s = 1, is what fixes the mobility. ∂_y²s = 1/(Cr²+D) exactly here, so
# this is an identity between two closed forms and holds to the last bit.
let rs = range(GS_RADIAL...; length = 13),
    e = maximum(abs(herrnegger_mobility(r) / (C * r^2 + D) - 1) for r in rs)

    check("eq:M-condition M ∂_y²s = 1 holds", e < 1e-15,
        @sprintf("max |M/(Cr²+D) − 1| = %.3e", e))
end

# The control on the profile: the manuscript ALSO prints (4π/c)J_φ = λ(Cr + D/r)ψ, which is the
# same equilibrium written for the current rather than for the state u = r(4π/c)J_φ. Reading
# that as f would be a plausible slip, and it does not invert ∂_y s.
let rs = range(GS_RADIAL...; length = 7), ys = range(0.5, 3.0; length = 5),
    e = maximum(abs(∂ys(r, (C * r + D / r) * y) - y) / abs(y) for r in rs, y in ys)

    check("CONTROL: the current profile (Cr + D/r) does NOT invert ∂_y s", e > 0.5,
        @sprintf("max relative defect %.4f — a factor of r, hence O(1) on [1,7]", e))
end

# `gs_stiffness` is -Δ* and not a weighted Laplacian: ∫∇ψ·∇v dμ = -∫ v Δ*ψ dμ, with Δ* applied
# analytically to a test function that vanishes on ∂Ω. This is the identity the whole μ-weighted
# formulation rests on, and it is what a stray factor of r in the operator breaks.
let s = rectangle(24, 3), x = quadrature_nodes(s), w = quadrature_weights(s)
    kr, kz = π / (GS_RADIAL[2] - GS_RADIAL[1]), π / (GS_AXIAL[2] - GS_AXIAL[1])
    ψf(p) = sin(kr * (p[1] - GS_RADIAL[1])) * sin(kz * (p[2] - GS_AXIAL[1]))
    # `vf` shares ψf's radial mode rather than being the next one up. A test function
    # orthogonal to ψf makes both sides of the identity vanish, which passes and says nothing;
    # the squared radial factor still vanishes at both ends and overlaps ψf.
    vf(p) = sin(kr * (p[1] - GS_RADIAL[1]))^2 * sin(kz * (p[2] - GS_AXIAL[1]))
    # Δ* = r ∂_r(r⁻¹∂_r) + ∂_z² = ∂_r² − r⁻¹∂_r + ∂_z², applied to ψf by hand.
    Δ✻ψ(p) =
        let (r, z) = p, sr = sin(kr * (r - GS_RADIAL[1])), cr = cos(kr * (r - GS_RADIAL[1]))
            (-kr^2 * sr - cr * kr / r - kz^2 * sr) * sin(kz * (z - GS_AXIAL[1]))
        end

    ψ̂, v̂ = project(s, ψf), project(s, vf)
    weak = dot(ψ̂, gs_stiffness(s), v̂)
    strong = -dot(w .* gs_density.(x), Δ✻ψ.(x) .* vf.(x))
    e = abs(weak - strong) / abs(strong)
    check("gs_stiffness IS −Δ*: ∫∇ψ·∇v dμ = −∫ v Δ*ψ dμ", e < 1e-6,
        @sprintf("weak %.12e   strong %.12e   relative %.3e", weak, strong, e))

    # The control: the same identity against the UNWEIGHTED stiffness matrix, i.e. against −Δ
    # instead of −Δ*. The two operators differ by the first-order term r⁻¹∂_r, which on [1,7] is
    # not a perturbation.
    plain = dot(ψ̂,
        weighted_matrix(s, ones(length(x)), (1, 0), (1, 0)) .+
        weighted_matrix(s, ones(length(x)), (0, 1), (0, 1)),
        v̂)
    check("CONTROL: the unweighted stiffness matrix is NOT −Δ*",
        abs(plain - strong) / abs(strong) > 0.5,
        @sprintf("plain %.6e against −∫vΔ*ψ dμ %.6e   relative %.3f",
            plain, strong, abs(plain - strong) / abs(strong)))
end

# =================================================================================================
header("2. Takeda's iteration converges to the smallest eigenvalue")

const g64 = TakedaGrid(64, 64)
const it64 = takeda_iterate(g64)

check("the iteration converges on the manuscript's 64×64 node grid", it64.converged,
    @sprintf("λ = %.12f after %d sweeps   last increment %.2e   %d interior unknowns",
        it64.λ, it64.iterations, it64.increment, length(g64)))

# λ comes out of a NORMALISATION (the axis value); the Rayleigh quotient is stationary at the
# fixed point and its error is second order in the eigenvector error where λ's is first order.
# The two agreeing is therefore a statement that the iteration converged, not a restatement of it.
let e = abs(it64.λ - it64.λ_rayleigh) / it64.λ
    check("the normalisation λ and the Rayleigh quotient agree", e < 1e-12,
        @sprintf("λ = %.14f   λ_R = %.14f   relative %.3e", it64.λ, it64.λ_rayleigh, e))
end

# An algorithmically independent solve of the same discrete problem: a dense symmetric
# generalised eigensolve, which shares nothing with the iteration but the two matrices.
for n in (16, 20, 24)
    g = TakedaGrid(n, n)
    λi, λe = takeda_iterate(g).λ, takeda_eigenvalue(g)
    check(
        @sprintf("at %d² nodes the fixed point IS the smallest generalised eigenvalue", n),
        abs(λi - λe) / λe < 1e-11,
        @sprintf("iteration %.14f   eigensolve %.14f   relative %.3e",
            λi, λe, abs(λi - λe) / λe))
end

# The relaxed state's axial structure is unambiguous only because the first two axial modes are
# far apart — a near-degeneracy would make "the" eigenfunction a two-parameter family.
let s1 = separable_eigenvalue(; mode = 1).λ, s2 = separable_eigenvalue(; mode = 2).λ

    check("the second axial mode is well separated from the first", s2 / s1 > 1.15,
        @sprintf("λ₁ = %.10f   λ₂ = %.10f   ratio %.4f", s1, s2, s2 / s1))
end

# =================================================================================================
header("3. the continuum eigenvalue, by three routes that share no assembly")

# Route 1: separation of variables plus a 1D spline Galerkin solve. Degree- and mesh-independent
# is the claim; the coarsest row is the control that the refinement is doing something.
let rows = [(n, p, separable_eigenvalue(; cells = n, degree = p).λ)
            for p in (2, 3, 4), n in (32, 64, 128)],
    fine = [λ for (n, p, λ) in rows if p ≥ 3], spread = maximum(fine) - minimum(fine)

    check(
        "the separable eigenvalue is degree- and mesh-independent at p ≥ 3", spread < 1e-9,
        @sprintf("spread over p ∈ {3,4} × n ∈ {32,64,128} = %.2e   λ = %.10f",
            spread, sum(fine) / length(fine)))
    check("it matches the recorded GS_LAMBDA_CONTINUUM",
        abs(sum(fine) / length(fine) - GS_LAMBDA_CONTINUUM) < 5e-10,
        @sprintf("measured %.12f   recorded %.10f   difference %.2e",
            sum(fine) / length(fine), GS_LAMBDA_CONTINUUM,
            sum(fine) / length(fine) - GS_LAMBDA_CONTINUUM))
end

# Route 2: a two-dimensional Galerkin solve that uses no separation at all. Agreement with
# route 1 is what makes the separation legitimate rather than assumed.
let λ2 = gs_eigenvalue(rectangle(24, 3)),
    e = abs(λ2 - GS_LAMBDA_CONTINUUM) / GS_LAMBDA_CONTINUUM

    check("the 2D Galerkin solve confirms the separation to 1e-8", e < 1e-8,
        @sprintf("2D λ = %.12f   separable %.10f   relative %.3e",
            λ2, GS_LAMBDA_CONTINUUM, e))
end

# Route 3: the second-order finite-volume iteration, refined. The ORDER is measured before the
# extrapolation is used, because Richardson at the wrong order manufactures digits.
const fd_ns = [16, 24, 32, 48, 64, 96]
const fd_λs = [takeda_iterate(TakedaGrid(n, n)).λ for n in fd_ns]

for (n, λ) in zip(fd_ns, fd_λs)
    println(@sprintf("      %3d² nodes   λ = %.12f   (λ − λ_∞)/λ_∞ = %+.3e",
        n, λ, (λ - GS_LAMBDA_CONTINUUM) / GS_LAMBDA_CONTINUUM))
end

let q = takeda_order(fd_λs, fd_ns)
    check("the finite-volume eigenvalue converges at second order", 1.8 < q < 2.2,
        @sprintf("fitted order %.4f over %d grids", q, length(fd_ns)))
end

let λ∞ = takeda_extrapolate(fd_λs[end - 1], fd_λs[end], fd_ns[end - 1], fd_ns[end]),
    e = abs(λ∞ - GS_LAMBDA_CONTINUUM) / GS_LAMBDA_CONTINUUM

    check("Richardson extrapolation lands on the separable value", e < 1e-5,
        @sprintf("extrapolated %.12f   separable %.10f   relative %.3e",
            λ∞, GS_LAMBDA_CONTINUUM, e))
end

# =================================================================================================
header("4. against the manuscript's printed λ = 0.030302")

# The claim under test is NOT that 0.030302 is the eigenvalue -- it is not, and Section 3 says so
# to nine digits. It is that 0.030302 is what a conforming discretisation of THIS problem on a
# 64 × 64 grid produces, and that no other reading of the problem comes anywhere near it.
let e = abs(it64.λ - GS_LAMBDA_RECTANGLE) / GS_LAMBDA_RECTANGLE
    check("the 64×64 finite-volume λ reproduces the printed λ to 0.3 %", e < 3e-3,
        @sprintf("takeda %.10f   printed %.6f   relative %+.3e   continuum %.10f",
            it64.λ, GS_LAMBDA_RECTANGLE,
            (it64.λ - GS_LAMBDA_RECTANGLE) /
            GS_LAMBDA_RECTANGLE, GS_LAMBDA_CONTINUUM))
end

# Where the last three digits come from: a conforming Galerkin eigenvalue approaches its limit
# FROM ABOVE, and the printed number is above the limit. Tensor-product Q₁ on this geometry
# brackets it -- coarser than 0.030302 at 16 cells, finer at 32 -- so 0.030302 is the number of a
# low-order element on a mesh of the stated order, and nothing sharper can be said without the
# authors' code.
let q1 = [(n, gs_eigenvalue(rectangle(n, 1))) for n in (16, 24, 32, 63)]
    for (n, λ) in q1
        println(@sprintf("      Q₁, %2d cells   λ = %.10f   (λ − λ_∞)/λ_∞ = %+.3e",
            n, λ, (λ - GS_LAMBDA_CONTINUUM) / GS_LAMBDA_CONTINUUM))
    end
    check("every Q₁ eigenvalue lies ABOVE the continuum limit, as a conforming one must",
        all(λ > GS_LAMBDA_CONTINUUM for (_, λ) in q1),
        @sprintf("smallest excess %+.3e at %d cells",
            minimum(λ - GS_LAMBDA_CONTINUUM for (_, λ) in q1), last(q1)[1]))
    check("the printed λ is bracketed by the Q₁ eigenvalues at 16 and 32 cells",
        q1[3][2] < GS_LAMBDA_RECTANGLE < q1[1][2],
        @sprintf("Q₁(32) = %.10f  <  %.6f  <  Q₁(16) = %.10f",
            q1[3][2], GS_LAMBDA_RECTANGLE, q1[1][2]))
end

# =================================================================================================
header("5. CONTROLS: the measure and the profile are load-bearing")

# `dμ = dr dz / r` has to appear in Δ* and in the profile matrix. Each of the three wrong
# combinations is O(1) wrong -- between 13 % and a factor of four -- so a stray factor of r
# cannot hide inside a tolerance. This is the check that the μ-weighting of `gs_stiffness` and
# `gs_profile_matrix` is not decoration.
let s = rectangle(20, 3), λ = variant_eigenvalue(s, true, true)
    check("the consistent pair reproduces the eigenvalue",
        abs(λ - GS_LAMBDA_CONTINUUM) / GS_LAMBDA_CONTINUUM < 1e-7,
        @sprintf("μ on both sides: λ = %.12f", λ))

    for (name, lhs, rhs) in (("no measure anywhere", false, false),
        ("μ in Δ* only", true, false),
        ("μ in the profile matrix only", false, true))
        λv = variant_eigenvalue(s, lhs, rhs)
        e = abs(λv - GS_LAMBDA_RECTANGLE) / GS_LAMBDA_RECTANGLE
        check("CONTROL: $(name) is O(1) wrong", e > 0.1,
            @sprintf("λ = %.8f against %.6f   relative %+.3f", λv, GS_LAMBDA_RECTANGLE, e))
    end

    λc = variant_eigenvalue(s, true, true, r -> C * r + D / r)
    let e = abs(λc - GS_LAMBDA_RECTANGLE) / GS_LAMBDA_RECTANGLE
        check("CONTROL: the current profile (Cr + D/r) is O(1) wrong", e > 0.1,
            @sprintf("λ = %.8f against %.6f   relative %+.3f",
                λc, GS_LAMBDA_RECTANGLE, e))
    end
end

# The iteration is written for a general profile, so the same control runs through it and not
# only through the eigensolver -- otherwise the controls would test `variant_eigenvalue` rather
# than `takeda_iterate`.
let g = TakedaGrid(32, 32), λ = takeda_iterate(g).λ,
    λc = takeda_iterate(g; profile = (r, y) -> (C * r + D / r) * y).λ

    check("CONTROL: the same substitution inside takeda_iterate is O(1) wrong",
        abs(λc - λ) / λ > 0.1,
        @sprintf("Herrnegger-Maschke %.10f   current profile %.10f   ratio %.3f",
            λ, λc, λc / λ))
end

# =================================================================================================
header("6. the eigenfunction is a single-signed interior mode")

# The manuscript's scatter plot is a line through the origin, which needs ψ single-signed: a
# node in ψ would make u/(Cr²+D) = λψ describe two branches. The lowest eigenfunction of a
# symmetric positive definite pair is single-signed, and this measures it rather than citing it.
let Ψ = takeda_field(g64, it64.ψ), interior = Ψ[2:(end - 1), 2:(end - 1)]
    check("ψ is single-signed on the interior", all(>(0), interior),
        @sprintf("min %.6e   max %.6f   at r = %.3f", minimum(interior), maximum(interior),
            g64.r[argmax(Ψ)[1]]))
    check("ψ vanishes on ∂Ω exactly",
        all(iszero, Ψ[1, :]) && all(iszero, Ψ[end, :]) && all(iszero, Ψ[:, 1]) &&
            all(iszero, Ψ[:, end]),
        "all four edges are identically zero by construction")
end

# =================================================================================================
header("7. C2's mapped domain and its reference eigenvalue λ = 0.002599")

# C2's RELAXATION is deferred -- the pole of eq:mapping needs a polar-spline space, and
# PoissonBrackets has no triangular DiscreteSpace either; `disk_eigenvalue` says both. Its
# reference EIGENVALUE is not deferred: a P₁ triangulation of the physical domain has the pole as
# an ordinary node, so the number the manuscript prints can be checked with nothing regularised.

# The map first. If it is mistranscribed the eigenvalue below is meaningless, so the geometry is
# measured against three independent consequences of the printed constants.
let t = DiskTriangulation(64, 128)
    check("the map's image is r ∈ [8, 16]",
        abs(minimum(t.r) - 8) < 1e-12 && abs(maximum(t.r) - 16) < 1e-12,
        @sprintf("r ∈ [%.12f, %.12f]", minimum(t.r), maximum(t.r)))
    check("its pole sits at r = 11.4129, just inboard of the Gaussian's r₀ = 12",
        abs(t.r[1] - 11.412924654785932) < 1e-9,
        @sprintf("pole at r = %.12f   z = %.3e   r₀ = 12", t.r[1], t.z[1]))
    # The true area is 114.77699, the Jacobian integrated over the parameter disk; the
    # triangulation approaches it from below at second order, hence the mesh-dependent bound.
    check("and its area is 114.777, comparable to C1's rectangle's 114",
        abs(disk_area(t) / 114.77699 - 1) < 5e-4 && disk_area(t) < 114.77699,
        @sprintf("|Ω_disk| = %.6f at 64×128 (exact 114.77699)   |Ω_rect| = %.1f   z ∈ [%.6f, %.6f]",
            disk_area(t), 6 * 19.0, minimum(t.z), maximum(t.z)))
end

# s = 0 collapses to a point: every node of the innermost ring of the PARAMETER disk maps to the
# same physical point. That is the obstruction, stated as a measurement rather than as prose.
let θs = range(0, 2π; length = 17)[1:(end - 1)], pts = [disk_map(0.0, θ) for θ in θs],
    spread = maximum(maximum(abs, p .- pts[1]) for p in pts)

    check("the map degenerates at s = 0 — the whole circle is one point", spread < 1e-14,
        @sprintf("spread of disk_map(0, θ) over 16 angles = %.3e — this is why C2's relaxation is deferred",
            spread))
end

const disk_ns = [(8, 16), (12, 24), (16, 32), (24, 48), (32, 64), (48, 96), (64, 128)]
const disk_λs = Float64[]

for (n, m) in disk_ns
    d = disk_eigenvalue(n, m)
    push!(disk_λs, d.λ)
    println(@sprintf("      %3d×%3d cells   dof = %5d   λ = %.10f   (λ − λ_∞)/λ_∞ = %+.3e",
        n, m, d.dof, d.λ, (d.λ - GS_LAMBDA_DISK_CONTINUUM) / GS_LAMBDA_DISK_CONTINUUM))
end

let q = takeda_order(disk_λs, [n for (n, _) in disk_ns])
    check(
        "the P₁ eigenvalue converges at second order on the mapped domain", 1.6 < q < 2.4,
        @sprintf("fitted order %.4f over %d refinements", q, length(disk_ns)))
end

let λ∞ = takeda_extrapolate(disk_λs[end - 1], disk_λs[end], 48, 64),
    e = abs(λ∞ - GS_LAMBDA_DISK_CONTINUUM) / GS_LAMBDA_DISK_CONTINUUM

    check("and extrapolates onto the recorded GS_LAMBDA_DISK_CONTINUUM", e < 1e-4,
        @sprintf("extrapolated %.10f   recorded %.7f   relative %+.3e",
            λ∞, GS_LAMBDA_DISK_CONTINUUM, e))
end

# The claim: 0.002599 is reproduced. As on the rectangle it is reproduced to its OWN
# discretisation error and not exactly, and it lies above the continuum limit because a
# conforming Galerkin eigenvalue does.
let e = abs(disk_λs[end] - GS_LAMBDA_DISK) / GS_LAMBDA_DISK
    check(
        "the 64×128 P₁ eigenvalue reproduces the printed λ = 0.002599 to 0.1 %", e < 1e-3,
        @sprintf("P₁ λ = %.10f   printed %.6f   relative %+.3e   continuum %.7f",
            disk_λs[end], GS_LAMBDA_DISK,
            (disk_λs[end] - GS_LAMBDA_DISK) / GS_LAMBDA_DISK,
            GS_LAMBDA_DISK_CONTINUUM))
    check("and the printed value lies above the continuum limit, from the right side",
        GS_LAMBDA_DISK > GS_LAMBDA_DISK_CONTINUUM,
        @sprintf("(printed − continuum)/continuum = %+.3e",
            (GS_LAMBDA_DISK - GS_LAMBDA_DISK_CONTINUUM) / GS_LAMBDA_DISK_CONTINUUM))
end

# The measure controls again, on the mapped domain, and they come out DIFFERENTLY here — which
# is a fact about the geometry and is reported rather than smoothed over. `disk_matrices` puts
# `dμ` in both matrices, so the pencil is rebuilt by hand with each half switchable.
function disk_variant(t, lhs::Bool, rhs::Bool)
    N = length(t.r)
    Is, Js, Ks, Bs = Int[], Int[], Float64[], Float64[]
    for e in t.triangles
        x = (t.r[e[1]], t.r[e[2]], t.r[e[3]])
        y = (t.z[e[1]], t.z[e[2]], t.z[e[3]])
        det = (x[2] - x[1]) * (y[3] - y[1]) - (x[3] - x[1]) * (y[2] - y[1])
        area = abs(det) / 2
        gx = (y[2] - y[3], y[3] - y[1], y[1] - y[2]) ./ det
        gy = (x[3] - x[2], x[1] - x[3], x[2] - x[1]) ./ det
        rc = (x[1] + x[2] + x[3]) / 3
        wk = lhs ? area / rc : area
        wb = rhs ? area * (C * rc^2 + D) / rc : area * (C * rc^2 + D)
        for p in 1:3, q in 1:3

            push!(Is, e[p])
            push!(Js, e[q])
            push!(Ks, wk * (gx[p] * gx[q] + gy[p] * gy[q]))
            push!(Bs, wb * (p == q ? 1 / 6 : 1 / 12))
        end
    end
    free = setdiff(1:N, t.boundary)
    K = Matrix(sparse(Is, Js, Ks, N, N)[free, free])
    B = Matrix(sparse(Is, Js, Bs, N, N)[free, free])
    minimum(real, eigvals(Symmetric(K), Symmetric(B)))
end

let t = DiskTriangulation(24, 48), λμ = disk_variant(t, true, true)
    for (name, lhs, rhs) in (("μ in Δ* only", true, false),
        ("μ in the profile matrix only", false, true))
        λv = disk_variant(t, lhs, rhs)
        e = abs(λv - λμ) / λμ
        check("CONTROL: on the mapped domain, $(name) is O(1) wrong", e > 0.1,
            @sprintf("λ = %.8f against %.8f   relative %+.3f", λv, λμ, e))
    end

    # FINDING, not a control that failed: dropping the measure from BOTH matrices is only 1.9 %
    # wrong here, against 13 % on the rectangle. The weight appears on both sides of the Rayleigh
    # quotient and largely cancels, and how much survives is set by how far r varies — a factor
    # of 2 on the disk (r ∈ [8,16]) against a factor of 7 on the rectangle (r ∈ [1,7]). So C1 is
    # the geometry that discriminates between the readings and C2 is not, which is worth knowing
    # before trusting a check run only on the disk.
    let λp = disk_variant(t, false, false), e = abs(λp - λμ) / λμ
        check(
            "FINDING: dropping dμ from BOTH matrices nearly cancels on the disk", e < 0.05,
            @sprintf("no measure: λ = %.8f   with it %.8f   relative %+.4f — against %+.3f on the rectangle, because r varies by 2 here and by 7 there",
                λp, λμ, (λp - λμ) / λμ, -0.128))
    end
end

summary("verify_takeda.jl")
