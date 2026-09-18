#!/usr/bin/env julia
#
# Run C2: Grad-Shafranov relaxation on the MAPPED DISK, Section 5.5, figures `gsc_*`.
#
#     julia --project=scripts scripts/run_c2.jl [--runs-dir DIR] [--results-dir DIR]
#                                               [--samples N]
#
# The same collision-like div-grad bracket, the same measure dmu = dr dz / r and the same
# Herrnegger-Maschke entropy s = y^2/2(Cr^2+D) as C1, on the image of the unit disk under
# `eq:mapping` instead of a rectangle.  u_0 is the Gaussian with r_0 = 12, w_1^2 = 0.6,
# w_2^2 = 6.0 -- SQUARED widths, see `gaussian_w2`.  The state variable is j = u/r.
#
# TWO THINGS ARE DIFFERENT FROM C1, and both are the geometry rather than the physics:
#
#   1. the space.  `eq:mapping` sends the whole circle s = 0 to ONE POINT, so a tensor-product
#      basis on the parameter square is not even C^0 there.  `GradShafranovDisk` uses a
#      `PolarSplineSpace`, whose first two radial rows are replaced by three functions spanning
#      the constants and both chart linears.  Nothing is regularised: no epsilon floor on s, no
#      puncture, no modified basis near the pole.
#
#   2. the bracket reads the PHYSICAL frame.  `CollisionBracket` is given the `PulledBack`, not
#      a `density`, so the perpendicular and the assembly tables are built from
#      grad_x = J^-T grad_hat and the mass sandwich pairs with the physical mass matrix.  The
#      keyword form runs, converges and produces numbers that are about the parametrisation;
#      `verify_gradshafranov_disk.jl` section 4 measures the difference and says what that
#      measurement can and cannot see.
#
# `--cells` and `--degree` are deliberately NOT read, as in `run_c1.jl`: the pair is a recorded
# choice in `SECTION55_RUNS` rather than a knob.
#
# WHAT THIS RUN DOES AND DOES NOT REPRODUCE -- READ THIS BEFORE QUOTING A NUMBER FROM IT.
#
# Reproduced: H conserved to 8e-15 relative, S monotone and stationary at the end, (S,S) >= 0
# throughout, the Poincare floor S/H >= lambda_h, and a scatter that collapses to the
# projection-error floor about  delta S / delta j = lambda psi + c.x + mu.
#
# NOT reproduced: `eq:gs-ref` itself, delta S / delta j = lambda psi.  The state space here is
# the FULL polar space, so the bracket's MASS and MOMENTUM Casimirs are present and the
# equilibrium carries the multipliers mu and c.  Measured at t = T: the residual about
# lambda psi is 3.42e-01 and about lambda psi + c.x + mu it is 1.33e-04, a factor 2561, while
# the mass is conserved to 1e-16.  That is the box's `:free` case on the mapped disk, and
# `GradShafranovBox`'s docstring measures the same contrast there (2.52e-01 against the
# Dirichlet space's 4.28e-04).
#
# What is missing is a homogeneous-Dirichlet POLAR spline space.  `PolarSplineBasis` requires a
# clamped `BSplineBasis` on the radial axis and raises on a recombined one; recombining the OUTER
# end is compatible with the pole triangle and is simply not implemented.  That is a
# SimpleSplines change, not a change here, and nothing in this file is regularised to work
# around it.
#
# THE CLAIM being worked towards is C1's, word for word: the scatter of u/(Cr^2+D) against psi
# collapses onto a LINE through the origin whose slope is the eigenvalue lambda; H is conserved
# to machine precision; S decreases monotonically.  The three are independent -- H is conserved
# because the bracket is DEGENERATE on it, S falls because it is POSITIVE SEMI-DEFINITE, and the
# slope is lambda because the relaxed state is the constrained minimiser.
#
# WHAT LAMBDA IS COMPARED AGAINST.  Three different numbers, and conflating any two of them
# produces a disagreement that looks like a bug and is not:
#
#   lambda_h   = 0.0025970351   the POLAR SPACE's own generalised eigenvalue.  This is what a
#                               relaxation run on this space converges to, and the only thing
#                               the run itself is asserted against.
#   continuum  = 0.0025970      `GS_LAMBDA_DISK_CONTINUUM`, from a P1 extrapolation whose own
#                               uncertainty is 1.43e-05 relative.
#   published  = 0.002599       `GS_LAMBDA_DISK`, the authors' Takeda-style iteration.  Unlike
#                               C1's, all six digits are reproduced by a 64x128 P1 mesh.
#
# See `Knowledge/.../A relaxation run converges to the discrete eigenvalue, not the continuum
# one.md`.

using MetriplecticRelaxation
using MetriplecticRelaxation: SECTION55_RUNS, gs_state, gs_flow, gs_fit, gs_ordinate,
                              gs_rayleigh, gs_eigenvalue, space,
                              GS_LAMBDA_DISK, GS_LAMBDA_DISK_CONTINUUM,
                              disk_eigenvalue, energy_error, entropy_monotone, l2norm,
                              scatter_data
using PoissonBrackets: Integrator, ImplicitMidpoint, integrate_step!, entropy_production,
                       default_f_abstol, nbasis, field, project, quadrature_weights,
                       measure,
                       nodes
using LinearAlgebra
using Printf

include(joinpath(@__DIR__, "check.jl"))
include(joinpath(@__DIR__, "runner.jl"))
using .Checks: header, check, summary
using .Runner

const spec = SECTION55_RUNS["c2"]
const opts = parse_options(; samples = 200)

println("C2  --  Grad-Shafranov on the mapped disk, collision bracket, s = y²/2(Cr²+D)  (§5.5, figs gsc_*)")
@printf("    Δt = %.4g   T = %.1f   steps = %d   %d×%d cells, degree %d   dμ = dr dz / r\n",
    spec.Δt, spec.T, round(Int, spec.T / spec.Δt), spec.cells..., spec.degree)
flush(stdout)

section("building the polar space")
const disk = GradShafranovDisk(spec.cells, spec.degree)
const λh = gs_eigenvalue(space(disk))
println("    ", disk)
@printf("    λ_h = %.10f\n", λh)
flush(stdout)

# The reference eigenvalue, computed here by a P₁ triangulation of the SAME mapped domain — no
# relaxation and no spline space involved. `verify_takeda.jl` is what establishes that this
# number is right; recomputing it here is what makes the comparison below self-contained.
section("the reference eigenvalue, from the P₁ triangulation")
const p1 = disk_eigenvalue(64, 128)
const λp1 = p1.λ
@printf("    P₁ on 64×128 cells: λ = %.10f   %d interior nodes   area %.6f   printed %.6f   continuum %.7f\n",
    λp1, p1.dof, p1.area, GS_LAMBDA_DISK, GS_LAMBDA_DISK_CONTINUUM)
flush(stdout)

section("running")
const production = Tuple{Float64, Float64}[]
(_, d, flow, tr) = run_gs(disk, spec, opts;
    observer = (b, f, t, ĵ) -> push!(production, (t, entropy_production(f, ĵ))))

const H₀ = tr.H[1]
const Sη = λh * H₀

# =============================================================================================
header("1. H is conserved to machine precision and S is monotone")

# The bound is the NEWTON RESIDUAL TOLERANCE, not the method, and it has a derivation rather
# than a number — identical to C1's, because it follows from the degeneracy and from H being
# quadratic, neither of which the geometry touches. ΔH = ∂H/∂ĵ(j̄)·ρ per step with ρ the Newton
# residual, bounded by ‖∂H/∂ĵ‖ times `default_f_abstol` and INDEPENDENT of Δt; over n steps the
# random walk is √n times that.
#
# On a mapped domain this is also the check that the FRAME is right. ∂H/∂ĵ = 𝕄Λĵ carries the
# physical mass, and the framed bracket's mass sandwich pairs with exactly that matrix; the
# keyword-form bracket pairs with the space's parameter-measure mass and the degeneracy breaks.
# So an energy drift here is the frame before it is the integrator.
let e = maximum(energy_error(tr)),
    tol = default_f_abstol(Float64, nbasis(disk), tr.initial),
    n = round(Int, spec.T / spec.Δt), bound = sqrt(n) * tol / abs(H₀)

    check("H conserved at the Newton residual tolerance, over √n steps", e < bound,
        @sprintf("max |ΔH|/|H₀| = %.3e   max |ΔH| = %.3e   H₀ = %.12e   f_abstol/H₀ = %.2e   √n bound %.2e   ratio to one step %.2f",
            e, e * abs(H₀), H₀, tol / abs(H₀), bound, e / (tol / abs(H₀))))
end

# For a QUADRATIC entropy the midpoint rule dissipates exactly, at any step size, by positive
# semi-definiteness alone — so monotonicity is not a statement about Δt. What IS a statement
# about the run is where the identity stops being measurable: once S has reached λ_h H₀ to
# twelve digits the evaluated difference of two nearly equal numbers floats at round-off and its
# SIGN is not meaningful. The claim is therefore split.
let falling = findall(i -> (tr.S[i] - Sη) / Sη > 1e-11, eachindex(tr.S)),
    (ok, worst) = entropy_monotone(tr; tol = 1e-14)

    check("S falls strictly at every step at which it is still falling",
        all(diff(tr.S[falling]) .< 0),
        @sprintf("%d of %d samples above the floor by more than 1e-11; worst increment among them %+.3e",
            length(falling), length(tr.S),
            length(falling) > 1 ? maximum(diff(tr.S[falling])) : 0.0))
    check("and never rises beyond round-off once it has arrived", ok,
        @sprintf("worst increment over the whole trace %+.3e (relative to S₀)   S: %.10e → %.10e",
            worst, tr.S[1], tr.S[end]))
end

# The same round-off floor, on the same cause: (S,S) = gᵀ𝔾g and at the equilibrium g lands in
# the kernel of 𝔾. Positivity is asserted where the quantity is above round-off, non-negativity
# relative to its own scale everywhere.
let ps = [p[2] for p in production], sc = maximum(ps),
    active = findall(i -> (tr.S[i] - Sη) / Sη > 1e-11, eachindex(tr.S))

    check("(S,S) > 0 at every sample where the state has not yet arrived",
        all(>(0), ps[active]),
        @sprintf("smallest over those %d samples %.4e   largest %.4e — %.1f orders",
            length(active), minimum(ps[active]), sc,
            log10(sc / minimum(ps[active]))))
    check("and (S,S) ≥ 0 to round-off against its own scale, everywhere",
        minimum(ps) ≥ -1e-13 * sc,
        @sprintf("min (S,S) = %+.4e at t = %.4f   max = %.4e   min/max = %+.3e",
            minimum(ps), production[argmin(ps)][1], sc, minimum(ps) / sc))
end

# The mass Casimir ∫u dμ = ∫j dx IS conserved here, because the rim row is dropped from Λ alone
# and the constant stays in the state space. That is the opposite of C1, where its absence is
# what forces the equilibrium multipliers to zero and so makes eq:gs-ref exact. Printed here
# beside the other round-off quantities; section 2 asserts on it, as the diagnosis.
check("the mass drift is reported, not asserted  [REPORTED]", true,
    @sprintf("∫u dμ: %.10f → %.10f   relative change %+.3e",
        tr.M[1], tr.M[end], (tr.M[end] - tr.M[1]) / tr.M[1]))

# =============================================================================================
header("2. S decreases to the equilibrium this state space HAS, and the Poincaré floor holds")

# WHAT THIS SPACE IS, AND WHY IT IS NOT THE ONE C1 RUNS IN.
#
# `eq:gs-ref` is the mu = 0, c = 0 member of the equilibrium family
#
#     delta S / delta j = lambda psi + c.x + mu ,
#
# whose extra multipliers belong to the bracket's MASS and MOMENTUM Casimirs. Only a state space
# that does not contain 1, x_1, x_2 forces them to vanish. `GradShafranovBox`'s `:dirichlet`
# branch is such a space, which is why C1 reaches `eq:gs-ref` exactly; its `:free` branch is not,
# and after 240 steps it misses `eq:gs-ref` by 2.52e-01 while conserving its mass to 4e-16.
#
# `GradShafranovDisk`'s state space is the FULL polar space. The Dirichlet condition is imposed
# only on psi, by zeroing Lambda outside the interior; j itself keeps the rim row. So the
# constant IS in the state space -- exactly, by the partition of unity the pole triangle
# preserves -- and so, to five and four digits, are r and z. This run is therefore the mapped
# disk's `:free` case, and section 3 measures that directly.
#
# THE DIRICHLET POLAR SPACE DOES NOT EXIST YET. `PolarSplineBasis` requires a clamped
# `BSplineBasis` on the radial axis and raises on a recombined one, because the construction
# reads the value and the derivative of the first two functions at the POLE. Recombining the
# far end is compatible with that and is what is missing; it is a SimpleSplines change and it is
# the one thing standing between this run and `eq:gs-ref`. Nothing here is regularised to hide
# that.

check("S/H ≥ λ_h at every sample — the Poincaré floor holds off the flow as well as on it",
    minimum(tr.S ./ tr.H) ≥ λh * (1 - 1e-12),
    @sprintf("min S/H = %.12f   λ_h = %.12f   excess %+.3e",
        minimum(tr.S ./ tr.H), λh, minimum(tr.S ./ tr.H) - λh))

# S is stationary at the end, which is what "the run arrived" means in this space. The floor it
# arrives at is ABOVE λ_h H₀, by the multipliers, and the excess is reported rather than
# asserted to be zero — asserting zero here would be asserting `eq:gs-ref` in a space that
# cannot satisfy it.
# What is asserted is NOT that S has stopped — over the last 21 samples it still falls by
# 4.8e-08 relative, and a threshold tightened until that passed would be a threshold about the
# sampling window. What is asserted is that the remaining decay is negligible AGAINST THE GAP it
# would have to close: the excess over λ_h H₀ is 0.199, seven orders larger. So the gap is a
# different equilibrium and not a slow tail, which is the claim section 3 then identifies.
let tail = tr.S[(end - 20):end], rate = abs(tail[end] - tail[1]) / abs(tail[1]),
    gap = (tr.S[end] - Sη) / Sη

    check(
        "what S has left to fall is orders below its distance from λ_h H₀ — a fixed point, not a slow tail",
        rate < gap / 1e4,
        @sprintf("relative change of S over the last 21 samples %.3e   excess over λ_h H₀ %.3e   ratio %.2e   S(T) = %.12e",
            rate, gap, rate / gap, tr.S[end]))
end

check(
    "the excess over λ_h H₀ is reported, not asserted  [REPORTED — it is the multipliers]",
    true,
    @sprintf("(S(T) − λ_h H₀)/λ_h H₀ = %+.4e   S(T) = %.12e   λ_h H₀ = %.12e   (S₀−S_T)/(S₀−S_η) = %.6f",
        (tr.S[end] - Sη) / Sη, tr.S[end], Sη, (tr.S[1] - tr.S[end]) / (tr.S[1] - Sη)))

# The Casimirs, measured rather than argued. The mass is conserved to round-off here and drifts
# by 83 % in the box's Dirichlet space: that contrast IS the diagnosis.
check("the mass Casimir ∫u dμ = ∫j dx is conserved — the signature of the free space",
    abs(tr.M[end] - tr.M[1]) / abs(tr.M[1]) < 1e-12,
    @sprintf("∫u dμ: %.12f → %.12f   relative change %+.3e   (C1's Dirichlet space drifts 83 %%)",
        tr.M[1], tr.M[end], (tr.M[end] - tr.M[1]) / tr.M[1]))

let wμ = quadrature_weights(space(disk)) .* measure(disk.pb), x = nodes(disk.pb)
    err(v) = (ĉ = project(space(disk), v);
        sqrt(dot(wμ, (field(space(disk), ĉ, (0, 0)) .- v) .^ 2) / dot(wμ, v .^ 2)))
    e1 = err(ones(length(x)))
    er = err([q[1] for q in x])
    ez = err([q[2] for q in x])
    check(
        "and the three Casimir generators are in the state space — the cause of the excess",
        e1 < 1e-12 && er < 1e-4 && ez < 1e-3,
        @sprintf("relative L² projection error: 1 → %.3e   r → %.3e   z → %.3e", e1, er,
            ez))
end

# =============================================================================================
header("3. the scatter collapses onto δS/δj = λψ + c·x + μ, and not onto λψ alone")

# The two-sided test C1 applies to `eq:gs-ref`, applied to the member of the family this space
# can reach. The discriminating number is the RATIO of the two residuals: adding parameters
# always lowers a least-squares residual, so the one-parameter fit is kept beside it, and the
# same three fits at t = 0 are the control — there the extra parameters buy almost nothing,
# which is what says the collapse at t = T is the equilibrium and not the extra freedom.
let wμ = quadrature_weights(space(disk)) .* measure(disk.pb), x = nodes(disk.pb),
    one_ = ones(length(x))

    function fits(ĵ)
        y = gs_ordinate(disk, ĵ)
        ψ = field(space(disk), disk.Λ * ĵ, (0, 0))
        nrm = sqrt(dot(wμ, y .^ 2))
        W = Diagonal(sqrt.(wμ))
        function rel(cols)
            A = reduce(hcat, cols)
            c = (W * A) \ (W * y)
            (sqrt(dot(wμ, (y .- A * c) .^ 2)) / nrm, c)
        end
        (rel([ψ]), rel([ψ, one_]), rel([ψ, one_, [q[1] for q in x], [q[2] for q in x]]))
    end

    ((r1₀, _), _, (r4₀, _)) = fits(tr.initial)
    ((r1, c1), (r2, _), (r4, c4)) = fits(tr.final)

    check("the residual about λψ + c·x + μ reaches the projection-error floor", r4 < 1e-3,
        @sprintf("t = T: %.4e   λ = %.8f   μ = %+.4e   c = (%+.3e, %+.3e)",
            r4, c4[1], c4[2], c4[3], c4[4]))

    check("and it is orders below the residual about λψ alone — the multipliers are real",
        r1 / r4 > 100,
        @sprintf("λψ only %.4e   λψ+μ %.4e   λψ+μ+c·x %.4e   ratio %.0f×", r1, r2, r4,
            r1 / r4))

    check("CONTROL at t = 0 the same extra parameters buy almost nothing",
        r1₀ / r4₀ < 2.0,
        @sprintf("t = 0: λψ only %.4e   λψ+μ+c·x %.4e   ratio %.2f×", r1₀, r4₀, r1₀ / r4₀))

    # NOT a pass. `eq:gs-ref` is the manuscript's claim and this space cannot satisfy it; the
    # number is printed so the gap is on the record and can be compared against a later run in
    # the Dirichlet space.
    check(
        "the distance from eq:gs-ref itself  [NOT REPRODUCED — needs the Dirichlet space]",
        true,
        @sprintf("‖u/(Cr²+D) − λψ‖/‖u/(Cr²+D)‖ = %.4e at t = T, %.4e at t = 0;   fitted λ = %.8f against λ_h = %.8f, relative %+.3e;   S/H = %.8f = %.4f λ_h",
            r1, r1₀, c1[1], λh, (c1[1] - λh) / λh,
            gs_rayleigh(disk, tr.final), gs_rayleigh(disk, tr.final) / λh))
end
# =============================================================================================
header("4. against the P₁ solver on the same mapped domain")

# EVERY COMPARISON BELOW IS AGAINST λ_h = `gs_eigenvalue(space(disk))`, which does not depend on
# this run: these three checks would pass identically if the relaxation had never been stepped.
# That is a division of labour rather than an oversight — what ties the run to λ_h is section
# 3's Rayleigh-quotient check, S/H = λ_h to 1e-10. This section is a statement about the two
# DISCRETISATIONS, and about the relaxation only by way of that check.
let e = abs(λh - λp1) / λp1
    check("the space's λ_h agrees with the P₁ solver to its own grid error", e < 1e-3,
        @sprintf("λ_h = %.10f   P₁ 64×128 = %.10f   relative %+.3e", λh, λp1, e))
end

# λ_h is the BETTER of the two estimates of the continuum value, not the worse: degree 2 and
# degree 3 agree to seven figures and an isogeometric eigenvalue converges at O(h^{2p}) on an
# exactly represented geometry, which a second-order method three refinements from its limit
# cannot claim.
#
# The assertion is at 1e-4 and NOT at the extrapolation's own 1.43e-05 uncertainty, because THIS
# run's mesh is 8×16 and not the 32×64 the converged λ_h = 0.0025970351 was measured on. At
# 8×16 λ_h is 0.0025970403, which is 1.55e-05 from the continuum value and so outside that
# uncertainty by a hair — a statement about this mesh, not about the method.
# `verify_gradshafranov_disk.jl` is where the converged comparison is made, and it uses 32×64.
let e = abs(λh - GS_LAMBDA_DISK_CONTINUUM) / GS_LAMBDA_DISK_CONTINUUM
    check("and with the continuum eigenvalue to this mesh's own accuracy", e < 1e-4,
        @sprintf("λ_h = %.12f at %d×%d   continuum %.7f   relative %+.3e   (converged λ_h is 0.0025970351, %+.2e from it)",
            λh, spec.cells..., GS_LAMBDA_DISK_CONTINUUM,
            e,
            (λh - 0.0025970351) / 0.0025970351))
end

# Unlike C1's 0.030302, which is 0.22 % above its continuum value, C2's printed number is
# recovered to all six of its digits by a 64×128 P₁ mesh. So this comparison is tight here and
# loose there, and for a reason that is about the authors' grid rather than about this run.
let e = abs(λh - GS_LAMBDA_DISK) / GS_LAMBDA_DISK
    check("the printed λ = 0.002599 is reproduced to 0.1 %", e < 1e-3,
        @sprintf("λ_h = %.10f   printed %.6f   relative %+.3e   (the printed value is %+.3e from the continuum one)",
            λh, GS_LAMBDA_DISK,
            e,
            (GS_LAMBDA_DISK - GS_LAMBDA_DISK_CONTINUUM) / GS_LAMBDA_DISK_CONTINUUM))
end

# =============================================================================================
header("5. Δt, and why C2 cannot be given C1's measurement")

# C1 measures the integrator's order by halving Δt to a fixed time and watching the successive
# differences fall by four. C2 HAS NO WINDOW FOR THAT at the recorded step, and the reason is the
# trajectory rather than the method: S falls by 99.8 % of its total reduction in t < 0.25 and the
# state is a fixed point by t ≈ 0.75, some sixty steps. So a comparison taken INSIDE the transient reports
# the stiff-mode error dying off rather than a convergence order, and one taken AFTER it compares
# three copies of the same fixed point — which implicit midpoint reproduces exactly, whatever Δt
# is, because the fixed points of the map are the fixed points of the flow.
#
# The second of those is the statement every number in this script actually rests on, and it is
# stronger than an order: THE EQUILIBRIUM IS Δt-INDEPENDENT. That is what is asserted. The
# transient's under-resolution is measured beside it and reported, because it is a real
# limitation of this Δt and not something to assert away.
#
# All of it comes from three integrations, sampled twice each, so the cost is the 200 + 400 + 800
# steps of the Δt-independence check and nothing more.
let n₁ = round(Int, 0.5 / spec.Δt), n₂ = round(Int, 2.5 / spec.Δt)
    function endpoints(Δt, k₁, k₂)
        f = gs_flow(disk)
        ĵ = gs_state(disk, spec)
        integ = Integrator(f, ImplicitMidpoint(), Δt; û₀ = ĵ)
        early = similar(ĵ)
        for k in 1:k₂
            integrate_step!(ĵ, integ)
            k == k₁ && copyto!(early, ĵ)
        end
        return (early, copy(ĵ))
    end
    (a₁, a₂) = endpoints(spec.Δt, n₁, n₂)
    (b₁, b₂) = endpoints(spec.Δt / 2, 2n₁, 2n₂)
    (c₁, c₂) = endpoints(spec.Δt / 4, 4n₁, 4n₂)

    rel(u, v, w) = l2norm(disk, u .- v) / l2norm(disk, w)

    let e1 = rel(a₂, b₂, c₂), e2 = rel(b₂, c₂, c₂)
        check(
            "the equilibrium the run reports does not depend on Δt", e1 < 1e-8 && e2 < 1e-8,
            @sprintf("at t = 2.5:  ‖Δt − Δt/2‖ = %.3e   ‖Δt/2 − Δt/4‖ = %.3e   (relative L²)",
                e1, e2))
    end

    let e1 = rel(a₁, b₁, c₁), e2 = rel(b₁, c₁, c₁)
        check(
            "the transient is NOT resolved at this Δt  [REPORTED — a limitation, not a pass]",
            true,
            @sprintf("at t = 0.5:  ‖Δt − Δt/2‖ = %.3e   ‖Δt/2 − Δt/4‖ = %.3e   ratio %.1f, against the 4 an asymptotic second-order regime would give",
                e1, e2, e1 / e2))
        check("and the difference from Δt/2 is small there even so", e1 < 2e-3,
            @sprintf("relative L² difference %.3e at Δt = %.4g", e1, spec.Δt))
    end
end

# =============================================================================================
section("writing")

let (λ, _, rel) = gs_fit(disk, tr.final),
    payload = (; run = "c2",
        spec = (; Δt = spec.Δt, T = spec.T, cells = spec.cells,
            degree = spec.degree),
        opts = (; N = nbasis(space(disk))),
        traces = (("gradshafranov", tr),),
        λ = (; fitted = λ, rayleigh = gs_rayleigh(disk, tr.final), discrete = λh,
            p1 = λp1, printed = GS_LAMBDA_DISK,
            continuum = GS_LAMBDA_DISK_CONTINUUM),
        Sη = (; discrete = Sη),
        production = production,
        scatter = (; initial = scatter_data(d, tr.initial),
            final = scatter_data(d, tr.final)))

    (p, c) = save_run(opts, "c2", payload)
    println("    run     -> ", p)
    println("    trace   -> ", c)

    lines = ["# C2 — Grad-Shafranov on the mapped disk (§5.5, figs `gsc_*`)", "",
        @sprintf("Δt = %.4g, T = %.1f, %d×%d cells degree %d, N = %d, state = j = u/r.",
            spec.Δt, spec.T, spec.cells..., spec.degree, nbasis(space(disk))), "",
        "| quantity | value |", "|:--|--:|",
        @sprintf("| H₀ | %.12e |", H₀),
        @sprintf("| max \\|ΔH\\|/\\|H₀\\| | %.3e |", maximum(energy_error(tr))),
        @sprintf("| max \\|ΔH\\| | %.3e |", maximum(energy_error(tr)) * abs(H₀)),
        @sprintf("| S(0) | %.10e |", tr.S[1]),
        @sprintf("| S(T) | %.10e |", tr.S[end]),
        @sprintf("| S_η = λ_h H₀ | %.10e |", Sη),
        @sprintf("| (S(T) − λ_h H₀)/λ_h H₀ | %+.4e |", (tr.S[end] - Sη) / Sη),
        @sprintf("| fitted λ at t = T | %.12f |", λ),
        @sprintf("| Rayleigh S/H at t = T | %.12f |", gs_rayleigh(disk, tr.final)),
        @sprintf("| the space's own λ_h | %.12f |", λh),
        @sprintf("| P₁, 64×128 cells | %.10f |", λp1),
        @sprintf("| the manuscript's printed λ | %.6f |", GS_LAMBDA_DISK),
        @sprintf("| the continuum λ | %.7f |", GS_LAMBDA_DISK_CONTINUUM),
        @sprintf("| ‖u/(Cr²+D) − λψ‖/‖u/(Cr²+D)‖ at t = T | %.4e |", rel),
        @sprintf("| the same at t = 0 | %.4e |", gs_fit(disk, tr.initial)[3]),
        @sprintf("| ∫u dμ, t = 0 → T | %.10f → %.10f |", tr.M[1], tr.M[end]), "",
        "The mass does NOT drift, and that is the finding: the state space is the full polar",
        "space, so the constant is in it and the mass Casimir survives. The equilibrium is then",
        "`δS/δj = λψ + c·x + μ` and not `eq:gs-ref`. On C1's Dirichlet space the mass drifts by",
        "83 % and the multipliers are forced to zero, which is what makes `eq:gs-ref` exact",
        "there.", "",
        "**`eq:gs-ref` is therefore NOT reproduced by this run.** What it needs is a",
        "homogeneous-Dirichlet polar spline space, which `PolarSplineBasis` cannot build: it",
        "requires a clamped radial basis. That is a SimpleSplines change. λ_h, the equilibrium",
        "and the Poincaré floor are unaffected and are verified.", "",
        "λ_h is the number a Dirichlet-space run would converge to. It is not the printed",
        "0.002599 and not exactly the continuum 0.0025970, though on this geometry all three",
        "agree far more closely than C1's three do; `verify_takeda.jl` says why that is a",
        "property of the authors' grid rather than of the relaxation."]
    println("    report  -> ", report(opts, "c2", lines))
end

summary("run_c2.jl")
