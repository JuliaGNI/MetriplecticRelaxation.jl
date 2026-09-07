#!/usr/bin/env julia
#
# The uniform resampling the Section 5.4 field maps are drawn on.
#
#     julia --project=scripts scripts/verify_euler_grid.jl
#
# `euler_grid` evaluates a Section 5.4 state off the space's own quadrature grid onto a uniform
# one.  It is four lines long and it is exactly the kind of code that looks right and is silently
# wrong -- a transposed axis, an off-by-one on the endpoints, a grid over the wrong interval --
# so every claim it makes is settled here in closed form rather than by looking at the picture
# it produces.
#
# Four sections, in the order a defect in one would invalidate the next:
#
#   1. the axis: N nodes, BOTH endpoints, uniform spacing, and the domain it covers -- checked
#      independently by integrating a resampled field with the trapezoidal rule and comparing
#      against the space's own quadrature;
#   2. exactness: a polynomial the space contains is reproduced to ROUND-OFF, at both degrees,
#      with a polynomial it does not contain as the control that says the claim is about the
#      degree rather than something the resampler returns for anything;
#   3. the homogeneous-Dirichlet boundary: all four edges vanish to round-off, for RANDOM
#      coefficients rather than for a projected function that happens to vanish there -- with the
#      `:free` space as the control that MUST fail;
#   4. the index convention.  The test function carries B1's own widths -- w_1^2 = 0.01 against
#      w_2^2 = 0.07, a factor of 2.6 and nothing else -- because a radially symmetric bump passes
#      a transposed implementation exactly.  The check is TWO-SIDED: the direct comparison
#      converges at the projection order under refinement, and the transposed comparison must be
#      wrong by order one.  B1's Gaussian itself does not vanish on the boundary, which is a
#      property of Section 5.4 rather than of the resampling, so the section reports that
#      separately and multiplies it by sin(pi x_1) sin(pi x_2) for the convergence rate.
#
# Everything here runs on spaces of at most 4096 degrees of freedom and the whole script is
# seconds.

using MetriplecticRelaxation
using MetriplecticRelaxation: SECTION5_RUNS, EulerSquare, euler_axis, euler_grid,
                              euler_state, initial_condition, integrate, state_extrema,
                              SQUARE_LENGTH
using PoissonBrackets: nbasis, project
using Printf
using Random

include(joinpath(@__DIR__, "check.jl"))
using .Checks: header, check, check_refined, summary

Random.seed!(0x9d31c07f)

println("verify_euler_grid.jl  --  the uniform resampling behind the §5.4 field maps")

# The number of samples per axis the figure script uses. Odd on purpose: B1, B2 and B3 share an
# initial condition centred at (1/2,1/2), and an odd node count puts a sample exactly on the
# peak rather than straddling it.
const SAMPLES = 129

# =============================================================================================
header("1. the axis is N nodes over [0,1], endpoints included")

let xs = euler_axis(SAMPLES)
    check("the axis has the requested number of nodes", length(xs) == SAMPLES,
        @sprintf("length = %d   requested %d", length(xs), SAMPLES))
    check("it starts at 0 and ends at the side length — BOTH endpoints",
        xs[1] == 0.0 && xs[end] == SQUARE_LENGTH,
        @sprintf("x₁ = %.17g   x_N = %.17g   L = %.17g", xs[1], xs[end], SQUARE_LENGTH))
    let d = diff(xs), h = SQUARE_LENGTH / (SAMPLES - 1)
        check("the spacing is uniform at L/(N-1)",
            maximum(abs, d .- h) < 1e-15,
            @sprintf("max |Δx - h| = %.3e   h = %.17g", maximum(abs, d .- h), h))
    end
    # The centre is a node, which is where the Gaussian's peak is. Not decoration: the field map
    # of an initial condition whose peak falls between two samples under-reports its own maximum,
    # and the colour range of `figure_fields` is shared between the two panels.
    check("the centre of the square is a node, for an odd N",
        any(x -> x == SQUARE_LENGTH / 2, xs),
        @sprintf("N = %d is odd: %s", SAMPLES, string(isodd(SAMPLES))))
end

# The independent statement about which domain the grid covers. If the axis ran over [0,2] or
# omitted an endpoint, the trapezoidal integral of a resampled field would not agree with the
# space's own quadrature — and that comparison shares no code with `euler_axis`.
let sq = EulerSquare(26, 2), spec = SECTION5_RUNS["b1"], ω̂ = euler_state(sq, spec)
    Z = euler_grid(sq, ω̂, SAMPLES)
    h = SQUARE_LENGTH / (SAMPLES - 1)
    w = fill(1.0, SAMPLES)
    w[1] = w[end] = 0.5
    trap = h^2 * sum(w[i] * w[j] * Z[i, j] for i in 1:SAMPLES, j in 1:SAMPLES)
    exact = integrate(sq, ω̂)
    e = abs(trap - exact) / abs(exact)
    check("the trapezoidal integral of the resampled field IS the space's own ∫ω",
        e < 2e-4,
        @sprintf("trapezoid %.10e   ∫ω = %.10e   relative %.3e", trap, exact, e))
end

# =============================================================================================
header("2. a polynomial the space contains is reproduced to round-off")

# The homogeneous-Dirichlet space of degree p contains, per axis, exactly the polynomials of
# degree ≤ p that vanish at both ends. At p = 2 that is span{x(1-x)}, which is one-dimensional
# and therefore SYMMETRIC under a transposition of the axes — so it settles exactness and says
# nothing about the index convention, which is section 4's business. At p = 3 it is
# span{x(1-x), x²(1-x)}, and the product below is asymmetric.
const p2_poly = (x, y) -> x * (1 - x) * y * (1 - y)
const p3_poly = (x, y) -> x * (1 - x) * y^2 * (1 - y)

function grid_error(sq, f, N)
    û = project(sq.space, x -> f(x[1], x[2]))
    Z = euler_grid(sq, û, N)
    xs = euler_axis(N)
    amp = maximum(abs(f(x, y)) for x in xs, y in xs)
    return (maximum(abs(Z[i, j] - f(xs[i], xs[j])) for i in 1:N, j in 1:N), amp)
end

let (e, amp) = grid_error(EulerSquare(6, 2), p2_poly, 41)
    check("degree 2 reproduces x(1-x)y(1-y) — bi-degree (2,2) — to round-off",
        e / amp < 1e-13, @sprintf("max error %.3e   amplitude %.4e", e, amp))
end

let (e, amp) = grid_error(EulerSquare(5, 3), p3_poly, 41)
    check("degree 3 reproduces x(1-x)y²(1-y) — bi-degree (3,3) — to round-off",
        e / amp < 1e-13, @sprintf("max error %.3e   amplitude %.4e", e, amp))
end

# THE CONTROL. The same bi-degree-(3,3) polynomial on the degree-2 space is not in it: a
# piecewise quadratic cannot be a global cubic. If this passed, "reproduced to round-off" would
# be a property of the resampler rather than of the space, and the two rows above would be
# vacuous.
let (e, amp) = grid_error(EulerSquare(6, 2), p3_poly, 41)
    check("degree 2 does NOT reproduce it — the control that must fail",
        e / amp > 1e-4, @sprintf("max error %.3e   relative %.3e   amplitude %.4e",
            e, e / amp, amp))
end

# =============================================================================================
header("3. the four edges vanish, because every basis function does")

# Random coefficients, not a projected function. A projection of something that already vanishes
# on ∂Ω would satisfy this whatever the space were; a random member of V_D satisfies it only
# because the basis is recombined.
function edge_maximum(sq, ω̂, N)
    Z = euler_grid(sq, ω̂, N)
    return maximum(abs, vcat(Z[1, :], Z[end, :], Z[:, 1], Z[:, end]))
end

let sq = EulerSquare(8, 2), ω̂ = randn(nbasis(sq.space))
    (lo, hi) = state_extrema(sq, ω̂)
    e = edge_maximum(sq, ω̂, 41)
    check("V_D: every edge is zero to round-off, for RANDOM degrees of freedom",
        e < 1e-13, @sprintf("max |ω| on ∂Ω = %.3e   ω ∈ [%.4f, %.4f] inside", e, lo, hi))
end

# THE CONTROL. The plain clamped basis of the `:free` space has functions that do not vanish on
# the boundary — that is the whole difference between the two state spaces — so the same
# measurement on it must be of order one. A resampler that clamped, wrapped or dropped its
# boundary samples would pass the row above and fail here.
let sq = EulerSquare(8, 2; state = :free), ω̂ = randn(nbasis(sq.space))
    (lo, hi) = state_extrema(sq, ω̂)
    e = edge_maximum(sq, ω̂, 41)
    check("V: the edges are NOT zero — the control that must fail",
        e > 0.1 * max(abs(lo), abs(hi)),
        @sprintf("max |ω| on ∂Ω = %.4f   ω ∈ [%.4f, %.4f] inside", e, lo, hi))
end

# =============================================================================================
header("4. Z[i,j] is ω(x_i, x_j) and not ω(x_j, x_i)")

# The test function carries B1's own widths -- w₁² = 0.01 against w₂² = 0.07, a factor of 2.6
# in the width and nothing else -- because that asymmetry is the only thing that makes this
# section possible: the radially symmetric bump a lazier test would reach for passes a
# transposed implementation exactly.
#
# It is B1's Gaussian TIMES sin(πx₁)sin(πx₂), and the factor is not decoration. B1's printed
# initial condition does not vanish on ∂Ω -- measured below -- so the max-norm error of its
# projection into V_D is a boundary mismatch that no mesh removes, and a convergence rate
# measured on it would be a rate of zero for a reason that has nothing to do with resampling.
# The factor makes the test function admissible while leaving the asymmetry intact, so the
# remaining error is projection error and the rate is a real tolerance.
const ω₀ = initial_condition(SECTION5_RUNS["b1"])
const f_asym = (x, y) -> sin(π * x) * sin(π * y) *
                         exp(-(x - 0.5)^2 / 0.01 - (y - 0.5)^2 / 0.07)

# REPORTED, not asserted: B1's Gaussian is 2.8 % of its peak at the midpoint of the x₂ = 0 edge,
# because w₂² = 0.07 is wide enough to reach it, while w₁² = 0.01 leaves 1e-11 on the x₁ edges.
# Every ω_h ∈ V_D is zero there, so the field map's initial panel shows the bump pinched to zero
# along two of the four edges. That is §5.4's own setup, not this reproduction's choice, and it
# is the same mismatch `interior_weights` exists to exclude for B3's reference fit.
check("B1's initial condition does NOT vanish on ∂Ω  [REPORTED]", true,
    @sprintf("ω₀(½,0) = %.4e   ω₀(0,½) = %.4e   peak = %.4f",
        ω₀(0.5, 0.0), ω₀(0.0, 0.5), ω₀(0.5, 0.5)))

# And the measurement that identifies it as a boundary mismatch rather than an interpolation
# defect, which is the distinction the row above is only half of: the max-norm error of the
# projected Gaussian is the SAME on every mesh and sits at the midpoint of an x₂ edge every time,
# because Z is exactly zero there and the error is exactly ω₀(½,0). An interpolation defect would
# converge; this cannot. The check asserts the non-convergence, which is why it is here and not
# in a scratch file.
function gaussian_edge_error(n)
    sq = EulerSquare(n, 2)
    û = project(sq.space, x -> ω₀(x[1], x[2]))
    N = 41
    Z = euler_grid(sq, û, N)
    xs = euler_axis(N)
    E = [abs(Z[i, j] - ω₀(xs[i], xs[j])) for i in 1:N, j in 1:N]
    k = argmax(E)
    return (E[k], xs[k[1]], xs[k[2]])
end

let (e16, x16, y16) = gaussian_edge_error(16), (e32, _, _) = gaussian_edge_error(32),
    (e64, x64, y64) = gaussian_edge_error(64)

    check("and its projection error is MESH-INDEPENDENT, on an x₂ edge",
        abs(e64 - e16) / e16 < 1e-12 && y16 in (0.0, SQUARE_LENGTH) &&
            y64 in (0.0, SQUARE_LENGTH),
        @sprintf("16 cells: %.4e at (%.2f, %.2f)   32: %.4e   64: %.4e at (%.2f, %.2f)",
            e16, x16, y16, e32, e64, x64, y64))
end

function asym_errors(n, p, N)
    sq = EulerSquare(n, p)
    û = project(sq.space, x -> f_asym(x[1], x[2]))
    Z = euler_grid(sq, û, N)
    xs = euler_axis(N)
    direct = maximum(abs(Z[i, j] - f_asym(xs[i], xs[j])) for i in 1:N, j in 1:N)
    swapped = maximum(abs(Z[i, j] - f_asym(xs[j], xs[i])) for i in 1:N, j in 1:N)
    return (direct, swapped)
end

let (d, s) = asym_errors(64, 2, 129)
    check("the direct comparison is small and the transposed one is order one",
        s > 100d, @sprintf("direct %.3e   transposed %.3e   ratio %.0fx", d, s, s / d))
    # The amplitude the transposed error has to be measured against: the test function's peak is
    # 1, so a transposed read is wrong by a fraction of the field itself and not by a residual.
    check("the transposed error is a fraction of the peak, not a residual", s > 0.2,
        @sprintf("transposed %.4f   peak = %.4f", s, f_asym(0.5, 0.5)))
end

# The other side of the two-sided check: the direct error is PROJECTION error and must therefore
# fall at the order of the space when the mesh is refined. A transposed or shifted resampling
# would leave a floor that does not converge, which is what makes a rate a sharper statement
# than any absolute threshold.
#
# Degree 2 is the runs' own degree, and `O(h^{p+1})` predicts 8x per halving of h. That IS the
# threshold, with no margin added: the law is the tolerance. Both measured factors come out
# ABOVE it -- 12x and 33x -- because w₁ = 0.1 against h = 1/16 leaves the coarsest mesh short of
# the asymptotic regime, so the error is still shedding the pre-asymptotic term as well as the
# leading one. A rate that fell BELOW 8 would say the resampling, not the projection, was
# setting the floor.
let (e16, _) = asym_errors(16, 2, 129), (e32, _) = asym_errors(32, 2, 129),
    (e64, _) = asym_errors(64, 2, 129)

    check_refined("the direct error falls at the degree-2 rate, 16 → 32 cells", e16, e32;
        atol = 1e-14, minrate = 8.0)
    check_refined("and again, 32 → 64 cells", e32, e64; atol = 1e-14, minrate = 8.0)
end

summary("verify_euler_grid.jl")
