#!/usr/bin/env julia
#
# The uniform resampling the Section 5.5 field maps are drawn on.
#
#     julia --project=scripts scripts/verify_gradshafranov_grid.jl
#
# `gs_grid` evaluates a Section 5.5 state off the space's own quadrature grid onto a uniform one,
# and `gs_ordinate_grid` weights it by sigma(r) to give the field the manuscript's colour plot
# actually shows.  Both are a handful of lines and both are the kind of code that looks right and
# is silently wrong, so every claim they make is settled here in closed form rather than by
# looking at the picture.
#
# This is `verify_euler_grid.jl` for the other geometry, and the differences from it are the
# point rather than an accident:
#
#   * the domain is [1,7] x [-9.5,9.5], not the unit square, so the two axes need separate node
#     counts and the script says which interval each covers;
#   * the measure is dmu = dr dz / r, not dx, so the independent statement about the domain
#     integrates in dmu -- a missing 1/r is exactly the defect a shared resampler would produce;
#   * the colour field is sigma(r) j and not the state, so the radial weight has its own section.
#
# Six sections, in the order a defect in one would invalidate the next:
#
#   1. the axes: the requested node counts, BOTH endpoints on each, uniform spacing, the right
#      interval per axis, and the Gaussian's centre on a node -- then the domain checked
#      independently by integrating a resampled field with the trapezoidal rule IN dmu and
#      comparing against the space's own quadrature, which shares no code with `gs_axes`;
#   2. exactness: a polynomial the space contains is reproduced to ROUND-OFF, at both degrees,
#      with a polynomial it does not contain as the control that says the claim is about the
#      degree rather than something the resampler returns for anything;
#   3. the homogeneous-Dirichlet boundary: all four edges vanish to round-off, for RANDOM
#      coefficients rather than for a projected function that happens to vanish there -- with the
#      `:free` space as the control that MUST fail;
#   4. the index convention, on equal node counts.  The figure samples 61x191, where a transposed
#      read is a dimension mismatch rather than a plausible picture -- a property of the sample
#      counts and not of the mesh; the convention and not the shape is what is under test, so the
#      row is measured where the comparison can be formed at all.
#   5. sigma(r) is applied along the RADIAL index.  The ratio of the two grids must be constant
#      along z, must equal sigma(r_i) along r, and must actually vary -- a weight that did not
#      would make the section vacuous;
#   6. C1's own initial condition, reported rather than asserted.
#
# Everything here runs on spaces of at most 2601 degrees of freedom and the whole script is
# seconds.

using MetriplecticRelaxation
using MetriplecticRelaxation: SECTION55_RUNS, GradShafranovBox, GS_RADIAL, GS_AXIAL,
                              gs_axes, gs_grid, gs_ordinate_grid, gs_entropy_weight,
                              gs_density, gs_state, initial_condition
using GeometricBrackets: nbasis, project, field, quadrature_nodes, quadrature_weights
using Printf
using Random

include(joinpath(@__DIR__, "check.jl"))
using .Checks: header, check, check_refined, summary

Random.seed!(0x5f1a20c3)

println("verify_gradshafranov_grid.jl  --  the uniform resampling behind the §5.5 field maps")

# The sample counts the figure script uses. Odd on both axes, and chosen so the two spacings
# match EXACTLY: 6/60 and 19/190 are both 0.1. An even count on either axis would move the
# Gaussian's centre off the grid, which section 1 asserts.
const NR = 61
const NZ = 191

const RADIAL_LENGTH = GS_RADIAL[2] - GS_RADIAL[1]
const AXIAL_LENGTH = GS_AXIAL[2] - GS_AXIAL[1]

# =============================================================================================
header("1. the axes are Nr and Nz nodes over their own intervals, endpoints included")

let (rs, zs) = gs_axes(NR, NZ)
    check("the axes have the requested node counts — two counts, not one",
        length(rs) == NR && length(zs) == NZ,
        @sprintf("length(rs) = %d   requested %d      length(zs) = %d   requested %d",
            length(rs), NR, length(zs), NZ))

    check("the radial axis spans GS_RADIAL, BOTH endpoints",
        rs[1] == GS_RADIAL[1] && rs[end] == GS_RADIAL[2],
        @sprintf("r₁ = %.17g   r_N = %.17g   GS_RADIAL = (%.17g, %.17g)",
            rs[1], rs[end], GS_RADIAL...))

    check("the axial axis spans GS_AXIAL, BOTH endpoints",
        zs[1] == GS_AXIAL[1] && zs[end] == GS_AXIAL[2],
        @sprintf("z₁ = %.17g   z_N = %.17g   GS_AXIAL = (%.17g, %.17g)",
            zs[1], zs[end], GS_AXIAL...))

    let dr = diff(rs), hr = RADIAL_LENGTH / (NR - 1)
        check("the radial spacing is uniform at (r₂-r₁)/(Nr-1)",
            maximum(abs, dr .- hr) < 1e-14,
            @sprintf("max |Δr - h| = %.3e   h = %.17g", maximum(abs, dr .- hr), hr))
    end
    let dz = diff(zs), hz = AXIAL_LENGTH / (NZ - 1)
        check("the axial spacing is uniform at (z₂-z₁)/(Nz-1)",
            maximum(abs, dz .- hz) < 1e-14,
            @sprintf("max |Δz - h| = %.3e   h = %.17g", maximum(abs, dz .- hz), hz))
    end

    # Not decoration, and the same reason `verify_euler_grid.jl` gives: C1's Gaussian is centred
    # at (r₀, z₀) = (4, 0), which is the midpoint of both intervals. An even count would straddle
    # it, the initial panel would under-report its own peak, and `figure_fields` shares its colour
    # range between the two panels — so the final panel would be drawn against a peak that is not
    # the run's.
    # The centre is read off the run rather than copied in, so that a change to `SECTION55_RUNS`
    # reddens this row instead of leaving it asserting a literal the run no longer has.
    let (r₀, z₀) = SECTION55_RUNS["c1"].gaussian.x₀
        check("C1's Gaussian centre is a node on both axes, for odd counts",
            any(==(r₀), rs) && any(==(z₀), zs) && isodd(NR) && isodd(NZ),
            @sprintf("Nr = %d, Nz = %d, both odd: %s   centre (%.4g, %.4g) on the grid: %s",
                NR, NZ,
                string(isodd(NR) && isodd(NZ)), r₀, z₀,
                string(any(==(r₀), rs) && any(==(z₀), zs))))
    end
end

# The independent statement about which intervals the grid covers, and the one place a missing
# 1/r would show. If either axis ran over the wrong interval, or the measure were taken as dx,
# the trapezoidal integral of a resampled field would not agree with the space's own quadrature
# — and that comparison shares no code with `gs_axes`.
let box = GradShafranovBox((18, 21), 2), spec = SECTION55_RUNS["c1"],
    ĵ = gs_state(box, spec)

    Z = gs_grid(box, ĵ, NR, NZ)
    (rs, _) = gs_axes(NR, NZ)
    hr = RADIAL_LENGTH / (NR - 1)
    hz = AXIAL_LENGTH / (NZ - 1)
    wr = fill(1.0, NR)
    wr[1] = wr[end] = 0.5
    wz = fill(1.0, NZ)
    wz[1] = wz[end] = 0.5
    trap = hr * hz * sum(wr[i] * wz[j] * Z[i, j] / rs[i] for i in 1:NR, j in 1:NZ)

    s = box.space
    wμ = quadrature_weights(s) .* gs_density.(quadrature_nodes(s))
    exact = sum(wμ .* field(s, ĵ, (0, 0)))

    e = abs(trap - exact) / abs(exact)
    check("the trapezoidal ∫j dμ of the resampled field IS the space's own ∫j dμ",
        e < 2e-4,
        @sprintf("trapezoid %.10e   ∫j dμ = %.10e   relative %.3e", trap, exact, e))
end

# =============================================================================================
header("2. a polynomial the space contains is reproduced to round-off")

# The homogeneous-Dirichlet space of degree p contains, per axis, exactly the polynomials of
# degree ≤ p that vanish at both ends of THAT axis — and the two axes have different ends here.
# At p = 2 each axis contributes a one-dimensional span, so the product is fixed up to a scale
# and settles exactness while saying nothing about the index convention, which is section 4's
# business. At p = 3 the axial factor can be asymmetric, and is.
const p2_poly = (r, z) -> (r - 1) * (7 - r) * (z + 9.5) * (9.5 - z)
const p3_poly = (r, z) -> (r - 1) * (7 - r) * (z + 9.5)^2 * (9.5 - z)

function grid_error(box, f, Nr, Nz)
    û = project(box.space, x -> f(x[1], x[2]))
    Z = gs_grid(box, û, Nr, Nz)
    (rs, zs) = gs_axes(Nr, Nz)
    amp = maximum(abs(f(r, z)) for r in rs, z in zs)
    return (maximum(abs(Z[i, j] - f(rs[i], zs[j])) for i in 1:Nr, j in 1:Nz), amp)
end

let (e, amp) = grid_error(GradShafranovBox((6, 8), 2), p2_poly, 41, 41)
    check("degree 2 reproduces (r-1)(7-r)(z+9.5)(9.5-z) — bi-degree (2,2) — to round-off",
        e / amp < 1e-13, @sprintf("max error %.3e   amplitude %.4e", e, amp))
end

let (e, amp) = grid_error(GradShafranovBox((5, 7), 3), p3_poly, 41, 41)
    check("degree 3 reproduces (r-1)(7-r)(z+9.5)²(9.5-z) — bi-degree (2,3) — to round-off",
        e / amp < 1e-13, @sprintf("max error %.3e   amplitude %.4e", e, amp))
end

# THE CONTROL. The same polynomial, cubic in z, on the degree-2 space is not in it: a piecewise
# quadratic cannot be a global cubic. If this passed, "reproduced to round-off" would be a
# property of the resampler rather than of the space, and the two rows above would be vacuous.
let (e, amp) = grid_error(GradShafranovBox((6, 8), 2), p3_poly, 41, 41)
    check("degree 2 does NOT reproduce it — the control that must fail",
        e / amp > 1e-4, @sprintf("max error %.3e   relative %.3e   amplitude %.4e",
            e, e / amp, amp))
end

# =============================================================================================
header("3. the four edges vanish, because every basis function does")

# Random coefficients, not a projected function. A projection of something that already vanishes
# on ∂Ω would satisfy this whatever the space were; a random member of V_D satisfies it only
# because the basis is recombined.
function edge_maximum(box, ĵ, Nr, Nz)
    Z = gs_grid(box, ĵ, Nr, Nz)
    return maximum(abs, vcat(Z[1, :], Z[end, :], Z[:, 1], Z[:, end]))
end

let box = GradShafranovBox((8, 9), 2), ĵ = randn(nbasis(box.space))
    (lo, hi) = extrema(field(box.space, ĵ, (0, 0)))
    e = edge_maximum(box, ĵ, 41, 41)
    check("V_D: every edge is zero to round-off, for RANDOM degrees of freedom",
        e < 1e-13, @sprintf("max |j| on ∂Ω = %.3e   j ∈ [%.4f, %.4f] inside", e, lo, hi))
end

# THE CONTROL. The plain clamped basis of the `:free` space has functions that do not vanish on
# the boundary — that is the whole difference between the two state spaces, and the difference
# C2's rim exists to remove — so the same measurement on it must be of order one. A resampler
# that clamped, wrapped or dropped its boundary samples would pass the row above and fail here.
let box = GradShafranovBox((8, 9), 2; state = :free), ĵ = randn(nbasis(box.space))
    (lo, hi) = extrema(field(box.space, ĵ, (0, 0)))
    e = edge_maximum(box, ĵ, 41, 41)
    check("V: the edges are NOT zero — the control that must fail",
        e > 0.1 * max(abs(lo), abs(hi)),
        @sprintf("max |j| on ∂Ω = %.4f   j ∈ [%.4f, %.4f] inside", e, lo, hi))
end

# =============================================================================================
header("4. Z[i,j] is j(r_i, z_j) and not j(r_j, z_i)")

# C1 is 18x21 cells and the figure samples 61x191, so a transposed `gs_grid` there returns a
# matrix of the wrong shape and dies on the next operation rather than drawing a wrong picture.
# That is a happy accident of the geometry and not a check, so the convention is measured where
# the transposed comparison can be formed at all: a SQUARE mesh and equal node counts. What is
# under test is the axis order inside `gs_grid`, which does not know the shape it was called with.
#
# The test function carries C1's own widths — w₁² = 0.5 against w₂² = 3.2 — times a sine
# envelope. The envelope is not decoration: the max-norm error of a projection into V_D is a
# boundary mismatch for any target that does not vanish on ∂Ω, and a convergence rate measured on
# one would be a rate of zero for a reason that has nothing to do with resampling. Section 6
# reports how much C1's printed Gaussian misses by, which is little, but "little" is not
# "admissible".
#
# THE CENTRE IS OFF C1's, AND IT HAS TO BE. C1's Gaussian sits at (r₀, z₀) = (4, 0), which is the
# midpoint of BOTH intervals, so on equal node counts it lands on the same index on both axes —
# and a transposition leaves the peak exactly where it found it. At C1's own centre the
# transposed error is 0.1541 against a peak of 1, because the discrepancy is pushed off-centre
# and only the shoulders disagree, and the row below does not hold at 0.2. What that asks for is
# a stronger control and not a looser bound: at (3, 2) the two centres fall on different indices,
# the transposition displaces the bump bodily, and the transposed error reaches the peak exactly
# — 0.8249 against a peak of 0.8249, a ratio of 1.00, and 965x the direct error against 148x at
# C1's own centre. The widths still differ, so the asymmetry the section needs is untouched.
const f_asym = (r, z) -> sin(π * (r - 1) / 6) * sin(π * (z + 9.5) / 19) *
                         exp(-(r - 3)^2 / 0.5 - (z - 2)^2 / 3.2)

function asym_errors(n, p, N)
    box = GradShafranovBox((n, n), p)
    û = project(box.space, x -> f_asym(x[1], x[2]))
    Z = gs_grid(box, û, N, N)
    (rs, zs) = gs_axes(N, N)
    direct = maximum(abs(Z[i, j] - f_asym(rs[i], zs[j])) for i in 1:N, j in 1:N)
    swapped = maximum(abs(Z[j, i] - f_asym(rs[i], zs[j])) for i in 1:N, j in 1:N)
    return (direct, swapped)
end

let (d, s) = asym_errors(32, 2, 65)
    check("the direct comparison is small and the transposed one is order one",
        s > 100d, @sprintf("direct %.3e   transposed %.3e   ratio %.0fx", d, s, s / d))
    # The amplitude the transposed error has to be measured against: a transposed read is wrong
    # by a fraction of the field itself and not by a residual.
    let peak = maximum(abs(f_asym(r, z))
        for r in range(GS_RADIAL...; length = 65),
        z in range(GS_AXIAL...; length = 65))
        check("the transposed error is a fraction of the peak, not a residual",
            s > 0.2 * peak,
            @sprintf("transposed %.4f   peak = %.4f   ratio %.2f", s, peak, s / peak))
    end
end

# The other side of the two-sided check: the direct error is PROJECTION error and must therefore
# fall at the order of the space when the mesh is refined. A transposed or shifted resampling
# would leave a floor that does not converge, which is what makes a rate a sharper statement than
# any absolute threshold.
#
# Degree 2 is C1's own degree, and `O(h^{p+1})` predicts 8x per halving of h. That IS the
# threshold, with no margin added: the law is the tolerance. A rate that fell BELOW 8 would say
# the resampling, not the projection, was setting the floor.
let (e8, _) = asym_errors(8, 2, 65), (e16, _) = asym_errors(16, 2, 65),
    (e32, _) = asym_errors(32, 2, 65)

    check_refined("the direct error falls at the degree-2 rate, 8 → 16 cells", e8, e16;
        atol = 1e-14, minrate = 8.0)
    check_refined("and again, 16 → 32 cells", e16, e32; atol = 1e-14, minrate = 8.0)
end

# =============================================================================================
header("5. σ(r) is applied along the radial index")

# `gs_ordinate_grid` is what the colour panels are drawn from, and the manuscript says so: the
# plot is u/(Cr²+D) = σ(r) j rather than the state, because only in that variable are the
# contours of the relaxed field the contours of ψ. σ depends on r ALONE, so putting it on the
# wrong index is not a small error and is also not a shape error — both grids are Nr x Nz — which
# is why it gets its own section rather than riding on section 4.
let box = GradShafranovBox((18, 21), 2), ĵ = randn(nbasis(box.space))
    Z = gs_grid(box, ĵ, NR, NZ)
    Y = gs_ordinate_grid(box, ĵ, NR, NZ)
    (rs, _) = gs_axes(NR, NZ)

    # Interior only: every edge of Z is zero, so the ratio there is 0/0.
    ii, jj = 2:(NR - 1), 2:(NZ - 1)
    ratio = [Y[i, j] / Z[i, j] for i in ii, j in jj]

    # (a) constant along z, which is what "depends on r alone" means pointwise.
    spread = maximum(maximum(abs, r .- r[1]) for r in eachrow(ratio))
    check("the weight is constant along the axial index", spread < 1e-12,
        @sprintf("max variation of Y/Z along z = %.3e", spread))

    # (b) and it is σ(r_i), not σ of something else.
    e = maximum(abs(ratio[k, 1] - gs_entropy_weight(rs[i])) for (k, i) in enumerate(ii))
    check("and it equals σ(r_i) on the radial index", e < 1e-12,
        @sprintf("max |Y/Z - σ(r)| = %.3e   σ(1) = %.6f   σ(7) = %.6f",
            e, gs_entropy_weight(1.0), gs_entropy_weight(7.0)))

    # (c) THE CONTROL. σ has to actually vary over [1,7], or (a) and (b) would both hold for a
    # resampler that applied a constant, and neither would say anything about the index.
    let lo = minimum(gs_entropy_weight, rs), hi = maximum(gs_entropy_weight, rs)
        check(
            "σ VARIES over the radial interval — otherwise the two rows above are vacuous",
            hi / lo > 1.5, @sprintf("σ ∈ [%.6f, %.6f]   ratio %.3f", lo, hi, hi / lo))
    end
end

# =============================================================================================
header("6. C1's initial condition on the boundary  [REPORTED]")

# §5.4's counterpart is a finding: B1's printed Gaussian is 2.8 % of its peak at the midpoint of
# an edge, so its projection into V_D carries a boundary mismatch no mesh removes. C1's is not —
# w₁² = 0.5 on a domain of width 6 and w₂² = 3.2 on one of width 19 both leave the Gaussian at
# round-off on ∂Ω. Reported and not asserted, because it is a property of §5.5's printed setup
# rather than of the resampling; it is here so that the contrast with B1 is measured rather than
# assumed, and so that section 4's sine envelope is visibly a precaution and not a fudge.
let ω₀ = initial_condition(SECTION55_RUNS["c1"])
    check("C1's initial condition is at round-off on ∂Ω, unlike B1's  [REPORTED]", true,
        @sprintf("j₀(1,0) = %.3e   j₀(7,0) = %.3e   j₀(4,±9.5) = %.3e   peak = %.4f",
            ω₀(1.0, 0.0), ω₀(7.0, 0.0), ω₀(4.0, 9.5), ω₀(4.0, 0.0)))
end

summary("verify_gradshafranov_grid.jl")
