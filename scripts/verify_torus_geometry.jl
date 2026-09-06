#!/usr/bin/env julia
#
# The contour geometry of Section 4.1, checked before any solver is built on it.
#
#     julia --project=scripts scripts/verify_torus_geometry.jl
#
# Three independent statements, each of which would silently corrupt A1 if it were wrong:
#
#   1. the closed form l_h = pi / (sqrt(h) agm(1,sqrt(h))) agrees with quadrature of the
#      arclength integral `eq:relaxation-time` actually writes;
#   2. l_h is the PERIOD of the X_h flow, which is the manuscript's own characterisation and
#      is what makes tau_h a relaxation time rather than a geometric curiosity;
#   3. the contour average u_infinity is the t -> infinity limit of parallel diffusion, i.e.
#      it is invariant along the X_h flow and reproduces the exact answer on a field that is
#      already a function of h.
#
# The control that can fail is in section 1: a WRONG closed form -- the same expression with
# agm(1,h) instead of agm(1,sqrt(h)) -- must disagree with the quadrature by order one. Without
# it, "the two agree" would only say that one of them was used to write the other.

using MetriplecticRelaxation
using MetriplecticRelaxation: agm, contour_length, contour_length_quadrature,
                              contour_average, relaxation_time, islands_h, Gaussian,
                              CENTRAL_ISLANDS, _gauss_legendre
using Printf

include(joinpath(@__DIR__, "check.jl"))
using .Checks: header, check, summary

# Down to 1e-5, i.e. contours a factor 10^5 closer to the separatrix than to the island
# centre. The tolerance is 1e-13 rather than a loose one because the parameterisation is
# cancellation-free: see `_contour_point`. Before that fix the h = 1e-4 residual was 1.3e-10
# and GREW with the node count, which is what identified it as round-off.
const HS = (0.9, 0.7, 0.5, 0.3, 0.1, 0.01, 1e-4, 1e-5)

# =============================================================================================
header("1. l_h = pi / (sqrt(h) agm(1,sqrt(h))) against quadrature of the arclength integral")

for h in HS
    exact = contour_length(h)
    quad = contour_length_quadrature(h; n = 400)
    check(@sprintf("h = %-8g   closed form vs quadrature", h),
        abs(exact - quad) / exact < 1e-13,
        @sprintf("%.12f vs %.12f   rel %.2e", exact, quad, abs(exact - quad) / exact))
end

# The control. If this passed, section 1 would be vacuous.
let h = 0.5
    wrong = π / (sqrt(h) * agm(1.0, h))
    quad = contour_length_quadrature(h; n = 400)
    check("CONTROL: agm(1,h) instead of agm(1,sqrt(h)) must FAIL",
        abs(wrong - quad) / quad > 1e-3,
        @sprintf("%.6f vs %.6f   rel %.2e", wrong, quad, abs(wrong - quad) / quad))
end

# =============================================================================================
header("2. l_h is the period of the X_h flow")

# X_h = (d_2 h, -d_1 h) for h = cos^2 x_1 sin^2 x_2, integrated with classical RK4 from a
# point on the contour. After exactly l_h the orbit must close.
∂₁h(x₁, x₂) = -2cos(x₁) * sin(x₁) * sin(x₂)^2
∂₂h(x₁, x₂) = 2cos(x₁)^2 * sin(x₂) * cos(x₂)
Xh(x) = (∂₂h(x[1], x[2]), -∂₁h(x[1], x[2]))

function flow(x, Δt, nsteps)
    for _ in 1:nsteps
        k1 = Xh(x)
        k2 = Xh((x[1] + Δt / 2 * k1[1], x[2] + Δt / 2 * k1[2]))
        k3 = Xh((x[1] + Δt / 2 * k2[1], x[2] + Δt / 2 * k2[2]))
        k4 = Xh((x[1] + Δt * k3[1], x[2] + Δt * k3[2]))
        x = (x[1] + Δt / 6 * (k1[1] + 2k2[1] + 2k3[1] + k4[1]),
            x[2] + Δt / 6 * (k1[2] + 2k2[2] + 2k3[2] + k4[2]))
    end
    return x
end

for h in (0.9, 0.7, 0.5, 0.3, 0.1)
    # Start on the contour at s = 0, i.e. directly above the island centre.
    c = CENTRAL_ISLANDS[1]
    x0 = (c[1], c[2] + acos(sqrt(h)))
    ℓ = contour_length(h)
    n = 20000
    x1 = flow(x0, ℓ / n, n)
    d = hypot(x1[1] - x0[1], x1[2] - x0[2])
    check(@sprintf("h = %-6g   orbit closes after xi = l_h", h), d < 1e-8,
        @sprintf("|x(l_h) - x(0)| = %.3e   (l_h = %.6f)", d, ℓ))
end

# The same statement read as a control: after HALF a period the orbit must be far from where
# it started, or "it closes" would hold of a flow that never moved.
let h = 0.5, c = CENTRAL_ISLANDS[1]
    x0 = (c[1], c[2] + acos(sqrt(h)))
    ℓ = contour_length(h)
    xh = flow(x0, ℓ / 2 / 20000, 20000)
    d = hypot(xh[1] - x0[1], xh[2] - x0[2])
    check("CONTROL: after l_h/2 the orbit is elsewhere", d > 0.1,
        @sprintf("|x(l_h/2) - x(0)| = %.4f", d))
end

# =============================================================================================
header("3. the contour average is the t -> infinity limit")

# (a) A field that is already a function of h is its own contour average, exactly.
for h in (0.8, 0.5, 0.2, 0.05)
    f(x₁, x₂) = 3islands_h(x₁, x₂)^2 - islands_h(x₁, x₂) + 1
    avg = contour_average(f, h, CENTRAL_ISLANDS[1]; n = 400)
    exact = 3h^2 - h + 1
    check(@sprintf("h = %-6g   u_0 = u_0(h) is its own average", h),
        abs(avg - exact) < 1e-11, @sprintf("%.12f vs %.12f", avg, exact))
end

# (b) For a general field the average is invariant along the flow, which is the property that
#     makes it the limit: parallel diffusion equalises along contours and does nothing else.
let h = 0.4, c = CENTRAL_ISLANDS[1]
    g = Gaussian((π, π + 0.1), (0.25, 0.4), 2π * 0.25 * 0.4)
    u₀(x₁, x₂) = g(x₁, x₂)
    a = contour_average(u₀, h, c; n = 400)
    # Same contour, sampled from a rotated starting point: the average cannot depend on it.
    ℓ = contour_length(h)
    x0 = (c[1], c[2] + acos(sqrt(h)))
    xs = flow(x0, ℓ / 3 / 5000, 5000)
    # Averaging u_0 composed with the flow over one period is the same integral, discretised
    # by the flow's own uniform time grid rather than by the theta substitution.
    n = 4000
    acc = 0.0
    x = xs
    for _ in 1:n
        acc += u₀(x[1], x[2])
        x = flow(x, ℓ / n, 1)
    end
    check("contour average = flow-time average of u_0", abs(a - acc / n) / abs(a) < 1e-6,
        @sprintf("%.10f vs %.10f   rel %.2e", a, acc / n, abs(a - acc / n) / abs(a)))
end

# =============================================================================================
header("4. tau_h and its two limits")

check("tau_h -> 1/4 at the island centre", abs(relaxation_time(1.0) - 0.25) < 1e-14,
    @sprintf("tau_1 = %.16f", relaxation_time(1.0)))

# The harmonic approximation h = 1 - s^2 - t^2 rotates at angular frequency 2, hence period pi.
check("l_h -> pi at the island centre", abs(contour_length(1.0) - π) < 1e-14,
    @sprintf("l_1 = %.16f   pi = %.16f", contour_length(1.0), π))

check("tau_h diverges at the separatrix", relaxation_time(1e-8) > 1e6,
    @sprintf("tau(1e-8) = %.4e", relaxation_time(1e-8)))

check("tau_h is monotone decreasing in h",
    all(relaxation_time(HS[i]) < relaxation_time(HS[i + 1]) for i in 1:(length(HS) - 1)),
    @sprintf("tau(0.9) = %.4f ... tau(1e-4) = %.4e", relaxation_time(0.9),
        relaxation_time(1e-4)))

# =============================================================================================
header("5. the Gauss-Legendre rule used above")

let (x, w) = _gauss_legendre(64, -1.0, 1.0)
    check("weights sum to the interval length", abs(sum(w) - 2) < 1e-14,
        @sprintf("%.16f", sum(w)))
    # Exact for polynomials of degree 2n-1; degree 127 is the last one it must get right.
    check("integrates x^126 exactly", abs(sum(w .* x .^ 126) - 2 / 127) < 1e-15,
        @sprintf("%.6e vs %.6e", sum(w .* x .^ 126), 2 / 127))
end

summary("verify_torus_geometry.jl")
