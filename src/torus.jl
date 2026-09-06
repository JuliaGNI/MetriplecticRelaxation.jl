#
# Section 4 of *Metriplectic relaxation to equilibria*: the periodic torus.
#
# This file is the problem definition and nothing else -- initial conditions, the prescribed
# Hamiltonian of the analytic test case, the closed-form entropy minimisers, and the contour
# geometry behind the relaxation time. It is deliberately free of any discretisation, so that
# `spectral.jl` and `spline.jl` are two independent solvers of the *same* stated problem and
# their agreement means something.
#

@doc raw"""
The side length of the periodic domain ``\Omega = \mathbb{T}^2 = [0, 2\pi]^2`` of Section 4.
"""
const DOMAIN_LENGTH = 2π

@doc raw"""
The area ``|\Omega| = 4\pi^2``, which divides the total mass to give the mean value
``u_\Omega``.
"""
const DOMAIN_AREA = 4π^2

## Initial conditions

@doc raw"""
    Gaussian(x₀, w, N)

The anisotropic Gaussian of `eq:initial_gaussian`,

```math
u_G(x) = \frac{1}{N} \exp\Big[
    - \frac{(x_1 - x_{0,1})^2}{w_1^2} - \frac{(x_2 - x_{0,2})^2}{w_2^2} \Big] .
```

Callable as `g(x₁, x₂)`. `images` selects between the two readings below; it is `0` by
default, which is the formula exactly as printed.

# The formula is not periodic, and for run A4 that is not a rounding detail

The manuscript writes this formula unmodified on ``\mathbb{T}^2``, with no summation over
periodic images. Taken literally it is therefore *discontinuous* across the boundary, by the
value it still has at the far edge of the domain. How much that matters depends entirely on
how far the centre sits from the wrap, and it differs by two orders of magnitude between the
four runs:

| run | ``x_{0,2}`` | ``w_2`` | amplitude | value at the wrap | jump / peak |
|:--|--:|--:|--:|--:|--:|
| A1 | ``\pi + 0.1`` | 0.25 | 1.59 | ``10^{-69}`` | ``10^{-69}`` |
| A2, A3 | ``\pi`` | 1 | 1 | ``5.2 \times 10^{-5}`` | ``5.2 \times 10^{-5}`` |
| **A4** | ``3\pi/2`` | 1 | **1.8** | **0.153** | **8.5 %** |

A4's centre is only ``\pi/2`` from the boundary, so ``1.8 \, e^{-(\pi/2)^2} = 0.153``. That
is a genuine discontinuity in the stated initial condition, not a discretisation artefact,
and the two discretisations resolve it differently — measured, the spectral and spline
initial states disagree by ``6.9 \times 10^{-2}`` for A4 against ``2.5 \times 10^{-4}`` for
A2, whose Gaussian is otherwise identical.

`images = n` sums over the periodic images ``x_0 + 2\pi k`` for ``k \in -n:n`` in each
direction, which is the reading under which ``u_G`` is a function on ``\mathbb{T}^2`` at all.
One image either way already converges to ``10^{-16}`` for every parameter set here.

Which reading the manuscript intends is not stated. Both are therefore run for A4 and both
are reported; `run_a4.jl` uses the literal form as the primary, being what is printed, and
the periodised form as the control that identifies the discontinuity as the cause of the
disagreement rather than the solver.
"""
struct Gaussian{T}
    x₀::NTuple{2, T}
    w::NTuple{2, T}
    N::T
    images::Int
end

# `x₀ = (π, π + 0.1)` is a `Tuple{Irrational, Float64}`, so the parameters have to be promoted
# rather than required to arrive already matching.
function Gaussian(x₀::Tuple, w::Tuple, N; images::Int = 0)
    T = promote_type(map(typeof ∘ float, (x₀..., w..., N))...)
    Gaussian{T}(map(T, x₀), map(T, w), T(N), images)
end

function (g::Gaussian{T})(x₁, x₂) where {T}
    acc = zero(T)
    for k₂ in (-g.images):(g.images), k₁ in (-g.images):(g.images)

        d₁ = x₁ - g.x₀[1] - k₁ * DOMAIN_LENGTH
        d₂ = x₂ - g.x₀[2] - k₂ * DOMAIN_LENGTH
        acc += exp(-d₁^2 / g.w[1]^2 - d₂^2 / g.w[2]^2)
    end
    return acc / g.N
end

(g::Gaussian)(x) = g(x[1], x[2])

"""
    periodise(g::Gaussian, images = 2)
    periodise(spec::RunSpec, images = 2)

The same object with the Gaussian summed over `images` periodic images in each direction —
the reading under which the initial condition is a function on the torus rather than a
formula restricted to ``[0,2\\pi]^2``. See [`Gaussian`](@ref).
"""
periodise(g::Gaussian, images::Int = 2) = Gaussian(g.x₀, g.w, g.N; images = images)

@doc raw"""
    islands_h(x₁, x₂)

The prescribed Hamiltonian ``h(x) = \cos^2(x_1) \sin^2(x_2)`` of `eq:islands-h`.

Its contours form a periodic array of four islands on ``\Omega``, centred at
[`ISLAND_CENTRES`](@ref) where ``h = 1``, and separated by the separatrix set ``h = 0``,
which is the union of the lines ``x_1 \in \{\pi/2, 3\pi/2\}`` and
``x_2 \in \{0, \pi, 2\pi\}``. On the separatrices ``\nabla h = 0``, hence ``X_h = 0``, and
the parallel diffusion of `eq:parallel-diffusion` stalls — which is what makes
[`relaxation_time`](@ref) diverge there.
"""
islands_h(x₁, x₂) = cos(x₁)^2 * sin(x₂)^2

islands_h(x) = islands_h(x[1], x[2])

@doc raw"""
The centres of the four islands of [`islands_h`](@ref), where ``h`` attains its maximum
``1``.

The two with ``x_1 = \pi`` are the "two central (full) islands" of Fig. 1: they lie wholly
inside ``[0, 2\pi]^2`` rather than being cut by the boundary, and they are the two the
manuscript evaluates ``\tau_h`` on. Run A1's Gaussian is centred at ``(\pi, \pi + 0.1)``,
just above the separatrix ``x_2 = \pi`` that divides them, which is what gives the scatter
plot its two distinct branches.
"""
const ISLAND_CENTRES = ((float(π), π / 2), (float(π), 3π / 2),
    (0.0, π / 2), (0.0, 3π / 2))

@doc raw"""
The two islands of [`ISLAND_CENTRES`](@ref) that lie wholly inside ``[0, 2\pi]^2``.
"""
const CENTRAL_ISLANDS = (ISLAND_CENTRES[1], ISLAND_CENTRES[2])

## Contour geometry: the relaxation time

@doc raw"""
    agm(a, b)

The arithmetic-geometric mean, iterated to convergence.

Present because [`contour_length`](@ref) needs the complete elliptic integral of the first
kind, ``K(m) = \pi / \big(2 \, \mathrm{agm}(1, \sqrt{1-m})\big)``, and the iteration
converges quadratically — six steps reach `eps()` for every argument used here. Writing it
out avoids a dependency on a special-function package for one closed form.
"""
function agm(a::T, b::T) where {T <: AbstractFloat}
    while abs(a - b) > 4eps(T) * abs(a)
        a, b = (a + b) / 2, sqrt(a * b)
    end
    return (a + b) / 2
end

agm(a, b) = agm(promote(float(a), float(b))...)

@doc raw"""
    contour_length(h)

The reduced contour length ``\ell_h = \oint_{C_h} ds / |\nabla h|`` of `eq:relaxation-time`,
on a contour of [`islands_h`](@ref) inside one island, in closed form.

# The closed form

``\ell_h`` is the *period* of the closed orbit of ``X_h`` on the contour: the flow has
``|\dot{x}| = |X_h| = |\nabla h|``, so ``ds = |\nabla h| \, d\xi`` and the manuscript's
``\xi`` is the flow's own time. The manuscript estimates it by integrating that flow
numerically; for this ``h`` it can be had exactly.

On the island centred at ``(\pi, \pi/2)``, write ``s = x_1 - \pi`` and ``t = x_2 - \pi/2``,
so that ``h = \cos^2 s \, \cos^2 t``. The contour ``h`` is ``\cos t = \sqrt{h}/\cos s``,
which exists for ``|s| \le s_{\max}``, ``\cos s_{\max} = \sqrt{h}``. Parameterising by ``s``
and eliminating ``t``,

```math
\ell_h = \frac{2}{\sqrt{h}} \int_0^{s_{\max}} \frac{ds}{\sqrt{\cos^2 s - h}} ,
```

and the substitution ``\sin s = \sqrt{1-h} \, \sin\theta`` maps this onto the complete
elliptic integral of the first kind with parameter ``m = 1 - h``:

```math
\ell_h = \frac{2}{\sqrt{h}} K(1-h) = \frac{\pi}{\sqrt{h} \, \mathrm{agm}(1, \sqrt{h})} .
```

The substitution is what removes the inverse-square-root singularity at the turning point,
which is also why [`contour_length_quadrature`](@ref) — the independent check — uses it
rather than integrating in ``s`` directly.

# The two limits

At the island centre ``h \to 1`` this gives ``\ell_h \to \pi``, which is the period of the
harmonic approximation ``h \approx 1 - s^2 - t^2``: that flow rotates at angular frequency
``2``. At the separatrix ``h \to 0`` both factors diverge and ``\ell_h \to \infty``,
consistently with ``\nabla h = 0`` there.
"""
function contour_length(h::T) where {T <: AbstractFloat}
    zero(T) < h <= one(T) ||
        throw(DomainError(h, "a contour inside an island has 0 < h <= 1"))
    return π / (sqrt(h) * agm(one(T), sqrt(h)))
end

contour_length(h) = contour_length(float(h))

@doc raw"""
    relaxation_time(h)

The relaxation time ``\tau_h = (\ell_h / 2\pi)^2`` of `eq:relaxation-time`, on a contour of
[`islands_h`](@ref).

With the closed form of [`contour_length`](@ref) this is
``\tau_h = 1 / \big(4 h \, \mathrm{agm}(1,\sqrt{h})^2\big)``, which is ``1/4`` at the island
centre and diverges logarithmically at the separatrix.
"""
relaxation_time(h) = (contour_length(h) / 2π)^2

@doc raw"""
    contour_length_quadrature(h; n = 200)

``\ell_h`` computed by quadrature of ``\oint ds / |\nabla h|`` along the contour, as an
independent check on the closed form of [`contour_length`](@ref).

This integrates the *arclength* form that `eq:relaxation-time` actually writes — it forms
``|\nabla h|`` from the derivatives of [`islands_h`](@ref) and accumulates ``ds`` from the
parameterised curve — so agreement with [`contour_length`](@ref) tests the closed-form
derivation rather than restating it. `n` is the number of Gauss-Legendre nodes per quarter
contour.
"""
function contour_length_quadrature(h; n::Int = 200)
    ξ, wq = _gauss_legendre(n, -π / 2, π / 2)
    total = 0.0
    for (θ, w) in zip(ξ, wq)
        cs, sn, ct, st = _contour_point(h, θ)
        # ds/dθ from the parameterisation, and dt/dθ from differentiating cos t = √h/cos s.
        # Both are bounded at the turning points θ = ±π/2, where ds/dθ vanishes at exactly
        # the rate that keeps dt/dθ finite.
        dsdθ = sqrt(1 - h) * cos(θ) / cs
        dtdθ = sqrt(h) * sqrt(1 - h) * sin(θ) * sign(cos(θ)) / cs^2
        arc = sqrt(dsdθ^2 + dtdθ^2)
        gradnorm = 2 * sqrt((sn * cs * ct^2)^2 + (cs^2 * ct * st)^2)
        total += w * arc / gradnorm
    end
    # The upper branch is half the closed contour; the lower is its mirror image.
    return 2 * total
end

"""
    _contour_point(h, θ)

`(cos s, sin s, cos t, sin t)` at the point of the contour `h` parameterised by `θ` through
`sin s = √(1-h) sin θ`, in the island-centred coordinates of [`contour_length`](@ref).

# Why these four, and not `cos`/`acos` of the obvious expressions

Written the direct way, ``\\sin t = \\sqrt{1 - h/\\cos^2 s}`` cancels catastrophically at the
turning points, where ``\\cos^2 s \\to h`` makes the bracket a difference of two nearly equal
numbers. Near the separatrix that costs most of the mantissa: at ``h = 10^{-4}`` the relative
error in [`contour_length_quadrature`](@ref) was ``1.3 \\times 10^{-10}`` and grew with the
node count instead of falling, which is the signature of round-off rather than of truncation.

Both identities below are exact and cancellation-free — every term is positive:

```math
\\cos^2 s = h + (1-h)\\cos^2\\theta , \\qquad
\\sin t = \\frac{\\sqrt{1-h} \\, |\\cos\\theta|}{\\cos s} .
```

The first follows from ``\\cos^2 s = 1 - (1-h)\\sin^2\\theta``, the second from
``\\sin^2 t = 1 - h/\\cos^2 s = (\\cos^2 s - h)/\\cos^2 s`` with
``\\cos^2 s - h = (1-h)\\cos^2\\theta``.
"""
function _contour_point(h, θ)
    cs = sqrt(h + (1 - h) * cos(θ)^2)
    sn = sqrt(1 - h) * sin(θ)
    ct = sqrt(h) / cs
    st = sqrt(1 - h) * abs(cos(θ)) / cs
    return cs, sn, ct, st
end

@doc raw"""
    contour_average(u₀, h, centre; n = 400)

The ``t \to \infty`` limit of parallel diffusion on the contour ``h`` of the island at
`centre`,

```math
u_\infty(h) = \frac{1}{\ell_h} \oint_{C_h} \frac{u_0 \, ds}{|\nabla h|} ,
```

the average of the initial condition over the contour, weighted by the flow's own time.

This is the reference the A1 check compares the relaxed state against. It is quadrature in
the ``\theta`` of [`contour_length`](@ref), where the weight ``d\xi = d\theta /
(2\sqrt{h}\cos s)`` is smooth and bounded — integrating in arclength instead would put an
inverse-square-root singularity at each of the two turning points.

Both branches ``t > 0`` and ``t < 0`` are averaged, which is what closes the contour.
"""
function contour_average(u₀, h, centre; n::Int = 400)
    pts, wts = contour_samples(h, centre; n = n)
    return sum(w * u₀(p[1], p[2]) for (p, w) in zip(pts, wts)) / sum(wts)
end

@doc raw"""
    contour_samples(h, centre; n = 400)

Points on the contour `h` of the island at `centre`, and the weights ``d\xi`` of the flow-time
measure, as `(points, weights)`.

Both branches ``t > 0`` and ``t < 0`` are returned, which is what closes the contour, so there
are `2n` points. `sum(weights)` is [`contour_length`](@ref)`(h)`, and
``\sum w_i f(x_i) / \sum w_i`` is [`contour_average`](@ref) — this is the rule that function
integrates with, exposed so that a *run* can be sampled on the same points.

The weight is ``d\xi = d\theta / (2\sqrt{h}\cos s)``, smooth and bounded, rather than the
arclength ``ds/|\nabla h|``, which has an inverse-square-root singularity at each turning
point.
"""
function contour_samples(h, centre; n::Int = 400)
    ξ, wq = _gauss_legendre(n, -π / 2, π / 2)
    pts = NTuple{2, Float64}[]
    wts = Float64[]
    for (θ, w) in zip(ξ, wq)
        cs, sn, ct, st = _contour_point(h, θ)
        # `atan(st, ct)` rather than `acos(ct)`: the latter loses half the mantissa where
        # ct → 1, which is exactly the turning point of a near-separatrix contour.
        t = atan(st, ct)
        s = asin(sn)
        dξ = w / (2 * sqrt(h) * cs)
        x₁ = centre[1] + s
        push!(pts, (x₁, centre[2] + t), (x₁, centre[2] - t))
        push!(wts, dξ, dξ)
    end
    return pts, wts
end

@doc raw"""
    contour_deviation(u, h, centre; n = 400)

The flow-time-weighted standard deviation of `u` along the contour `h`,
``\big(\overline{u^2} - \bar{u}^2\big)^{1/2}``, with the average of
[`contour_average`](@ref).

This is what decays at the rate `eq:relaxation-time` predicts. The manuscript's own solution
of the reduced equation is
``v(t,\theta,h) = \sum_n \hat{v}_n(0) e^{-n^2 \kappa_h t + in\theta}`` with
``\kappa_h = 1/\tau_h``: the ``n = 0`` mode is the contour average and is constant, and
everything else decays, the slowest at ``e^{-t/\tau_h}``. So the deviation from the average
decays at rate ``1/\tau_h`` asymptotically, and measuring that rate against the closed form is
a direct test of `eq:relaxation-time` on the run rather than a restatement of it.
"""
function contour_deviation(u, h, centre; n::Int = 400)
    pts, wts = contour_samples(h, centre; n = n)
    W = sum(wts)
    vals = [u(p[1], p[2]) for p in pts]
    m = sum(w * v for (w, v) in zip(wts, vals)) / W
    return sqrt(max(sum(w * (v - m)^2 for (w, v) in zip(wts, vals)) / W, 0.0))
end

"""
    _gauss_legendre(n, a, b)

`n`-point Gauss-Legendre nodes and weights on `[a, b]`, by Newton iteration on the Legendre
polynomial.

Written out rather than taken from `QuadratureRules` because the contour integrals above are
the only quadrature in this file that is not already supplied by a discretisation, and the
rule is needed at a resolution (`n` in the hundreds) chosen for the integrand rather than
for a basis.
"""
function _gauss_legendre(n::Int, a::Real, b::Real)
    x = zeros(n)
    w = zeros(n)
    for i in 1:n
        # Chebyshev starting guess, then Newton on P_n via the standard recurrence.
        z = cos(π * (i - 0.25) / (n + 0.5))
        for _ in 1:100
            p0, p1 = 1.0, 0.0
            for j in 1:n
                p0, p1 = ((2j - 1) * z * p0 - (j - 1) * p1) / j, p0
            end
            dp = n * (z * p0 - p1) / (z^2 - 1)
            dz = p0 / dp
            z -= dz
            abs(dz) < 1e-15 && break
        end
        p0, p1 = 1.0, 0.0
        for j in 1:n
            p0, p1 = ((2j - 1) * z * p0 - (j - 1) * p1) / j, p0
        end
        dp = n * (z * p0 - p1) / (z^2 - 1)
        x[i] = (a + b) / 2 + (b - a) / 2 * z
        w[i] = (b - a) / (1 - z^2) / dp^2
    end
    return x, w
end

## Closed-form entropy minimisers

@doc raw"""
    analytic_minimiser(h, M₀, H₀, hΩ, h_norm²)

The unique constrained entropy minimiser `eq:u-eta_analytic` of the analytic test case,

```math
u_\eta = \frac{M_0}{4\pi^2} +
    \frac{H_0}{\|h - h_\Omega\|^2_{L^2}} (h - h_\Omega) ,
```

as a function of position, given the mass `M₀` and energy `H₀` of the initial condition, the
mean `hΩ` of ``h``, and ``\|h - h_\Omega\|^2_{L^2}`` as `h_norm²`.

The three scalars are passed in rather than computed here because they are quadratures, and
which quadrature is the discretisation's business: the spectral and the spline solver
compute them in their own inner products, and the whole point of A1 is that the two then
agree.
"""
function analytic_minimiser(h, M₀, H₀, hΩ, h_norm²)
    (x₁, x₂) -> M₀ / DOMAIN_AREA + H₀ / h_norm² * (h(x₁, x₂) - hΩ)
end

@doc raw"""
    analytic_entropy_minimum(H₀, h_norm²)

``S_\eta = (H_0^2/2) / \|h - h_\Omega\|^2_{L^2}``, `eq:Seta-analytic`, the value of the
entropy at [`analytic_minimiser`](@ref).
"""
analytic_entropy_minimum(H₀, h_norm²) = H₀^2 / 2 / h_norm²

@doc raw"""
    euler_minimiser(H₀, θ₀, θ₁, θ₂)

One member of the three-parameter family `eq:u-eta_Euler_periodic` of constrained entropy
minimisers of the reduced Euler test case,

```math
\omega(x) = \frac{\sqrt{H_0}}{\pi} \big[
    \cos\theta_0 \cos(x_1 + \theta_1) + \sin\theta_0 \cos(x_2 + \theta_2) \big] ,
```

as a function of position.

The minimum is attained on a whole family rather than at a point, because the lowest
nontrivial eigenvalue ``\lambda_2 = 1`` of ``-\Delta`` on ``\mathbb{T}^2`` is fourfold
degenerate. That is why comparing a relaxed state against "the" minimiser requires a fit
over the three phases — see [`best_fit_euler`](@ref) — rather than a subtraction.
"""
function euler_minimiser(H₀, θ₀, θ₁, θ₂)
    a = sqrt(H₀) / π
    (x₁, x₂) -> a * (cos(θ₀) * cos(x₁ + θ₁) + sin(θ₀) * cos(x₂ + θ₂))
end

@doc raw"""
    euler_entropy_minimum(H₀)

``S_\eta = H_0``, `eq:Seta-Euler_periodic`: the constrained entropy minimum of the reduced
Euler test case is the initial energy itself, since ``\lambda_2 = 1``.

This is the horizontal line every entropy trace in Section 4 is measured against. Under the
metric double bracket the entropy plateaus strictly *above* it — the manuscript's incomplete
relaxation — and under the projector bracket it reaches it.
"""
euler_entropy_minimum(H₀) = H₀

## The four runs of Section 4

@doc raw"""
    RunSpec

The specification of one of the four Section 4 runs: the bracket to use, the Gaussian and the
`background` field added to it, the prescribed `h` where there is one, and the time step and
final time.

The initial condition is held as its two pieces rather than as one closure, so that
[`periodise`](@ref) can reach the [`Gaussian`](@ref) inside it — A4's ``u_0 = \\cos(2x_2) +
u_G`` is the one case where the two readings of ``u_G`` differ materially. Use
[`initial_condition`](@ref) to get the callable.

`Δt` is the manuscript's own; `T` is not, and is recorded here as a choice of this
reproduction. See [`SECTION4_RUNS`](@ref).
"""
struct RunSpec{B, H}
    name::String
    section::String
    bracket::Symbol
    gaussian::Gaussian{Float64}
    background::B
    h::H
    Δt::Float64
    T::Float64
end

"""
    initial_condition(spec)

The initial condition of run `spec` as a callable `(x₁, x₂)`: the Gaussian, plus the
`background` field where the run has one.
"""
function initial_condition(spec::RunSpec)
    spec.background === nothing ? spec.gaussian :
    (x₁, x₂) -> spec.background(x₁, x₂) + spec.gaussian(x₁, x₂)
end

function periodise(spec::RunSpec, images::Int = 2)
    RunSpec(spec.name * "-periodic", spec.section, spec.bracket,
        periodise(spec.gaussian, images), spec.background, spec.h, spec.Δt, spec.T)
end

@doc raw"""
The four runs of Section 4, keyed by name.

| run | § | bracket | initial condition | ``\Delta t`` |
|:--|:--|:--|:--|:--|
| `a1` | 4.1 | double | ``u_G``, ``x_0 = (\pi, \pi{+}0.1)``, ``w = (0.25, 0.4)``, ``N = 2\pi w_1 w_2`` | ``10^{-4}`` |
| `a2` | 4.1 | double | ``u_G``, ``x_0 = (\pi,\pi)``, ``w = (0.3, 1)``, ``N = 1`` | ``10^{-3}`` |
| `a3` | 4.2 | projector | as `a2` | ``10^{-3}`` |
| `a4` | 4.2 | projector | ``\cos(2x_2) + u_G``, ``1/N = 1.8``, ``x_0 = (\pi, 3\pi/2)`` | ``10^{-3}`` |

Everything in this table is the manuscript's, **except the final time ``T``**, which it never
states.

**A1 takes ``T = 10``.** It relaxes on the scale of [`relaxation_time`](@ref), which is
``1/4`` at an island centre and diverges at the separatrix, so this is 40 relaxation times at
the centre: enough to resolve the island interiors while leaving the separatrix stalled, which
is the behaviour Fig. 1 shows. A1 still costs the most of the four, because the manuscript's
``\Delta t = 10^{-4}`` is set by the explicit stability limit rather than by accuracy and buys
it ten times the step count of the others.

**A2, A3 and A4 take ``T = 20``, and the reason is the rates rather than the figures.** Fig. 6
puts the vertex of the cone at ``t \approx 5``, so ``T = 10`` would be ample to *draw* every
panel of §4.2 — but not to *measure* its exponents. The linearised spectrum of the projector
flow is exactly ``\{1 - 1/\lambda\}`` over the eigenvalues of ``-\Delta``
(`scripts/verify_projector_rates.jl`), so the manuscript's rate ``\approx 1/2`` is the
``\lambda = 2`` mode exactly, and the next mode up, ``\lambda = 4``, decays at ``3/4`` — only
``1/4`` faster. It therefore contaminates a fitted rate as ``e^{-t/4}``, which is still
**8.2 %** of the signal at ``t = 10`` and biases the fit to ``0.519``. At ``t = 20`` it is
**0.67 %** and the bias is ``0.0017``. A ``T = 10`` run measured ``0.597`` for what is exactly
``1/2``, which is a pre-asymptotic window and not a defect.
"""
const SECTION4_RUNS = Dict(
    "a1" => RunSpec("a1", "4.1", :double,
        Gaussian((π, π + 0.1), (0.25, 0.4), 2π * 0.25 * 0.4), nothing,
        islands_h, 1e-4, 10.0),
    "a2" => RunSpec("a2", "4.1", :double,
        Gaussian((π, π), (0.3, 1.0), 1.0), nothing, nothing, 1e-3, 20.0),
    "a3" => RunSpec("a3", "4.2", :projector,
        Gaussian((π, π), (0.3, 1.0), 1.0), nothing, nothing, 1e-3, 20.0),
    "a4" => RunSpec("a4", "4.2", :projector,
        Gaussian((π, 3π / 2), (0.3, 1.0), 1 / 1.8), (x₁, x₂) -> cos(2x₂),
        nothing, 1e-3, 20.0))

"The four runs in the order Section 4 presents them."
const SECTION4_ORDER = ("a1", "a2", "a3", "a4")
