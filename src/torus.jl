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

Callable as `g(x₁, x₂)`.

# It is not periodic, and that is faithful

The manuscript writes this formula unmodified on ``\mathbb{T}^2``, with no periodic
summation over images, so this reproduction does the same. The consequence is a jump across
the boundary of size ``u_G`` evaluated at distance ``\pi`` from the centre: negligible at
``10^{-69}`` for Section 4.1's ``w_2 = 0.4``, but ``5.2 \times 10^{-5}`` for the
``w_2 = 1`` of the reduced Euler runs. Both discretisations absorb it the same way — the
spectral one by sampling on a grid that is periodic by construction, the spline one by
``L^2`` projection onto a periodic basis — so it does not bias the comparison between them.
It does put a floor under how well either can represent the stated initial condition, which
is why the run reports use the *projected* initial condition rather than the formula when
they quote conserved quantities.
"""
struct Gaussian{T}
    x₀::NTuple{2, T}
    w::NTuple{2, T}
    N::T
end

# `x₀ = (π, π + 0.1)` is a `Tuple{Irrational, Float64}`, so the parameters have to be promoted
# rather than required to arrive already matching.
function Gaussian(x₀::Tuple, w::Tuple, N)
    T = promote_type(map(typeof ∘ float, (x₀..., w..., N))...)
    Gaussian{T}(map(T, x₀), map(T, w), T(N))
end

function (g::Gaussian)(x₁, x₂)
    exp(-(x₁ - g.x₀[1])^2 / g.w[1]^2 - (x₂ - g.x₀[2])^2 / g.w[2]^2) / g.N
end

(g::Gaussian)(x) = g(x[1], x[2])

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
    ξ, wq = _gauss_legendre(n, -π / 2, π / 2)
    num = 0.0
    den = 0.0
    for (θ, w) in zip(ξ, wq)
        cs, sn, ct, st = _contour_point(h, θ)
        # `atan(st, ct)` rather than `acos(ct)`: the latter loses half the mantissa where
        # ct → 1, which is exactly the turning point of a near-separatrix contour.
        t = atan(st, ct)
        s = asin(sn)
        dξ = w / (2 * sqrt(h) * cs)
        x₁ = centre[1] + s
        num += dξ * (u₀(x₁, centre[2] + t) + u₀(x₁, centre[2] - t))
        den += 2dξ
    end
    return num / den
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

The specification of one of the four Section 4 runs: its initial condition `u₀`, the bracket
to use, the prescribed `h` where there is one, and the time step and final time.

`Δt` is the manuscript's own; `T` is not, and is recorded here as a choice of this
reproduction. See [`SECTION4_RUNS`](@ref).
"""
struct RunSpec{F, H}
    name::String
    section::String
    bracket::Symbol
    u₀::F
    h::H
    Δt::Float64
    T::Float64
end

@doc raw"""
The four runs of Section 4, keyed by name.

| run | § | bracket | initial condition | ``\Delta t`` |
|:--|:--|:--|:--|:--|
| `a1` | 4.1 | double | ``u_G``, ``x_0 = (\pi, \pi{+}0.1)``, ``w = (0.25, 0.4)``, ``N = 2\pi w_1 w_2`` | ``10^{-4}`` |
| `a2` | 4.1 | double | ``u_G``, ``x_0 = (\pi,\pi)``, ``w = (0.3, 1)``, ``N = 1`` | ``10^{-3}`` |
| `a3` | 4.2 | projector | as `a2` | ``10^{-3}`` |
| `a4` | 4.2 | projector | ``\cos(2x_2) + u_G``, ``1/N = 1.8``, ``x_0 = (\pi, 3\pi/2)`` | ``10^{-3}`` |

Everything in this table is the manuscript's, **except the final time ``T``**, which it
never states. Fig. 6 places the vertex of the cone at ``t \approx 5``, so ``T = 10`` is
taken for the reduced Euler runs, leaving as much trajectory again past the vertex to
measure the asymptotic rates on. A1 has no such landmark and relaxes on the scale of
[`relaxation_time`](@ref), which diverges at the separatrix and is ``1/4`` at an island
centre; ``T = 20`` resolves the island interiors while leaving the separatrix stalled, which
is the behaviour Fig. 1 shows.
"""
const SECTION4_RUNS = Dict(
    "a1" => RunSpec("a1", "4.1", :double,
        Gaussian((π, π + 0.1), (0.25, 0.4), 2π * 0.25 * 0.4), islands_h, 1e-4, 20.0),
    "a2" => RunSpec("a2", "4.1", :double,
        Gaussian((π, π), (0.3, 1.0), 1.0), nothing, 1e-3, 10.0),
    "a3" => RunSpec("a3", "4.2", :projector,
        Gaussian((π, π), (0.3, 1.0), 1.0), nothing, 1e-3, 10.0),
    "a4" => RunSpec("a4", "4.2", :projector,
        let g = Gaussian((π, 3π / 2), (0.3, 1.0), 1 / 1.8)
            (x₁, x₂) -> cos(2x₂) + g(x₁, x₂)
        end, nothing, 1e-3, 10.0))
