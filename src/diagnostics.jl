#
# The diagnostics every Section 4 run reports, written once against both discretisations.
#
# `Diagnostics` is the one place that knows which energy a run has -- linear for the analytic
# test case, the elliptic ½(φ,ω) for the reduced Euler ones -- so the drivers and the figures
# never branch on it. Everything below dispatches on the solver, so the spectral and the spline
# run produce numbers that are comparable by construction rather than by convention.
#

@doc raw"""
    Diagnostics(solver, spec)

The conserved and monitored quantities of run `spec` on `solver`, which may be a
[`SpectralTorus`](@ref) or a [`SplineTorus`](@ref).

Holds whatever the run's energy needs precomputed: ``h - h_\Omega`` for the analytic test
case, nothing for the reduced Euler ones, where the generating field is recomputed from the
state.
"""
struct Diagnostics{S, X}
    solver::S
    spec::RunSpec
    hz::X
end

function Diagnostics(solver, spec::RunSpec)
    hz = spec.h === nothing ? nothing : _centred_h(solver, spec.h)
    Diagnostics{typeof(solver), typeof(hz)}(solver, spec, hz)
end

function _centred_h(g::SpectralTorus, h)
    H = torus_field(g, h)
    return H .- mean_value(g, H)
end

function _centred_h(t::SplineTorus, h)
    ĥ = project(t.space, x -> h(x[1], x[2]))
    return ĥ .- mean_value(t, ĥ)
end

@doc raw"""
    potential(d, ω)

The stream function ``\phi`` of `eq:Poisson-eq-periodic`, or ``h - h_\Omega`` for the
analytic test case — in both cases ``\delta H / \delta u``, which is what the metric bracket
is degenerate on and what the scatter plots are drawn against.
"""
potential(d::Diagnostics{<:SpectralTorus}, ω) = d.hz === nothing ?
                                                poisson_periodic(d.solver, ω) : d.hz

potential(d::Diagnostics{<:SplineTorus}, ω̂) = d.hz === nothing ? d.solver.Λ * ω̂ : d.hz

@doc raw"""
    energy(d, ω)

``H``: the linear ``(h - h_\Omega, u)_{L^2}`` of `eq:analytical_H` for the analytic test
case, and the quadratic ``\tfrac12 (\phi, u)_{L^2}`` of `eq:Euler_H_periodic` otherwise.

Both are conserved exactly by the semi-discrete flow, because the metric bracket is
degenerate on them; what a run measures is how much of that survives the time discretisation.
"""
function energy(d::Diagnostics, ω)
    d.hz === nothing ? l2inner(d.solver, potential(d, ω), ω) / 2 :
    l2inner(d.solver, d.hz, ω)
end

@doc raw"""
    entropy(d, ω)

``S = \tfrac12 \int_\Omega \omega^2 dx``, the entropy `Som2` with ``s(y) = y^2/2``.
"""
entropy(d::Diagnostics, ω) = l2inner(d.solver, ω, ω) / 2

"``\\int_\\Omega \\omega \\, dx``, which is zero for the whole run and is monitored as a check."
vorticity_mass(d::Diagnostics, ω) = integrate(d.solver, ω)

"``\\|\\phi\\|^2_{L^2}``, the second coordinate of the cone diagram of Fig. 6."
function potential_norm²(d::Diagnostics, ω)
    φ = potential(d, ω)
    return l2inner(d.solver, φ, φ)
end

@doc raw"""
    Trace

The time series a run records: the time, energy, entropy, ``\|\phi\|^2`` and vorticity mass
at each sample, plus the state at the first and last sample.

Samples are taken every `stride` steps rather than every step — a Section 4 run takes 10⁴ to
10⁵ of them and nothing in the manuscript's figures needs that resolution.
"""
mutable struct Trace{X}
    t::Vector{Float64}
    H::Vector{Float64}
    S::Vector{Float64}
    φ²::Vector{Float64}
    M::Vector{Float64}
    initial::X
    final::X
end

function Trace(initial)
    Trace(Float64[], Float64[], Float64[], Float64[], Float64[],
        copy(initial), copy(initial))
end

"""
    record!(tr, d, t, ω)

Append one sample to `tr`, and keep `ω` as the running final state.
"""
function record!(tr::Trace, d::Diagnostics, t, ω)
    push!(tr.t, t)
    push!(tr.H, energy(d, ω))
    push!(tr.S, entropy(d, ω))
    push!(tr.φ², potential_norm²(d, ω))
    push!(tr.M, vorticity_mass(d, ω))
    tr.final = copy(ω)
    return tr
end

@doc raw"""
    energy_error(tr)

The relative energy error ``|H(t) - H_0| / |H_0|`` along the trace, which is the quantity the
lower panels of Figs. 2 and 4 plot.
"""
energy_error(tr::Trace) = abs.(tr.H .- tr.H[1]) ./ abs(tr.H[1])

@doc raw"""
    entropy_monotone(tr; tol = 0)

Whether the entropy never increases along the trace, and by how much it does where it does.

Returns `(ok, worst)` with `worst` the largest positive increment, normalised by ``S_0``. The
manuscript's claim is monotone dissipation, and asserting only that ``S`` *ends* lower would
pass for a trace that rose and fell.
"""
function entropy_monotone(tr::Trace; tol = 0.0)
    Δ = diff(tr.S) ./ abs(tr.S[1])
    worst = isempty(Δ) ? 0.0 : maximum(Δ)
    return (worst <= tol, worst)
end

## The reduced Euler entropy minimiser, fitted

@doc raw"""
    best_fit_euler(d, ω, H₀)

The member of the three-phase family `eq:u-eta_Euler_periodic` closest to `ω` in ``L^2``,
returned as `(fit, residual, coefficients)`.

# It is a projection, not a search

The manuscript computes this "varying the three phases ``\theta_0``, ``\theta_1``,
``\theta_2``". Expanding the family,

```math
\omega_\eta = A \big( a_1 \cos x_1 + b_1 \sin x_1 + a_2 \cos x_2 + b_2 \sin x_2 \big) ,
\qquad A = \frac{\sqrt{H_0}}{\pi} ,
```

with ``a_1 = \cos\theta_0\cos\theta_1``, ``b_1 = -\cos\theta_0\sin\theta_1`` and likewise for
the second pair, so that ``a_1^2 + b_1^2 + a_2^2 + b_2^2 = 1``. The three phases parameterise
exactly the unit sphere in those four coefficients, and the four modes are ``L^2``-orthogonal
on ``\mathbb{T}^2`` with ``\|\cos x_1\|^2 = 2\pi^2``. So

```math
\|\omega - \omega_\eta\|^2 = \|\omega\|^2 - 2A \, v \cdot p + 2\pi^2 A^2 ,
\qquad p_i = (\omega, e_i) ,
```

which over ``\|v\| = 1`` is minimised at ``v = p / \|p\|`` — the normalised projection. No
optimisation is needed and none is done; `verify_diagnostics.jl` checks the closed form
against a brute-force scan over the three phases.

The amplitude ``A`` is **fixed** by the energy and is not fitted. That is the manuscript's own
statement — the family is parameterised by phases alone — and it is what makes the residual a
measure of distance from the constraint set ``\mathfrak{C}_\eta`` rather than a plain
four-mode truncation error.
"""
function best_fit_euler(d::Diagnostics, ω, H₀)
    A = sqrt(H₀) / π
    modes = _euler_modes(d.solver)
    p = [l2inner(d.solver, ω, e) for e in modes]
    n = sqrt(sum(abs2, p))
    v = n > 0 ? p ./ n : [1.0, 0.0, 0.0, 0.0]
    fit = sum(A * v[i] .* modes[i] for i in 1:4)
    res = sqrt(max(l2inner(d.solver, ω .- fit, ω .- fit), 0.0))
    return (fit, res, A .* v)
end

"The four lowest nontrivial eigenmodes of ``-\\Delta`` on ``\\mathbb{T}^2``, in the solver's
own representation."
function _euler_modes(g::SpectralTorus)
    (torus_field(g, (a, b) -> cos(a)), torus_field(g, (a, b) -> sin(a)),
        torus_field(g, (a, b) -> cos(b)), torus_field(g, (a, b) -> sin(b)))
end

function _euler_modes(t::SplineTorus)
    (project(t.space, x -> cos(x[1])), project(t.space, x -> sin(x[1])),
        project(t.space, x -> cos(x[2])), project(t.space, x -> sin(x[2])))
end

## The cone of eq:theoretical-limits

@doc raw"""
    cone_coordinates(tr)

The trajectory in the plane ``(S(u), 1/\|\phi\|^2_{L^2})`` of Fig. 6.
"""
cone_coordinates(tr::Trace) = (tr.S, 1 ./ tr.φ²)

@doc raw"""
    cone_residual(tr, H₀)

How far the trajectory violates `eq:theoretical-limits`, as
`(below_Sη, below_lower, above_upper)` — the worst violation of each of the three
inequalities, normalised, and negative where the inequality holds with room to spare.

The cone is

```math
S(u) \ge H_0 = S_\eta , \qquad
\frac{1}{2H_0} \le \frac{1}{\|\phi\|^2} \le
\frac{1}{2H_0} + \frac{S(u) - S_\eta}{2H_0^2} .
```

The lower bound is the Poincaré inequality `eq:tb3` and the upper one `eq:tb1`, which is
Cauchy-Schwarz. Both are statements about *any* state of energy ``H_0``, not about the
dynamics, so a violation is a defect in the run rather than a falsification of the
manuscript — which is why they are checked as a *precondition* of the rate fits.
"""
function cone_residual(tr::Trace, H₀)
    S, inv² = tr.S, 1 ./ tr.φ²
    lower = 1 / (2H₀)
    upper = lower .+ (S .- H₀) ./ (2H₀^2)
    return (maximum((H₀ .- S) ./ abs(H₀)),
        maximum((lower .- inv²) .* (2H₀)),
        maximum((inv² .- upper) .* (2H₀)))
end

## Rate fitting

@doc raw"""
    fit_rate(t, y; window = (0.5, 0.9), floor = 1e-12)

The exponential decay rate of `y` against `t`, as `(rate, r², n)`.

A least-squares straight line through ``\log y`` over the sub-interval
`window .* (t[end] - t[1]) .+ t[1]`, returning the negated slope. `r²` is the coefficient of
determination of that fit and is reported alongside, because a rate quoted without it says
nothing about whether the decay was exponential at all.

Samples with `y <= floor * maximum(y)` are dropped: once the quantity reaches round-off its
logarithm is noise, and including it drags the slope toward zero. Which samples were used is
returned as `n`.

# Why a window and not the whole trace

The manuscript's own reading of Fig. 6 is that the rate is asymptotic — for the second initial
condition "the initial entropy relaxation rate is slower, but approaches ≈ 1 as the trajectory
approaches the vertex of the cone". Fitting from ``t = 0`` would therefore measure a mixture
of the transient and the asymptote and report neither. The window is a recorded choice.
"""
function fit_rate(t::AbstractVector, y::AbstractVector; window = (0.5, 0.9), floor = 1e-12)
    t₀, t₁ = t[1] + window[1] * (t[end] - t[1]), t[1] + window[2] * (t[end] - t[1])
    ymax = maximum(y)
    idx = [i
           for i in eachindex(t)
           if t[i] >= t₀ && t[i] <= t₁ && y[i] > floor * ymax && isfinite(y[i])]
    length(idx) < 3 && return (NaN, NaN, length(idx))

    ts = t[idx]
    ls = log.(y[idx])
    t̄, l̄ = sum(ts) / length(ts), sum(ls) / length(ls)
    Stt = sum((ts .- t̄) .^ 2)
    slope = sum((ts .- t̄) .* (ls .- l̄)) / Stt
    intercept = l̄ - slope * t̄
    ss_res = sum((ls .- (intercept .+ slope .* ts)) .^ 2)
    ss_tot = sum((ls .- l̄) .^ 2)
    r² = ss_tot > 0 ? 1 - ss_res / ss_tot : NaN
    return (-slope, r², length(idx))
end

## Scatter data

@doc raw"""
    scatter_data(d, ω, N)

The points ``(\phi_{ij}, \omega_{ij})`` of the scatter plots of Figs. 1, 3 and 5, sampled on
an `N`-by-`N` grid.

For the analytic test case ``\phi`` is ``h`` itself, and the plot is the manuscript's
``h``-``u`` plane; for the reduced Euler runs it is the stream function. A functional relation
appearing in this cloud is the manuscript's evidence that the final state is an equilibrium.
"""
function scatter_data(d::Diagnostics{<:SpectralTorus}, ω, N::Int = 0)
    g = d.solver
    φ = d.hz === nothing ? poisson_periodic(g, ω) : torus_field(g, d.spec.h)
    return (vec(φ), vec(ω))
end

function scatter_data(d::Diagnostics{<:SplineTorus}, ω̂, N::Int = 128)
    t = d.solver
    φ̂ = d.hz === nothing ? t.Λ * ω̂ : project(t.space, x -> d.spec.h(x[1], x[2]))
    return (vec(spline_grid(t, φ̂, N)), vec(spline_grid(t, ω̂, N)))
end
