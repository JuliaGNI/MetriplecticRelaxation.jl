#
# The Fourier spectral solver: the manuscript's own discretisation for Section 4.
#
# This is the reference half of the pair. The manuscript solves every Section 4 run "with a
# standard spectral method with Fourier basis, on a 256 x 256 uniform grid" (line 3874), so this
# reproduces the published method, and `spline.jl` reproduces the same *equations* in a B-spline
# Galerkin discretisation that shares no code path with it. Agreement between the two is the
# control: two unrelated discretisations agreeing is evidence about the equation rather than
# about either method.
#
# The machinery is the same as `Papers/Metriplectic Relaxation to Equilibria/scripts/torus.jl`,
# which was written for the verification pass of that manuscript -- `poisson_periodic`,
# `canonical_bracket`, `hamiltonian_field` and the trapezoidal `integrate`. It could not be
# depended on as a module: it lives in a manuscript repository with no package boundary and no
# UUID, so `using` it here would mean either a hard-coded absolute path in committed code or a
# copy of the file. What is done instead is to build the same operators on PoissonBrackets'
# `TorusGrid`, and to CHECK the two agree -- `verify_spectral.jl` compares every derivative
# against `spectral_grid`'s differentiation matrix, which is the package's independent
# implementation of the same operator.
#

@doc raw"""
    SpectralTorus(N)

An `N`-by-`N` Fourier collocation grid on ``\Omega = [0, 2\pi]^2``, with the wavenumbers and
FFT plans the spectral operators need.

Fields are `N`-by-`N` matrices indexed `[i,j]` with `i` running over ``x_1``, matching
PoissonBrackets' [`TorusGrid`](@ref) convention so that the two can be compared elementwise.

# The Nyquist mode

Wavenumbers are wrapped as `0:N÷2` then `-N÷2+1:-1`, so for even `N` the mode `N÷2` is
assigned ``+N/2``. For a real field the coefficient there is real, so its contribution to a
first derivative is purely imaginary and is discarded by taking the real part — which is what
leaves the differentiation antisymmetric. Nothing in Section 4 is band-limited exactly at
Nyquist in any case: the fields are Gaussians and low trigonometric modes.
"""
struct SpectralTorus{T, P, IP}
    N::Int
    x₁::Matrix{T}
    x₂::Matrix{T}
    k₁::Matrix{T}
    k₂::Matrix{T}
    Δ⁻¹::Matrix{T}
    plan::P
    iplan::IP
end

function SpectralTorus(N::Int)
    N >= 4 ||
        throw(ArgumentError("a spectral grid needs at least four points, got N = $(N)"))
    xs = collect(0:(N - 1)) .* (DOMAIN_LENGTH / N)
    ks = float.(vcat(collect(0:(N ÷ 2)), collect((-N ÷ 2 + 1):-1)))
    x₁ = [xs[i] for i in 1:N, j in 1:N]
    x₂ = [xs[j] for i in 1:N, j in 1:N]
    k₁ = [ks[i] for i in 1:N, j in 1:N]
    k₂ = [ks[j] for i in 1:N, j in 1:N]
    # The inverse of -Δ, with the zero mode annihilated rather than treated as an error: on the
    # torus the Poisson problem is solvable only for mean-zero data, and setting the mean of φ
    # to zero is `eq:Poisson-eq-periodic`'s own normalisation.
    Δ⁻¹ = [(k₁[i, j]^2 + k₂[i, j]^2) == 0 ? 0.0 : 1 / (k₁[i, j]^2 + k₂[i, j]^2)
           for i in 1:N, j in 1:N]
    buf = zeros(ComplexF64, N, N)
    plan = plan_fft(buf)
    iplan = plan_ifft(buf)
    SpectralTorus(N, x₁, x₂, k₁, k₂, Δ⁻¹, plan, iplan)
end

Base.size(g::SpectralTorus) = (g.N, g.N)
Base.length(g::SpectralTorus) = g.N^2

function Base.show(io::IO, g::SpectralTorus)
    print(io, "SpectralTorus(N=", g.N, ")")
end

"""
    torus_field(g, f)

Evaluate `f(x₁, x₂)` at every node of `g`.
"""
torus_field(g::SpectralTorus, f) = f.(g.x₁, g.x₂)

_apply(g::SpectralTorus, u, mult) = real.(g.iplan * (mult .* (g.plan * complex.(u))))

"Partial derivative in the first coordinate."
∂₁(g::SpectralTorus, u) = _apply(g, u, im .* g.k₁)

"Partial derivative in the second coordinate."
∂₂(g::SpectralTorus, u) = _apply(g, u, im .* g.k₂)

"The Laplacian ``\\Delta u``."
laplacian(g::SpectralTorus, u) = _apply(g, u, .-(g.k₁ .^ 2 .+ g.k₂ .^ 2))

@doc raw"""
    poisson_periodic(g, ω)

Solve ``-\Delta \phi = \omega`` on the torus with ``\phi_\Omega = 0``, i.e.
`eq:Poisson-eq-periodic`.

The zero mode of `ω` is discarded rather than raising: the equation is solvable only for
mean-zero data, and the manuscript enforces that by evolving ``\omega = u - u_\Omega``
throughout. Passing a field with a nonzero mean therefore silently solves for its
fluctuation.
"""
poisson_periodic(g::SpectralTorus, ω) = _apply(g, ω, g.Δ⁻¹)

@doc raw"""
    canonical_bracket(g, f, k)

``[f, k] = \partial_1 f \, \partial_2 k - \partial_1 k \, \partial_2 f``, the canonical
Poisson bracket of the plane, defined after `eq:Euler-equilibrium`.
"""
function canonical_bracket(g::SpectralTorus, f, k)
    ∂₁(g, f) .* ∂₂(g, k) .- ∂₁(g, k) .* ∂₂(g, f)
end

@doc raw"""
    hamiltonian_field(g, h)

``X_h = (\partial_2 h, -\partial_1 h)``, so that ``X_h \cdot \nabla f = [f, h]`` and
``\nabla \cdot X_h = 0``.
"""
hamiltonian_field(g::SpectralTorus, h) = (∂₂(g, h), .-∂₁(g, h))

@doc raw"""
    integrate(g, u)

``\int_\Omega u \, dx`` by the rectangle rule, which on a periodic grid *is* the trapezoidal
rule and is spectrally accurate: every non-constant Fourier mode sums to zero over a full
period, so for a band-limited field it is exact.
"""
integrate(g::SpectralTorus, u) = sum(u) * (DOMAIN_LENGTH / g.N)^2

l2inner(g::SpectralTorus, u, v) = integrate(g, u .* v)
l2norm(g::SpectralTorus, u) = sqrt(max(l2inner(g, u, u), 0.0))
mean_value(g::SpectralTorus, u) = integrate(g, u) / DOMAIN_AREA

## The two vector fields of Section 4

@doc raw"""
    double_bracket_field(g, ω, h)

The right-hand side of `eq:parallel-diffusion`,

```math
\partial_t u = [h, [h, u]] = \nabla \cdot \big( X_h \otimes X_h \, \nabla u \big) ,
```

evaluated as the nested canonical bracket rather than as the divergence form.

# Why the nested bracket, and why they are the same

``[f,h] = \nabla f \cdot X_h``, so ``[h,[h,u]] = X_h \cdot \nabla (X_h \cdot \nabla u)``. The
divergence form expands to the same thing plus ``(\nabla \cdot X_h)(X_h \cdot \nabla u)``, and
``\nabla \cdot X_h = \partial_1 \partial_2 h - \partial_2 \partial_1 h`` vanishes identically.
So the two are equal in the continuum, and on a spectral grid they are equal to round-off,
because mixed partials of the DFT commute exactly.

The nested form is used because it costs four derivative applications per evaluation against
the divergence form's six, and because it is the form `eq:parallel-diffusion` is stated in.
`verify_spectral.jl` checks the two agree.
"""
function double_bracket_field(g::SpectralTorus, ω, h)
    canonical_bracket(g, h, canonical_bracket(g, h, ω))
end

@doc raw"""
    parallel_diffusion(g, ω, X)

The same right-hand side as [`double_bracket_field`](@ref), evaluated from an
**already-computed** Hamiltonian field `X` as ``X \cdot \nabla (X \cdot \nabla \omega)``.

Run A1's ``h`` is prescribed and fixed, so ``X_h`` is a constant of the whole simulation. The
nested-bracket form recomputes it anyway: each [`canonical_bracket`](@ref) differentiates both
of its arguments, so ``\partial_1 h`` and ``\partial_2 h`` are transformed four times per
evaluation out of eight derivative applications. Hoisting them halves the cost exactly, which
on A1 — the manuscript's ``\Delta t = 10^{-4}`` over ``10^5`` steps at ``256^2`` — is the
difference between 69 and 35 minutes.

`verify_spectral.jl` checks the two forms agree to round-off, so this is a hoist rather than a
second discretisation.
"""
function parallel_diffusion(g::SpectralTorus, ω, X)
    q = X[1] .* ∂₁(g, ω) .+ X[2] .* ∂₂(g, ω)
    return X[1] .* ∂₁(g, q) .+ X[2] .* ∂₂(g, q)
end

@doc raw"""
    projector_bracket_field(g, ω, φ)

The right-hand side of the projector-bracket evolution of §4.2 for the reduced Euler case,

```math
\partial_t u = - \Big[ \omega - \frac{2 H(u)}{\|\phi\|^2_{L^2}} \phi \Big] ,
\qquad H(u) = \tfrac{1}{2}(\phi, u)_{L^2} ,
```

with ``\phi`` the solution of `eq:Poisson-eq-periodic`.

!!! warning "The manuscript's prose drops the factor 2 here"
    The equation as printed just below `eq:projector-brackets` reads
    ``\partial_t u = -[u - u_\Omega - \fun{H}(u)\|\phi\|^{-2}\phi]``, with ``H`` rather than
    ``2H``. That is inconsistent with the manuscript's own `eq:SS-projector`. The coefficient
    is fixed by `eq:L2-projector`: ``c(u,v) = \|\phi\|^{-2}(\phi, v)`` with
    ``v = \delta S/\delta u = \omega`` gives ``c = \|\phi\|^{-2}(\phi,\omega) = 2H/\|\phi\|^2``,
    since `eq:Euler_H_periodic` defines ``H = \tfrac12(\phi,u)``. Substituting that same ``c``
    into ``(S,S) = (\omega, \Pi_H \omega)`` reproduces `eq:SS-projector`,
    ``(S,S) = 2S - 4H_0^2/\|\phi\|^2``, exactly — with ``H/\|\phi\|^2`` it would give
    ``2S - 2H_0^2/\|\phi\|^2`` and contradict the printed equation.

    The analytic test case a page earlier is *not* affected, because there
    `eq:analytical_H` defines ``H = (h - h_\Omega, u)`` with no factor ``\tfrac12``, so
    ``c = H/\|h-h_\Omega\|^2`` and the printed equation is right.

    `verify_projector_factor.jl` settles it numerically: only the factor 2 conserves the energy
    and reproduces `eq:SS-projector`.
"""
function projector_bracket_field(g::SpectralTorus, ω, φ)
    nφ² = l2inner(g, φ, φ)
    c = l2inner(g, φ, ω) / nφ²
    return .-(ω .- c .* φ)
end

@doc raw"""
    spectral_state(g, spec)

The initial state of run `spec` on the grid `g`, as the pair `(ω, u_Ω)`.

The state variable is the **vorticity** ``\omega = u - u_\Omega`` throughout, not ``u``. That
is the manuscript's own alternative — "equivalently, we could have chosen the phase space to
be the subspace of functions satisfying ``u_\Omega = 0`` and ``u = \omega``", stated at
`eq:Poisson-eq-periodic` — and it is what makes the projector bracket well posed: its
generating field ``\phi`` has zero mean by construction, so a state carrying a mean would
have that mean projected inconsistently.

Both vector fields preserve the zero mean exactly. The double bracket does because
``\nabla \cdot (X_h \otimes X_h \nabla u)`` is a divergence; the projector bracket because
both ``\omega`` and ``\phi`` are mean-free. ``u_\Omega`` is therefore a constant of the
motion, carried alongside and added back only for plotting.
"""
function spectral_state(g::SpectralTorus, spec::RunSpec)
    u = torus_field(g, initial_condition(spec))
    uΩ = mean_value(g, u)
    return (u .- uΩ, uΩ)
end

@doc raw"""
    spectral_step!(g, ω, spec, ĥ)

One classical fourth-order Runge-Kutta step of size `spec.Δt`, the manuscript's own
integrator ("the standard 4th order explicit Runge-Kutta method").

`ĥ` is the prescribed Hamiltonian sampled on the grid for the analytic test case, and
`nothing` for the reduced Euler runs, where the generating field is recomputed from
``\omega`` at every stage through `eq:Poisson-eq-periodic`.

Returns the new ``\omega``; the input is not modified.
"""
function spectral_step!(g::SpectralTorus, ω, spec::RunSpec, ĥ)
    f = _spectral_rhs(g, spec, ĥ)
    Δt = spec.Δt
    k1 = f(ω)
    k2 = f(ω .+ (Δt / 2) .* k1)
    k3 = f(ω .+ (Δt / 2) .* k2)
    k4 = f(ω .+ Δt .* k3)
    return ω .+ (Δt / 6) .* (k1 .+ 2 .* k2 .+ 2 .* k3 .+ k4)
end

"""
    _spectral_rhs(g, spec, ĥ)

The right-hand side of run `spec` as a closure of ``\\omega`` alone, which is what the
Runge-Kutta stages need.
"""
function _spectral_rhs(g::SpectralTorus, spec::RunSpec, ĥ)
    if spec.bracket === :double && ĥ !== nothing
        # A1: h is prescribed and fixed, so X_h is hoisted out of the time loop entirely.
        X = hamiltonian_field(g, ĥ)
        return ω -> parallel_diffusion(g, ω, X)
    elseif spec.bracket === :double
        # A2: h is the stream function φ, recomputed from ω at every stage.
        return ω -> double_bracket_field(g, ω, poisson_periodic(g, ω))
    elseif spec.bracket === :projector
        return ω -> projector_bracket_field(g, ω, poisson_periodic(g, ω))
    else
        throw(ArgumentError("unknown bracket $(spec.bracket)"))
    end
end
