#
# The B-spline Galerkin solver: the deliberate deviation from the manuscript.
#
# The manuscript uses a Fourier spectral method for Section 4 (`spectral.jl`). This half solves
# the same equations in a periodic tensor-product B-spline space, through PoissonBrackets'
# `TensorSplineSpace`, `DoubleBracket`, `ProjectorBracket` and `MetriplecticFlow`.
#
# THE NUMBERS WILL NOT MATCH BIT FOR BIT, and no write-up should suggest otherwise. What is
# reproduced is the *claims*: energy conserved to machine precision, monotone entropy decay,
# convergence to the stated closed-form references, the rates ~1 and ~1/2, and
# tau_h = (l_h/2pi)^2. Agreement between two discretisations that share no code path is evidence
# that those claims are about the equation rather than about either method; a disagreement is a
# finding.
#

@doc raw"""
    PoissonMap(space)

The solution operator of `eq:Poisson-eq-periodic`, ``\hat{\phi} = \Lambda \hat{\omega}`` with
``-\Delta \phi = \omega`` and ``\phi_\Omega = 0``, as a **lazy** matrix.

[`DoubleBracket`](@ref) and [`ProjectorBracket`](@ref) take their generating field either as a
coefficient vector or as the matrix of the linear map that produces it, and apply it as
`Λ * û`. That is the whole interface this type has to satisfy, so it can be a factorisation
rather than an assembled matrix.

# Why lazy, and why bordered

``\Lambda = \mathbb{K}^{-1} \mathbb{M}`` is **dense**: the inverse of the stiffness matrix has
no sparsity left. At the ``64^2`` cells of a Section 4 run that is a 4096-by-4096 matrix, 134
MB, and every Runge-Kutta stage would stream all of it — memory bandwidth, not arithmetic,
would set the run time, and a ``128^2`` run would need 2.1 GB and be impossible. Applying a
cached sparse factorisation instead is one triangular solve per stage.

``\mathbb{K}`` is singular on the torus: constants are in its kernel, which is
`eq:Poisson-eq-periodic`'s reason for fixing ``\phi_\Omega = 0``. Regularising it as
``\mathbb{K} + b b^T / |\Omega|`` would restore invertibility but destroy the sparsity, since
``b_i = \int_\Omega \Phi_i`` is dense. The **bordered** system keeps it:

```math
\begin{pmatrix} \mathbb{K} & b \\ b^T & 0 \end{pmatrix}
\begin{pmatrix} \hat{\phi} \\ \alpha \end{pmatrix} =
\begin{pmatrix} \mathbb{M} \hat{\omega} \\ 0 \end{pmatrix} ,
```

which adds one row and one column rather than ``N^2`` entries. It is nonsingular because
``b^T e = |\Omega| \neq 0`` for the constant coefficient vector ``e``, and the multiplier
``\alpha`` comes out at round-off whenever the right-hand side is mean-free — which it always
is here, since ``b^T \hat{\omega} = \int_\Omega \omega = 0``. `verify_spline.jl` measures it.
"""
struct PoissonMap{T, FT, MT} <: AbstractMatrix{T}
    F::FT
    M::MT
    b::Vector{T}
    N::Int
end

function PoissonMap(s::TensorSplineSpace{T}) where {T}
    N = nbasis(s)
    K = sparse(stiffness_matrix(s))
    M = sparse(mass_matrix(s))
    b = Vector{T}(basis_integrals(s))
    B = [K sparse(reshape(b, N, 1)); sparse(reshape(b, 1, N)) spzeros(T, 1, 1)]
    F = lu(B)
    return PoissonMap{T, typeof(F), typeof(M)}(F, M, b, N)
end

Base.size(P::PoissonMap) = (P.N, P.N)
Base.eltype(::PoissonMap{T}) where {T} = T

function Base.:*(P::PoissonMap{T}, û::AbstractVector) where {T}
    rhs = Vector{T}(undef, P.N + 1)
    @views mul!(rhs[1:(P.N)], P.M, û)
    rhs[P.N + 1] = zero(T)
    sol = P.F \ rhs
    return sol[1:(P.N)]
end

"""
    getindex(P::PoissonMap, i, j)

One entry of the assembled ``\\Lambda``, at the cost of a full solve.

Present only because `AbstractMatrix` promises it. Nothing on the hot path uses it, and
anything that does — assembling ``\\Lambda`` densely, or an implicit integrator's Jacobian —
should be reaching for a different formulation instead. See [`PoissonMap`](@ref).
"""
function Base.getindex(P::PoissonMap{T}, i::Int, j::Int) where {T}
    e = zeros(T, P.N)
    e[j] = one(T)
    return (P * e)[i]
end

## Two discrete Hamiltonians the package does not carry

@doc raw"""
    LinearHamiltonian(g)

The linear functional ``H(\hat{u}) = g \cdot \hat{u}``, with constant gradient `g` and
vanishing Hessian.

This is `eq:analytical_H`, ``H(u) = (h - h_\Omega, u)_{L^2}``, with
``g = \mathbb{M}(\hat{h} - h_\Omega e)``. PoissonBrackets has no linear
`DiscreteHamiltonian`: `MassCasimir` has exactly this structure, but its own docstring defines
it as the total mass ``\int_\Omega u \, dx``, so using it for a different linear functional
would put a false statement in the code. A candidate for the package, once a second caller
turns up.
"""
struct LinearHamiltonian{T} <: DiscreteHamiltonian{T}
    g::Vector{T}
end

function hamiltonian(H::LinearHamiltonian, ::DiscreteSpace, û::AbstractVector)
    dot(H.g, û)
end

gradient(H::LinearHamiltonian, ::DiscreteSpace, ::AbstractVector) = H.g

function hessian(H::LinearHamiltonian{T}, ::DiscreteSpace,
        û::AbstractVector) where {T}
    zeros(T, length(û), length(û))
end

@doc raw"""
    EllipticEnergy(Λ, M)

The reduced Euler energy `eq:Euler_H_periodic`,

```math
H(u) = \tfrac{1}{2} \int_\Omega |\nabla \phi|^2 dx = \tfrac{1}{2} (\phi, u)_{L^2}
     = \tfrac{1}{2} \hat{u}^T \mathbb{M} \Lambda \hat{u} ,
```

with ``\phi`` from `eq:Poisson-eq-periodic` and gradient ``\mathbb{M} \Lambda \hat{u} =
\mathbb{M} \hat{\phi}``.

`QuadraticHamiltonian(M * Λ)` would express the same thing, and is what an implicit
integrator would want, but its constructor checks symmetry with `isapprox(A, A')` — which
would assemble the dense ``\mathbb{M}\Lambda`` that [`PoissonMap`](@ref) exists to avoid. The
matrix *is* symmetric, ``\mathbb{M}\mathbb{K}^{-1}\mathbb{M}``, and `verify_spline.jl` checks
that on a small space where assembling it is affordable.

The Hessian is that same dense matrix and is not formed: the Section 4 runs use explicit
Runge-Kutta, as the manuscript does, so nothing asks for it.
"""
struct EllipticEnergy{T, LT, MT} <: DiscreteHamiltonian{T}
    Λ::LT
    M::MT
end

function EllipticEnergy(Λ::PoissonMap{T}, M) where {T}
    EllipticEnergy{T, typeof(Λ), typeof(M)}(Λ, M)
end

function hamiltonian(H::EllipticEnergy, ::DiscreteSpace, û::AbstractVector)
    dot(H.M * û, H.Λ * û) / 2
end

function gradient(H::EllipticEnergy, ::DiscreteSpace, û::AbstractVector)
    H.M * (H.Λ * û)
end

function hessian(::EllipticEnergy, ::DiscreteSpace, ::AbstractVector)
    throw(ArgumentError("the Hessian of an EllipticEnergy is dense and is not formed; " *
                        "the Section 4 runs are explicit, see EllipticEnergy"))
end

## The solver

@doc raw"""
    SplineTorus(n, p = 3)

A periodic tensor-product B-spline space of `n` cells and degree `p` per direction on
``\Omega = [0,2\pi]^2``, with the operators the Section 4 runs need assembled once.

`Λ` is the [`PoissonMap`](@ref); `M` the mass matrix and `Mfac` its factorisation; `b` the
basis integrals, so that ``\int_\Omega u_h = b \cdot \hat{u}``.

The number of degrees of freedom is `n^2` — a periodic B-spline basis has one function per
cell, independent of the degree.
"""
struct SplineTorus{T, ST <: TensorSplineSpace{T, 2}, LT, MT, FT}
    space::ST
    Λ::LT
    M::MT
    Mfac::FT
    b::Vector{T}
end

function SplineTorus(n::Int, p::Int = 3)
    s = TensorSplineSpace((n, n), p; L = DOMAIN_LENGTH)
    Λ = PoissonMap(s)
    M = sparse(mass_matrix(s))
    SplineTorus(s, Λ, M, mass_factorization(s), Vector{Float64}(basis_integrals(s)))
end

space(t::SplineTorus) = t.space
nbasis(t::SplineTorus) = nbasis(t.space)

function Base.show(io::IO, t::SplineTorus)
    print(io, "SplineTorus(n=", ncells(t.space)[1], ", p=", degree(t.space)[1],
        ", N=", nbasis(t.space), ")")
end

"``\\int_\\Omega u_h \\, dx``."
integrate(t::SplineTorus, û::AbstractVector) = dot(t.b, û)

"``(u_h, v_h)_{L^2}``, through the mass matrix."
l2inner(t::SplineTorus, û::AbstractVector, v̂::AbstractVector) = dot(û, t.M, v̂)

l2norm(t::SplineTorus, û::AbstractVector) = sqrt(max(l2inner(t, û, û), 0.0))
mean_value(t::SplineTorus, û::AbstractVector) = integrate(t, û) / DOMAIN_AREA

@doc raw"""
    spline_state(t, spec)

The initial state of run `spec`, as the pair `(ω̂, u_Ω)`.

The initial condition is ``L^2``-projected onto the space and its mean is then subtracted
exactly, in the same convention as [`spectral_state`](@ref): the state is the vorticity
``\omega = u - u_\Omega`` and ``u_\Omega`` is a constant of the motion carried alongside.

Subtracting the mean means subtracting ``u_\Omega e`` with `e` the constant coefficient
vector, which a periodic B-spline basis represents exactly by partition of unity.
"""
function spline_state(t::SplineTorus, spec::RunSpec)
    u₀ = initial_condition(spec)
    û = project(t.space, x -> u₀(x[1], x[2]))
    uΩ = mean_value(t, û)
    return (û .- uΩ, uΩ)
end

@doc raw"""
    spline_flow(t, spec)

The [`MetriplecticFlow`](@ref) of run `spec`: the metric bracket, the energy it must be
degenerate on, and the entropy ``S = \tfrac12 \int_\Omega \omega^2 dx``.

| run | bracket | generating field | energy |
|:--|:--|:--|:--|
| A1 | [`DoubleBracket`](@ref) | ``\hat{h}``, prescribed and fixed | [`LinearHamiltonian`](@ref) |
| A2 | [`DoubleBracket`](@ref) | ``\Lambda``, the elliptic solve | [`EllipticEnergy`](@ref) |
| A3, A4 | [`ProjectorBracket`](@ref) | ``\Lambda`` | [`EllipticEnergy`](@ref) |

The Poisson half is absent — `bracket = nothing` — which is the case every Section 4
experiment runs: the manuscript drops it in `eq:metric-system`.
"""
function spline_flow(t::SplineTorus, spec::RunSpec)
    S = QuadraticHamiltonian(Matrix(t.M))
    if spec.h !== nothing
        ĥ = project(t.space, x -> spec.h(x[1], x[2]))
        hΩ = mean_value(t, ĥ)
        H = LinearHamiltonian(Vector(t.M * (ĥ .- hΩ)))
        return MetriplecticFlow(t.space, DoubleBracket(t.space, ĥ), H, S)
    end
    H = EllipticEnergy(t.Λ, t.M)
    G = spec.bracket === :double ? DoubleBracket(t.space, t.Λ) :
        spec.bracket === :projector ? ProjectorBracket(t.space, t.Λ) :
        throw(ArgumentError("unknown bracket $(spec.bracket)"))
    return MetriplecticFlow(t.space, G, H, S)
end

@doc raw"""
    fixed_double_operator(t, spec)

The assembled weak-form operator ``\mathbb{A}_{KL} = \int_\Omega \partial_k \Phi_K
(X_h \otimes X_h)_{kl} \partial_l \Phi_L dx`` of run A1, where ``h`` is prescribed and
therefore ``X_h`` never changes.

# Why A1 does not go through `vectorfield`

The generic path re-forms ``X_h`` and runs the quadrature loop at every Runge-Kutta stage. It
has to, because for A2-A4 the generating field depends on the state. For A1 it does not, so
the whole operator is a constant sparse matrix and the right-hand side collapses to
``-\mathbb{M}^{-1} \mathbb{A} \hat{\omega}`` — two sparse operations instead of four
quadrature passes.

That is not a micro-optimisation. A1 runs at the manuscript's ``\Delta t = 10^{-4}``, which is
set by the explicit stability limit rather than by accuracy, so it takes 10⁵ steps against the
10⁴ of the other three. At roughly 10 ms per generic evaluation those 4·10⁵ stages would be
over an hour; assembled, they are minutes.

`verify_spline.jl` checks the two paths agree to round-off, which is what makes the
substitution legitimate rather than a second implementation.
"""
function fixed_double_operator(t::SplineTorus, spec::RunSpec)
    spec.h === nothing &&
        throw(ArgumentError("run $(spec.name) has no prescribed h"))
    ĥ = project(t.space, x -> spec.h(x[1], x[2]))
    X = hamiltonian_field(t.space, ĥ)
    𝔻 = Matrix{typeof(X[1])}(undef, 2, 2)
    for l in 1:2, k in 1:2

        𝔻[k, l] = X[k] .* X[l]
    end
    return sparse(tensor_weighted_matrix(t.space, 𝔻))
end

@doc raw"""
    spline_rhs(t, spec)

The right-hand side of run `spec` as a closure of ``\hat{\omega}`` alone.

A1 takes the assembled path of [`fixed_double_operator`](@ref); A2-A4 go through
[`MetriplecticFlow`](@ref)'s `vectorfield`, which is
``-\mathbb{G}(\hat{u}) \, \partial S / \partial \hat{u}``.
"""
function spline_rhs(t::SplineTorus, spec::RunSpec)
    if spec.h !== nothing
        A = fixed_double_operator(t, spec)
        return ω̂ -> -(t.Mfac \ (A * ω̂))
    end
    f = spline_flow(t, spec)
    return ω̂ -> vectorfield(f, ω̂)
end

@doc raw"""
    spline_step(rhs, ω̂, Δt)

One classical fourth-order Runge-Kutta step, matching the manuscript's own integrator.

Returns the new coefficient vector; the input is not modified.
"""
function spline_step(rhs, ω̂::AbstractVector, Δt)
    k1 = rhs(ω̂)
    k2 = rhs(ω̂ .+ (Δt / 2) .* k1)
    k3 = rhs(ω̂ .+ (Δt / 2) .* k2)
    k4 = rhs(ω̂ .+ Δt .* k3)
    return ω̂ .+ (Δt / 6) .* (k1 .+ 2 .* k2 .+ 2 .* k3 .+ k4)
end

@doc raw"""
    spline_grid(t, û, N)

The spline field sampled on the same ``N``-by-``N`` collocation grid a
[`SpectralTorus`](@ref) of size `N` uses, so that the two discretisations can be compared
pointwise.
"""
function spline_grid(t::SplineTorus, û::AbstractVector, N::Int)
    xs = collect(0:(N - 1)) .* (DOMAIN_LENGTH / N)
    return [evaluate(t.space, û, (xs[i], xs[j])) for i in 1:N, j in 1:N]
end
