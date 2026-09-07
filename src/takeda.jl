#
# Section 5.5 of *Metriplectic relaxation to equilibria*: the Grad-Shafranov problem itself,
# and the classical iterative solver the manuscript's reference eigenvalues come from.
#
# THIS FILE CARRIES NO RELAXATION.  Its whole reason for existing separately is that
# `lambda = 0.030302` and `lambda = 0.002599` are quoted in the manuscript as coming from "a
# standard Grad-Shafranov solver", i.e. from a computation that has nothing to do with the
# metriplectic flow.  Reproducing them from the flow would settle nothing: the flow is what is
# being checked.  So the problem, the measure, the Delta-star operator and Takeda's iteration
# live here, `gradshafranov.jl` builds the relaxation on top, and the two meet only where the
# run is compared against the number this file produced.
#
# WHY THERE ARE THREE ROUTES TO THE SAME EIGENVALUE, and why that is not redundancy:
#
#   * `takeda_iterate` is the manuscript's own kind of computation -- the classical iteration on
#     a uniform finite-difference grid, second order, the paper's `64 x 64 nodes`.  It is the
#     route that can be compared against the printed number on its own terms;
#   * `gs_eigenvalue` is a spline Galerkin solve of the same eigenvalue problem, on a completely
#     different discretisation.  Agreement between the two is evidence about the PROBLEM;
#   * `separable_eigenvalue` exploits the fact that on the rectangle the problem SEPARATES, and
#     reduces it to a one-dimensional Sturm-Liouville problem which is then solved to nine
#     digits.  That is what pins the continuum value, which neither of the other two can.
#

@doc raw"""
The profile constant ``C`` of the Herrnegger-Maschke entropy `eq:mm-entropy`, ``C = 0.6``.
"""
const HERRNEGGER_C = 0.6

@doc raw"""
The profile constant ``D`` of the Herrnegger-Maschke entropy `eq:mm-entropy`, ``D = 0.2``.
"""
const HERRNEGGER_D = 0.2

"The Section 5.5 rectangular domain in ``r``, ``[1, 7]``."
const GS_RADIAL = (1.0, 7.0)

"The Section 5.5 rectangular domain in ``z``, ``[-9.5, 9.5]``."
const GS_AXIAL = (-9.5, 9.5)

@doc raw"""
The eigenvalue ``\lambda = 0.030302`` the manuscript quotes for the rectangular Grad-Shafranov
case C1, "obtained from the standard iterative solver of the Grad-Shafranov equation, which has
also been implemented in FEniCS".

**It is a discrete number, not the continuum eigenvalue, and the difference is measurable.**
[`separable_eigenvalue`](@ref) pins the continuum value at
[`GS_LAMBDA_CONTINUUM`](@ref)` = 0.0302346260`, nine digits confirmed by three independent
routes; the printed number sits ``0.22\,\%`` **above** it. A conforming Galerkin eigenvalue
approaches its continuum limit from above, so the sign is the expected one and the size is the
discretisation error of the authors' own reference solver on their own ``64 \times 64`` grid —
tensor-product ``Q_1`` on that grid gives ``0.0302469`` and on 32 cells ``0.0302823``, which
brackets ``0.030302`` from either side. `verify_takeda.jl` measures that progression.

So a reproduction cannot recover the last three digits of this number, and claiming to would be
claiming to have recovered the authors' mesh, element and stopping tolerance from a
one-sentence description. What it can do, and what is checked, is reproduce it **to within its
own discretisation error** while pinning the continuum value it converges to.
"""
const GS_LAMBDA_RECTANGLE = 0.030302

@doc raw"""
The continuum eigenvalue of C1's problem, ``\lambda = 0.0302346260``.

Nine digits, and each one is measured rather than quoted: [`separable_eigenvalue`](@ref) is
degree- and mesh-independent to ``10^{-9}``, the two-dimensional Galerkin solve
[`gs_eigenvalue`](@ref) agrees without using separation at all, and the second-order finite
difference [`takeda_iterate`](@ref) Richardson-extrapolates onto it. See
[`GS_LAMBDA_RECTANGLE`](@ref) for why this is not the manuscript's printed number.
"""
const GS_LAMBDA_CONTINUUM = 0.0302346260

@doc raw"""
    gs_density(x)

The density ``m(r) = 1/r`` of the Section 5.5 measure ``d\mu = r^{-1} dr\,dz``
(`eq:mm-entropy`), as a function of a coordinate tuple `x = (r, z)`.

This is the one place the measure is written down. Everything that integrates over ``\Omega``
in Section 5.5 — the ``\Delta^*`` operator of [`gs_stiffness`](@ref), the inner moment
quadrature and the outer divergence of [`gs_flow`](@ref)'s bracket, the entropy — reaches it
through this function or through the weight it produces, so that a stray factor of ``r`` in one
of them cannot be a local slip. `verify_gradshafranov.jl` shows that each of the three would be
an ``O(1)`` error and not a tolerance-scale one.
"""
gs_density(x) = inv(x[1])

@doc raw"""
    herrnegger_mobility(r)

The mobility ``M(r) = C r^2 + D`` of the Herrnegger-Maschke entropy.

`eq:M-condition` fixes it: with ``s(r,y) = y^2 / 2(Cr^2+D)`` the second derivative is
``\partial_y^2 s = (Cr^2+D)^{-1}``, so ``M = 1/\partial_y^2 s = Cr^2 + D``. It depends on
``r`` alone and **not on the state**, which is what makes the dissipative residual of
[`gs_flow`](@ref) an exact cubic polynomial in the degrees of freedom.
"""
herrnegger_mobility(r) = HERRNEGGER_C * r^2 + HERRNEGGER_D

@doc raw"""
    herrnegger_profile(r, y)

The Grad-Shafranov profile ``f(r,y) = (C r^2 + D) \, y`` of the Herrnegger-Maschke equilibrium.

The manuscript defines ``f(r, \cdot) = \partial_y s(r, \cdot)^{-1}``, the inverse of the
entropy derivative in its second argument. For ``s = y^2/2(Cr^2+D)`` that derivative is
``y/(Cr^2+D)``, whose inverse is ``(Cr^2+D) y``, so ``f`` is **linear in ``y``** and
`eq:Grad-Shafranov-equation`

```math
-\Delta^* \psi = \lambda \, (C r^2 + D) \, \psi , \qquad \psi|_{\partial\Omega} = 0 ,
```

is a linear generalised **eigenvalue problem**. That is a property of this particular profile
and not of the method: [`takeda_iterate`](@ref) is written for a general `profile`, and it is
only because this one is linear that its fixed point is the smallest generalised eigenvalue and
can be cross-checked by an eigensolver.

Note that the manuscript also prints the *current* of the same equilibrium as
``(4\pi/c) J_\varphi = \lambda (C r + D/r) \psi``. That is the same statement — the state
variable is ``u = r (4\pi/c) J_\varphi``, so ``u = \lambda (Cr^2+D)\psi`` — and mistaking one
for the other is a plausible transcription error with an ``O(1)`` consequence:
`verify_takeda.jl` runs it as a control and it moves ``\lambda`` from ``0.0302`` to ``0.1404``.
"""
herrnegger_profile(r, y) = herrnegger_mobility(r) * y

## The Delta-star operator as a Galerkin matrix

@doc raw"""
    gs_stiffness(space)

The Galerkin matrix of ``-\Delta^*`` in the measure ``\mu``,

```math
\mathbb{K}^\mu_{IJ} = \int_\Omega \nabla \Phi_I \cdot \nabla \Phi_J \, d\mu
    = \int_\Omega \frac{1}{r} \big( \partial_r \Phi_I \partial_r \Phi_J
      + \partial_z \Phi_I \partial_z \Phi_J \big) \, dr \, dz ,
```

symmetric positive definite on a homogeneous-Dirichlet space.

**This is ``-\Delta^*`` and not a weighted Laplacian dressed up as one.** The Grad-Shafranov
operator ``\Delta^* = r \, \partial_r(r^{-1}\partial_r) + \partial_z^2`` is exactly the
divergence associated with ``d\mu = r^{-1} dr\,dz`` composed with the ordinary gradient,
``\Delta^* = \operatorname{div}_\mu \nabla`` with
``\operatorname{div}_\mu \varpi = m^{-1} \partial_i (m \varpi^i)`` and ``m = 1/r``; so for
``v`` vanishing on ``\partial\Omega``,

```math
\int_\Omega \nabla \psi \cdot \nabla v \, d\mu
    = - \int_\Omega v \operatorname{div}_\mu \nabla \psi \, d\mu
    = - \int_\Omega v \, \Delta^* \psi \, d\mu ,
```

which is the identity that makes the weighted stiffness matrix above the right operator with no
integration-by-parts boundary term and no non-symmetric first-order remainder. Dropping the
``1/r`` gives the ordinary Laplacian, which is a different operator with a different spectrum —
`verify_takeda.jl` measures the difference as a control.
"""
function gs_stiffness(s::TensorSplineSpace{T, 2}) where {T}
    ρ = gs_density.(quadrature_nodes(s))
    weighted_matrix(s, ρ, (1, 0), (1, 0)) .+ weighted_matrix(s, ρ, (0, 1), (0, 1))
end

@doc raw"""
    gs_profile_matrix(space)

The Galerkin matrix of the right-hand side of `eq:Grad-Shafranov-equation`,

```math
\mathbb{B}_{IJ} = \int_\Omega (C r^2 + D) \, \Phi_I \Phi_J \, d\mu ,
```

so that the eigenvalue problem ``-\Delta^*\psi = \lambda (Cr^2+D)\psi`` is
``\mathbb{K}^\mu \hat\psi = \lambda \mathbb{B} \hat\psi``.

Both matrices carry the measure. That is not a convention that cancels: ``\mathbb{K}^\mu`` with
``\mathbb{B}`` unweighted gives ``\lambda = 0.00607`` and the other way round ``0.120``, against
the ``0.0302`` of the consistent pair.
"""
function gs_profile_matrix(s::TensorSplineSpace{T, 2}) where {T}
    w = [herrnegger_mobility(x[1]) * gs_density(x) for x in quadrature_nodes(s)]
    weighted_matrix(s, w, (0, 0), (0, 0))
end

@doc raw"""
    gs_eigenvalue(space)

The smallest eigenvalue of ``\mathbb{K}^\mu \hat\psi = \lambda \mathbb{B} \hat\psi`` on
`space`, i.e. the Grad-Shafranov eigenvalue of the space — a two-dimensional Galerkin solve
that uses no separation of variables.

This is the number a relaxation run on the same space can actually reach, and it plays the part
`dirichlet_eigenvalue` plays for Section 5.4: a fitted ``\lambda`` that agreed with the
continuum value more closely than the space's own eigenvalue does would be a coincidence.

Dense, so it is a diagnostic rather than something to call in a loop; at the resolutions
Section 5.5 runs at, ``\mathbb{K}^\mu`` is a few hundred square.
"""
function gs_eigenvalue(s::TensorSplineSpace{T, 2}) where {T}
    K = Symmetric(Matrix(gs_stiffness(s)))
    B = Symmetric(Matrix(gs_profile_matrix(s)))
    minimum(real, eigvals(K, B))
end

@doc raw"""
    separable_eigenvalue(; cells = 64, degree = 3, mode = 1)

The continuum Grad-Shafranov eigenvalue of the **rectangle**, through separation of variables.

On ``[1,7] \times [-9.5, 9.5]`` the operator ``\Delta^*`` has no mixed term and the profile
``C r^2 + D`` depends on ``r`` alone, so ``\psi = R(r) Z(z)`` separates exactly. With
``Z_k = \sin\big(k\pi(z - z_a)/L_z\big)``, hence ``\partial_z^2 Z = -\kappa_k^2 Z`` and
``\kappa_k = k\pi/L_z``, the two-dimensional problem reduces to the Sturm-Liouville problem

```math
- \frac{d}{dr}\Big( \frac{1}{r} \frac{dR}{dr} \Big) + \frac{\kappa_k^2}{r} R
    = \lambda \, \frac{C r^2 + D}{r} \, R , \qquad R(1) = R(7) = 0 .
```

``\lambda`` is increasing in ``\kappa_k^2``, so ``k = 1`` gives the smallest eigenvalue and
`mode` exists only to demonstrate that — the second axial mode is ``0.03647``, far from the
first, which is what makes the relaxed state's axial structure unambiguous.

Solved by a one-dimensional spline Galerkin method, whose convergence in `degree` is what makes
this the sharp route: at degree 3 the answer is unchanged to ``10^{-9}`` between 32 and 128
cells, where the second-order finite differences of [`takeda_iterate`](@ref) need Richardson
extrapolation to reach five digits.

**Valid on a rectangle with an ``r``-only profile, and nowhere else.** C2's mapped domain has
neither property, which is why it has a `takeda_iterate` and no separable reference.
"""
function separable_eigenvalue(; cells::Int = 64, degree::Int = 3, mode::Int = 1)
    zmesh = UniformMesh(cells, GS_AXIAL)
    zs = TensorSplineSpace((zmesh,), degree, Dirichlet())
    κ² = sort(real.(eigvals(Symmetric(Matrix(stiffness_matrix(zs))),
        Symmetric(Matrix(mass_matrix(zs))))))[mode]

    rmesh = UniformMesh(cells, GS_RADIAL)
    rs = TensorSplineSpace((rmesh,), degree, Dirichlet())
    ρ = gs_density.(quadrature_nodes(rs))
    A = weighted_matrix(rs, ρ, (1,), (1,)) .+ κ² .* weighted_matrix(rs, ρ, (0,), (0,))
    B = weighted_matrix(
        rs, [herrnegger_mobility(x[1]) * gs_density(x)
             for x in quadrature_nodes(rs)],
        (0,), (0,))
    return (; λ = minimum(real, eigvals(Symmetric(Matrix(A)), Symmetric(Matrix(B)))), κ²)
end

## Takeda's iteration, on a uniform finite-difference grid

@doc raw"""
    TakedaGrid(nr, nz; radial = GS_RADIAL, axial = GS_AXIAL)

A uniform ``n_r \times n_z`` **node** grid on a rectangle, carrying the second-order
finite-volume discretisation of ``-\Delta^*`` in the measure ``\mu`` on its interior nodes.

The manuscript's rectangular case is "discretized by a uniform grid of ``64 \times 64``
nodes", so `TakedaGrid(64, 64)` is the grid the printed eigenvalue was computed on, up to the
element and the solver — which the manuscript does not state.

# Why finite volumes rather than the spline space the runs use

Because agreement between two discretisations of the same equation is evidence about the
equation, and agreement of a discretisation with itself is not. [`gs_eigenvalue`](@ref) solves
the same eigenvalue problem in the spline space; this solves it on a five-point stencil, and
the two share no assembly code.

The stencil is the conservative one for ``\operatorname{div}_\mu \nabla``, taken from the
weak form rather than from the differential operator: with ``\rho = 1/r`` evaluated at the
**cell faces** ``r_{i\pm1/2}``,

```math
(\mathbb{K}\psi)_{ij} = \frac{h_z}{h_r}
    \Big[ \frac{\psi_{ij} - \psi_{i+1,j}}{r_{i+1/2}}
        + \frac{\psi_{ij} - \psi_{i-1,j}}{r_{i-1/2}} \Big]
    + \frac{h_r}{h_z} \frac{2\psi_{ij} - \psi_{i,j+1} - \psi_{i,j-1}}{r_i} ,
```

which is symmetric positive definite exactly — the face values are shared between neighbours —
so the generalised eigenvalue problem it produces is symmetric and its Rayleigh quotient is
meaningful. Taking ``\rho`` at the nodes instead would break the symmetry and turn a real
eigenvalue problem into one whose imaginary parts have to be argued away.

`b` holds the lumped measure weights ``h_r h_z / r_i``, which is the mass matrix of the same
finite-volume scheme. Lumping is second order, matching the stencil, and it makes the profile
matrix diagonal.
"""
struct TakedaGrid{T}
    r::Vector{T}
    z::Vector{T}
    hr::T
    hz::T
    K::SparseMatrixCSC{T, Int}
    b::Vector{T}
    rint::Vector{T}
end

function TakedaGrid(nr::Int, nz::Int; radial = GS_RADIAL, axial = GS_AXIAL)
    (nr ≥ 3 && nz ≥ 3) || throw(ArgumentError(
        "a grid needs at least one interior node per direction, got $(nr) × $(nz)"))
    radial[1] > 0 || throw(ArgumentError(
        "the measure dμ = dr dz / r needs r > 0 on the closure, but the domain starts at " *
        "$(radial[1])"))
    hr = (radial[2] - radial[1]) / (nr - 1)
    hz = (axial[2] - axial[1]) / (nz - 1)
    r = [radial[1] + (i - 1) * hr for i in 1:nr]
    z = [axial[1] + (j - 1) * hz for j in 1:nz]

    ni, nj = nr - 2, nz - 2
    N = ni * nj
    index(i, j) = (j - 1) * ni + i          # interior (i,j) ↔ node (i+1, j+1)

    Is, Js, Vs = Int[], Int[], Float64[]
    push_entry!(a, c, v) = (push!(Is, a); push!(Js, c); push!(Vs, v))
    for j in 1:nj, i in 1:ni

        ii = i + 1
        rp, rm = (r[ii] + r[ii + 1]) / 2, (r[ii] + r[ii - 1]) / 2
        push_entry!(index(i, j), index(i, j),
            (hz / hr) * (1 / rp + 1 / rm) + (hr / hz) * 2 / r[ii])
        i < ni && push_entry!(index(i, j), index(i + 1, j), -(hz / hr) / rp)
        i > 1 && push_entry!(index(i, j), index(i - 1, j), -(hz / hr) / rm)
        j < nj && push_entry!(index(i, j), index(i, j + 1), -(hr / hz) / r[ii])
        j > 1 && push_entry!(index(i, j), index(i, j - 1), -(hr / hz) / r[ii])
    end

    rint = [r[i + 1] for j in 1:nj for i in 1:ni]
    TakedaGrid(r, z, hr, hz, sparse(Is, Js, Vs, N, N), (hr * hz) ./ rint, rint)
end

Base.length(g::TakedaGrid) = length(g.b)

function Base.show(io::IO, g::TakedaGrid)
    print(io, "TakedaGrid(", length(g.r), "×", length(g.z),
        " nodes, ", length(g.b), " interior)")
end

@doc raw"""
    takeda_iterate(g::TakedaGrid; profile = herrnegger_profile, axis = 1.0,
                   tol = 1e-13, maxiter = 2000)

The classical Grad-Shafranov iteration of Takeda and Tokuda, Eqs. (2.111)-(2.112), on the grid
`g`, returning `(; λ, ψ, iterations, λ_rayleigh, increment, converged)`.

Two steps per sweep, and the second is what determines ``\lambda``:

```math
-\Delta^* \tilde\psi^{k} = f(r, \psi^{k}) , \qquad
\lambda^{k+1} = \frac{\psi_{\mathrm{axis}}}{\max_\Omega \tilde\psi^{k}} , \qquad
\psi^{k+1} = \lambda^{k+1} \tilde\psi^{k} ,
```

so that ``-\Delta^*\psi^{k+1} = \lambda^{k+1} f(r, \psi^{k})`` holds at every sweep and
`eq:Grad-Shafranov-equation` at the fixed point. The eigenvalue is a *consequence* of a
normalisation — the value ``\psi_{\mathrm{axis}}`` on the magnetic axis — rather than something
solved for, which is exactly the structure of the classical scheme and the reason it is written
for a general `profile` rather than for the linear one C1 happens to have.

The elliptic solve is one sparse Cholesky factorisation of ``\mathbb{K}``, taken once before
the loop. ``\mathbb{K}`` does not depend on the state, so factorising it inside would be the
whole cost of the iteration.

`λ_rayleigh` is ``(\hat\psi^T \mathbb{K} \hat\psi) / (\hat\psi^T \mathbb{B} \hat\psi)`` with
``\mathbb{B}`` the lumped profile matrix, reported alongside because for a **linear** profile
the fixed point is the smallest generalised eigenvalue and the Rayleigh quotient is stationary
there: its error is second order in the eigenvector error where ``\lambda``'s is first order.
The two agreeing is therefore a convergence statement about the iteration and not a
restatement of it.

`increment` is the last ``|\lambda^{k+1} - \lambda^{k}|/\lambda^{k+1}``, and `converged` says
whether `tol` was reached. A run that hit `maxiter` returns rather than raising, because the
number it returns is still the honest answer to "where did the iteration get to".
"""
function takeda_iterate(g::TakedaGrid{T}; profile = herrnegger_profile, axis::T = one(T),
        tol::T = T(1e-13), maxiter::Int = 2000) where {T}
    F = cholesky(Symmetric(g.K))
    ψ = fill(axis, length(g))
    λ = zero(T)
    increment = T(Inf)
    converged = false
    iterations = maxiter

    for k in 1:maxiter
        ψ̃ = F \ (g.b .* profile.(g.rint, ψ))
        λnew = axis / maximum(ψ̃)
        increment = abs(λnew - λ) / abs(λnew)
        λ = λnew
        ψ = λ .* ψ̃
        if increment < tol
            iterations = k
            converged = true
            break
        end
    end

    B = g.b .* herrnegger_mobility.(g.rint)
    λ_rayleigh = dot(ψ, g.K, ψ) / dot(ψ, B .* ψ)
    return (; λ, ψ, iterations, λ_rayleigh, increment, converged)
end

@doc raw"""
    takeda_eigenvalue(g::TakedaGrid)

The smallest eigenvalue of ``\mathbb{K}\hat\psi = \lambda \mathbb{B}\hat\psi`` on `g` by a
**dense** symmetric eigensolve — the algorithmically independent check on
[`takeda_iterate`](@ref)'s fixed point.

The two must agree to round-off for the linear Herrnegger-Maschke profile and there is no
reason they should for any other, which is why this is a check on the iteration rather than a
replacement for it. Dense, hence usable on the coarse grids of a convergence study and not on
the ``64 \times 64`` one.
"""
function takeda_eigenvalue(g::TakedaGrid{T}) where {T}
    B = Diagonal(g.b .* herrnegger_mobility.(g.rint))
    minimum(real, eigvals(Symmetric(Matrix(g.K)), Matrix(B)))
end

@doc raw"""
    takeda_field(g::TakedaGrid, ψ)

The interior solution vector `ψ` written back onto the full ``n_r \times n_z`` node grid, with
zeros on the boundary — the shape a contour plot of ``\psi`` wants.
"""
function takeda_field(g::TakedaGrid{T}, ψ::AbstractVector) where {T}
    Ψ = zeros(T, length(g.r), length(g.z))
    ni = length(g.r) - 2
    for k in eachindex(ψ)
        i, j = mod1(k, ni), (k - 1) ÷ ni + 1
        Ψ[i + 1, j + 1] = ψ[k]
    end
    return Ψ
end

@doc raw"""
    takeda_order(λs, ns)

The observed convergence order of the eigenvalues `λs` computed on the grids `ns`, as the
least-squares slope of ``\log|\lambda_n - \lambda_{\mathrm{ref}}|`` against ``\log n`` with
``\lambda_{\mathrm{ref}}`` the Richardson extrapolation of the two finest.

Second order is what the stencil of [`TakedaGrid`](@ref) has, and measuring it is what says the
extrapolation below is legitimate rather than a way of getting a digit that is not there.
"""
function takeda_order(λs::AbstractVector, ns::AbstractVector)
    length(λs) == length(ns) ≥ 3 || throw(ArgumentError(
        "an order needs at least three grids, got $(length(ns))"))
    ref = takeda_extrapolate(λs[end - 1], λs[end], ns[end - 1], ns[end])
    x = log.(Float64.(ns[1:(end - 1)]))
    y = log.(abs.(λs[1:(end - 1)] .- ref))
    x̄, ȳ = sum(x) / length(x), sum(y) / length(y)
    return -sum((x .- x̄) .* (y .- ȳ)) / sum((x .- x̄) .^ 2)
end

@doc raw"""
    takeda_extrapolate(λa, λb, na, nb)

Richardson extrapolation of a second-order quantity from the grids `na` and `nb`,
``\lambda \approx \lambda_b + (\lambda_b - \lambda_a)/((n_b/n_a)^2 - 1)``.

Second order is [`TakedaGrid`](@ref)'s stencil, and [`takeda_order`](@ref) measures it rather
than assuming it — an extrapolation applied at the wrong order manufactures digits.
"""
function takeda_extrapolate(λa, λb, na, nb)
    q = (nb / na)^2
    return λb + (λb - λa) / (q - 1)
end

## C2's geometry, and its reference eigenvalue

@doc raw"""
The constants of `eq:mapping`, the map that carries the unit disk onto C2's domain:
``e = 1.4``, ``\varepsilon = 0.3``, ``a = 4``, ``b = 3``, ``c = 6.3``, and
``\xi = 1/\sqrt{1 - \varepsilon^2/4}``.

A slightly modified form of the map of Zoni and Güçlü, which the manuscript cites — the same
paper whose subject is the ``C^1`` polar-spline construction C2's relaxation would need.
"""
const DISK_MAP = (e = 1.4, ε = 0.3, a = 4.0, b = 3.0, c = 6.3, ξ = 1 / sqrt(1 - 0.3^2 / 4))

@doc raw"""
    disk_map(s, θ)

`eq:mapping`: the point of C2's domain that the unit-disk point ``s e^{i\theta}`` maps to,
as ``(r, z)``.

```math
r = a \Big[ b + \frac{1}{\varepsilon}
    \Big( 1 - \sqrt{1 + \varepsilon(\varepsilon + 2 s \cos\theta)} \Big) \Big] , \qquad
z = c \, \frac{e \, \xi \, s \sin\theta}
             {2 - \sqrt{1 + \varepsilon(\varepsilon + 2 s \cos\theta)}} .
```

The image is ``r \in [8, 16]``, ``z \in [-9.749139, 9.749139]``, of area ``114.777`` — the
Jacobian integrated over the parameter disk, and comparable to C1's rectangle, which has area
114 — with the pole ``s = 0`` at ``r = 11.412925``. That last number is what makes C2's initial
condition consistent: its Gaussian is centred at ``r_0 = 12``, just outboard of the axis.
None of the four appears in the manuscript; all four are what its constants produce.

**The map degenerates at ``s = 0``**: the whole circle ``s = 0`` collapses to one point, so
the parametrisation is not a diffeomorphism there, and that is the whole reason C2's relaxation
is deferred rather than run. See [`disk_eigenvalue`](@ref) for why the *eigenvalue* is
nevertheless computable without regularising anything.
"""
function disk_map(s, θ)
    (; e, ε, a, b, c, ξ) = DISK_MAP
    q = sqrt(1 + ε * (ε + 2s * cos(θ)))
    return (a * (b + (1 - q) / ε), c * e * ξ * s * sin(θ) / (2 - q))
end

@doc raw"""
The eigenvalue ``\lambda = 0.002599`` the manuscript quotes for C2, "computed by a standard
Grad-Shafranov solver".

As for [`GS_LAMBDA_RECTANGLE`](@ref) this is a discrete number: the continuum value is
[`GS_LAMBDA_DISK_CONTINUUM`](@ref)` = 0.0025970`, and the printed one sits ``0.075\,\%`` above
it, from above as a conforming Galerkin eigenvalue must. A ``P_1`` triangulation of 64 radial
by 128 angular cells gives ``0.00259909``, which rounds to the printed value.
"""
const GS_LAMBDA_DISK = 0.002599

@doc raw"""
The continuum Grad-Shafranov eigenvalue of C2's mapped domain, ``\lambda = 0.0025970``.

Measured, not quoted: [`disk_eigenvalue`](@ref) converges onto it at second order over six
refinements and Richardson-extrapolates to this value from the two finest. Five digits rather
than the rectangle's nine, because there is no separation of variables here and no
one-dimensional reference to be had — see [`separable_eigenvalue`](@ref).
"""
const GS_LAMBDA_DISK_CONTINUUM = 0.0025970

@doc raw"""
    DiskTriangulation(n, m)

A triangulation of C2's domain: `n` rings by `m` angular sectors of the unit disk, pushed
through [`disk_map`](@ref), with the pole as node 1.

`triangles` are index triples, `boundary` the node indices on ``\partial\Omega`` (the image of
``s = 1``). The innermost ring is `m` triangles with a vertex at the pole; every outer ring is
`2m` triangles.
"""
struct DiskTriangulation{T}
    r::Vector{T}
    z::Vector{T}
    triangles::Vector{NTuple{3, Int}}
    boundary::Vector{Int}
end

function DiskTriangulation(n::Int, m::Int)
    (n ≥ 1 && m ≥ 3) || throw(ArgumentError(
        "a disk needs at least one ring and three sectors, got $(n) × $(m)"))
    (r₀, z₀) = disk_map(0.0, 0.0)
    r, z = [r₀], [z₀]
    for i in 1:n, j in 1:m

        (rj, zj) = disk_map(i / n, 2π * (j - 1) / m)
        push!(r, rj)
        push!(z, zj)
    end
    node(i, j) = 1 + (i - 1) * m + mod(j - 1, m) + 1
    tris = NTuple{3, Int}[(1, node(1, j), node(1, j + 1)) for j in 1:m]
    for i in 1:(n - 1), j in 1:m

        push!(tris, (node(i, j), node(i + 1, j), node(i + 1, j + 1)))
        push!(tris, (node(i, j), node(i + 1, j + 1), node(i, j + 1)))
    end
    DiskTriangulation(r, z, tris, [node(n, j) for j in 1:m])
end

Base.length(t::DiskTriangulation) = length(t.r)

function Base.show(io::IO, t::DiskTriangulation)
    print(io, "DiskTriangulation(", length(t.r), " nodes, ", length(t.triangles),
        " triangles, ", length(t.boundary), " on ∂Ω)")
end

@doc raw"""
    disk_area(t::DiskTriangulation)

The area of the triangulated domain, ``\sum_e |T_e|`` — the check that the map has been
transcribed correctly.

It approaches the true area ``114.777`` **from below**, at second order, because a polygon
inscribed in a convex boundary is smaller than the region it approximates: measured, `114.605`
at 32 × 64, `114.734` at 64 × 128 and `114.766` at 128 × 256. So a check against the exact
value needs a tolerance that says which mesh it was measured on.
"""
function disk_area(t::DiskTriangulation{T}) where {T}
    s = zero(T)
    for (i, j, k) in t.triangles
        s += abs((t.r[j] - t.r[i]) * (t.z[k] - t.z[i]) -
                 (t.r[k] - t.r[i]) * (t.z[j] - t.z[i])) / 2
    end
    return s
end

@doc raw"""
    disk_matrices(t::DiskTriangulation)

The ``P_1`` Galerkin matrices of C2's eigenvalue problem on the interior nodes of `t`,

```math
\mathbb{K}_{IJ} = \int_\Omega \nabla \phi_I \cdot \nabla \phi_J \, d\mu , \qquad
\mathbb{B}_{IJ} = \int_\Omega (C r^2 + D) \, \phi_I \phi_J \, d\mu ,
```

returned as `(K, B, free)` with `free` the interior node indices.

Both weights are taken at the element **centroid**, which is exact for a linear weight and
second order in general — the same order as ``P_1``'s eigenvalue error, so it is not the term
that limits the result. `verify_takeda.jl` measures the order rather than taking that on trust.
"""
function disk_matrices(t::DiskTriangulation{T}) where {T}
    N = length(t.r)
    Is, Js, Ks, Bs = Int[], Int[], T[], T[]
    for e in t.triangles
        x = (t.r[e[1]], t.r[e[2]], t.r[e[3]])
        y = (t.z[e[1]], t.z[e[2]], t.z[e[3]])
        det = (x[2] - x[1]) * (y[3] - y[1]) - (x[3] - x[1]) * (y[2] - y[1])
        area = abs(det) / 2
        iszero(area) && continue
        # ∇φ of the three linear shape functions, constant on the element
        gx = (y[2] - y[3], y[3] - y[1], y[1] - y[2]) ./ det
        gy = (x[3] - x[2], x[1] - x[3], x[2] - x[1]) ./ det
        rc = (x[1] + x[2] + x[3]) / 3
        # ∫ dμ over the element, and the same against the profile: dμ = dr dz / r
        wk = area / rc
        wb = area * herrnegger_mobility(rc) / rc
        for p in 1:3, q in 1:3

            push!(Is, e[p])
            push!(Js, e[q])
            push!(Ks, wk * (gx[p] * gx[q] + gy[p] * gy[q]))
            # the exact P₁ element mass matrix, area/12 × (1 + δ_pq), scaled by the weight
            push!(Bs, wb * (p == q ? 1 / 6 : 1 / 12))
        end
    end
    free = setdiff(1:N, t.boundary)
    return (sparse(Is, Js, Ks, N, N)[free, free],
        sparse(Is, Js, Bs, N, N)[free, free], free)
end

@doc raw"""
    disk_eigenvalue(n, m; tol = 1e-13, maxiter = 2000)

The smallest eigenvalue of ``-\Delta^*\psi = \lambda (Cr^2+D)\psi`` on C2's mapped domain, by
inverse iteration on the ``P_1`` matrices of [`disk_matrices`](@ref).

# Why a triangulation of the physical domain, and why this is not what C2's relaxation needs

The obvious route is the isogeometric one: solve on the parameter square
``(s,\theta) \in [0,1] \times [0,2\pi)``, periodic in ``\theta``, pulling the metric back
through [`disk_map`](@ref). It cannot be taken here, and the obstruction is not a matter of
effort: at ``s = 0`` the map collapses the whole circle to a point, so a function on the
parameter square is single-valued at the pole only if its ``\theta``-dependence there is
constrained, and a tensor-product spline space provides no such constraint. Enforcing it needs
a **polar-spline** construction — the first two rows of the ``s``-basis replaced by a
three-function pole triangle — which is what Zoni and Güçlü build and what SimpleSplines does
not have.

A triangulation of the **physical** domain has no such problem: the pole is an ordinary node,
``P_1`` needs only ``C^0``, and nothing is regularised, punctured or floored. So the reference
eigenvalue is computable exactly as the manuscript's own solver computes it.

What this does *not* give is C2's relaxation. That needs ``\nabla\psi`` on the space the state
lives in, sampled at the quadrature points of a [`CollisionBracket`](@ref) — and
PoissonBrackets' only two-dimensional space is the tensor-product
[`TensorSplineSpace`](@ref). There is no ``P_1`` triangular `DiscreteSpace`, so this mesh
cannot carry the flow, and the parameter square cannot carry the pole. Either obstruction alone
defers C2; `CHANGELOG.md` records both.
"""
function disk_eigenvalue(n::Int, m::Int; tol = 1e-13, maxiter::Int = 2000)
    t = DiskTriangulation(n, m)
    (K, B, free) = disk_matrices(t)
    F = cholesky(Symmetric(K))
    x = ones(length(free))
    λ = 0.0
    for _ in 1:maxiter
        y = F \ (B * x)
        λnew = dot(x, B, x) / dot(x, B * y)
        conv = abs(λnew - λ) < tol * abs(λnew)
        λ = λnew
        x = y ./ maximum(abs, y)
        conv && break
    end
    return (; λ, dof = length(free), area = disk_area(t))
end
