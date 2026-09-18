@doc raw"""
    GradShafranovDisk(cells, degree = 3)

Section 5.5's **second** geometry, C2: the image of the unit disk under [`disk_map`](@ref),
discretised isogeometrically on a `PolarSplineSpace`.

The counterpart of [`GradShafranovBox`](@ref). It differs in three ways, and the third is the
one that decides what a relaxation on it converges to: the pole, the two measures, and the state
space.

# The pole, and why a tensor-product space cannot be used

`eq:mapping` sends the whole circle ``s = 0`` to one point, so a tensor-product basis on the
parameter square is **not even ``C^0``** there — a strictly stronger obstruction than the loss
of ``C^1`` one expects from a polar singularity, and one no refinement removes. A
`PolarSplineSpace` replaces the first two radial rows by three functions spanning the constants
and the two linear functions of the pseudo-Cartesian chart, which is ``C^0`` and ``C^1`` at the
pole by construction.

Nothing is regularised. There is no ``\varepsilon`` floor on ``s``, no puncture at the origin
and no modified basis near the pole; each of those would produce a run whose numbers are about
the regularisation.

# Two measures, and they are not interchangeable

The geometry enters through [`PulledBack`](@ref), twice, because the box's assembly uses two
different measures and the distinction is invisible on an unmapped domain:

| matrix | measure | why |
|:--|:--|:--|
| ``\mathbb{K}^\mu`` and ``\mathbb{B}`` | ``d\mu = dr\,dz/r`` | ``-\Delta^* = \operatorname{div}_\mu \nabla``, and the profile must carry the same weight or the eigenvalue is wrong by an order of magnitude |
| ``\mathbb{M}`` | ``dr\,dz`` | the source is ``j = u/r``, so ``\int u \varphi \, d\mu = \int j \varphi \, dr\,dz`` — the ``1/r`` is already in the state variable |

Both come from the same map and the same analytic Jacobian, and each `PulledBack` yields the
measure, the metric and the physical coordinates together, so the three places a weight is
needed cannot disagree. `scripts/verify_gradshafranov_disk.jl` measures what happens when they
do: swapping the two measures moves ``\lambda_h`` from ``0.00259704`` to ``0.000208`` one way
and ``0.0309`` the other.

# The boundary, and the state space it does NOT constrain

``\psi = 0`` on the rim ``s = 1``. The radial basis is clamped rather than Dirichlet-recombined
— `PolarSplineBasis` takes a plain clamped basis — so the condition is imposed by **dropping**
the last radial row from ``\Lambda``, the only functions that do not vanish at ``s = 1``. That
is the same elimination `GradShafranovBox`'s `:free` branch performs with a recombination
matrix.

**The state ``\hat{j}`` is not restricted, and this matters more than it looks.** Only
``\psi = \Lambda\hat{j}`` is; ``\hat{j}`` keeps the rim row and lives in the full polar
space. So the constant is in the state space — *exactly*, by the partition of unity the pole
triangle preserves — and so, to five and four digits, are ``r`` and ``z``. The bracket's mass
and momentum **Casimirs are therefore present**, and a relaxation converges to
``\delta S/\delta j = \lambda\psi + c\cdot x + \mu`` rather than to `eq:gs-ref`. This is
`GradShafranovBox`'s `:free` case, on a geometry that has no `:dirichlet` option, because a
homogeneous-Dirichlet *polar* spline space does not exist yet.

`scripts/run_c2.jl` measures it: residual `3.4162e-01` about ``\lambda\psi`` and `1.3341e-04`
about ``\lambda\psi + c\cdot x + \mu``, with ``\int u \, d\mu`` conserved to `1.1e-16`
where the box's Dirichlet space drifts 83 %. **An eigenvalue check cannot see any of it** — the
pencil ``(\mathbb{M}\Lambda, \mathbb{W})`` carries no mass constraint — which is why
[`gs_eigenvalue`](@ref) and the equilibrium check are unaffected and correct.
"""
struct GradShafranovDisk{T, ST <: PolarSplineSpace{T}, MT, PT}
    space::ST
    Λ::Matrix{T}
    M::MT
    MΛ::Matrix{T}
    W::Matrix{T}
    b::Vector{T}
    interior::Vector{Int}
    pb::PT
end

function GradShafranovDisk(cells::Tuple{Int, Int}, degree::Int = 3)
    s = PolarSplineSpace(cells, degree)

    F(x) = disk_map(x[1], x[2])
    DF(x) = disk_jacobian(x[1], x[2])

    # The two measures. `density` is a function of the *physical* point, so the Grad-Shafranov
    # weight is written in the coordinates it belongs to and composed with the map here.
    pbμ = PulledBack(s, F, DF; density = x -> gs_density(x))
    pbx = PulledBack(s, F, DF)

    Kμ = sparse(tensor_weighted_matrix(s, metric(pbμ)))
    M = sparse(weighted_matrix(s, measure(pbx), (0, 0), (0, 0)))

    keep = disk_interior(s)
    E = sparse(1:length(keep), keep, ones(length(keep)), length(keep), nbasis(s))'

    # Λ maps the current ĵ to the potential ψ̂, solving −Δ*ψ = j on the interior and extending
    # by zero. The factorisation is of the interior block, which is where the Dirichlet
    # condition lives; outside it Λ is zero, which is the condition itself.
    Λ = E * (cholesky(Symmetric(Matrix(E' * Kμ * E))) \ Matrix(E' * M))

    A = M * Λ

    # `gs_entropy_weight` already absorbs the measure: substituting u = rj into
    # S = ∫ u²/2(Cr²+D) dμ turns dμ and the Jacobian of the substitution into one power of r,
    # leaving σ = r/(Cr²+D) against the **plain** dr dz. So W takes `pbx` and not `pbμ`.
    # Using the μ measure here double-counts the 1/r; it does not raise, it moves the Rayleigh
    # quotient to 0.000208 and puts the initial state *below* λ_h, breaking the Poincaré floor.
    σ = [gs_entropy_weight(x[1]) for x in nodes(pbx)]
    W = Matrix(weighted_matrix(s, σ .* measure(pbx), (0, 0), (0, 0)))

    # `pbμ` is kept whole rather than unpacked into its measure and its nodes: `gs_flow` hands
    # it to `CollisionBracket`, which reads the frame `J⁻ᵀ` and the volume element from it as
    # well. Rebuilding it there would be a second object that could disagree with this one.
    GradShafranovDisk(s, Matrix(Λ), M, Matrix((A .+ A') ./ 2), (W .+ W') ./ 2,
        Vector{Float64}(M * ones(nbasis(s))), keep, pbμ)
end

GradShafranovDisk(cells::Int, degree::Int = 3) = GradShafranovDisk((cells, cells), degree)

space(disk::GradShafranovDisk) = disk.space
nbasis(disk::GradShafranovDisk) = nbasis(disk.space)

function Base.show(io::IO, disk::GradShafranovDisk)
    print(
        io, "GradShafranovDisk(cells=", ncells(disk.space), ", p=", degree(disk.space)[1],
        ", N=", nbasis(disk), ", interior=", length(disk.interior), ")")
end

@doc raw"""
    disk_interior(space::PolarSplineSpace)

The indices of the basis functions that vanish on the rim ``s = 1``, i.e. everything but the
last radial row.

A clamped radial basis has exactly one function that is nonzero at its right endpoint, so the
functions to drop are ``\Psi_K`` with ``K = 3 + (N_s-2) + (j-1)(N_s-2)``, one per angular
index. The three pole functions are never among them: they are supported on the first two
radial cells, which do not touch the rim.
"""
function disk_interior(s::PolarSplineSpace)
    radial, angular = bases(basis(s))
    Ns, Nθ = nbasis(radial), nbasis(angular)
    rim = Set(3 + (Ns - 2) + (j - 1) * (Ns - 2) for j in 1:Nθ)
    return [k for k in 1:nbasis(s) if !(k in rim)]
end

@doc raw"""
    gs_stiffness(space::PolarSplineSpace)
    gs_profile_matrix(space::PolarSplineSpace)
    gs_eigenvalue(space::PolarSplineSpace)

The mapped-disk counterparts of the box's three, assembled through [`PulledBack`](@ref) on
[`disk_map`](@ref):

```math
\mathbb{K}^\mu_{KL} = \int_\Omega \nabla \Psi_K \cdot \nabla \Psi_L \, d\mu ,
\qquad
\mathbb{B}_{KL} = \int_\Omega (Cr^2+D) \, \Psi_K \Psi_L \, d\mu ,
```

with the gradient the **physical** one — ``\nabla_x = J^{-T} \hat\nabla`` — and
``\lambda_h`` the smallest eigenvalue of ``\mathbb{K}^\mu \hat\psi = \lambda \mathbb{B}
\hat\psi`` on the interior.

``\lambda_h`` is **the number a relaxation run on this space converges to**, and it is not the
published ``0.002599`` nor exactly the continuum ``0.0025970``; see
[`GS_LAMBDA_DISK_CONTINUUM`](@ref). Measured here it is ``0.0025970351``, already to eight
figures at ``8 \times 16`` cubic cells — eigenvalues of a degree-``p`` isogeometric
discretisation converge at ``O(h^{2p})`` and the geometry is exact, since the Jacobian is
evaluated analytically rather than interpolated.
"""
function gs_stiffness(s::PolarSplineSpace)
    pb = PulledBack(s, x -> disk_map(x[1], x[2]), x -> disk_jacobian(x[1], x[2]);
        density = x -> gs_density(x))
    tensor_weighted_matrix(s, metric(pb))
end

function gs_profile_matrix(s::PolarSplineSpace)
    pb = PulledBack(s, x -> disk_map(x[1], x[2]), x -> disk_jacobian(x[1], x[2]);
        density = x -> gs_density(x))
    w = [herrnegger_mobility(x[1]) for x in nodes(pb)] .* measure(pb)
    weighted_matrix(s, w, (0, 0), (0, 0))
end

function gs_eigenvalue(s::PolarSplineSpace)
    keep = disk_interior(s)
    K = Symmetric(Matrix(gs_stiffness(s)[keep, keep]))
    B = Symmetric(Matrix(gs_profile_matrix(s)[keep, keep]))
    minimum(real, eigvals(K, B))
end

@doc raw"""
    gs_state(disk::GradShafranovDisk, spec::GSSpec)

The initial degrees of freedom of run `spec` on the mapped disk: the ``L^2`` projection of
``j = u_0/r``, with ``u_0`` the manuscript's Gaussian in the **physical** ``(r,z)`` coordinates
and the projection taken in the parameter square's own measure.

The Gaussian is a function of position, so it is sampled at [`nodes`](@ref) of the map and not
at the parameter nodes. Those are different points and nothing would complain.
"""
function gs_state(disk::GradShafranovDisk, spec::GSSpec)
    u₀ = initial_condition(spec)
    pb = PulledBack(disk.space, x -> disk_map(x[1], x[2]), x -> disk_jacobian(x[1], x[2]))
    project(disk.space, [u₀(x[1], x[2]) / x[1] for x in nodes(pb)])
end

@doc raw"""
    gs_flow(disk::GradShafranovDisk)

The [`MetriplecticFlow`](@ref) of the Section 5.5 relaxation on the mapped disk: the same three
objects as [`gs_flow`](@ref)`(::GradShafranovBox)` — the collision-like bracket, the
``\Delta^*`` energy and the Herrnegger-Maschke entropy — with the geometry supplied to the
bracket as a [`PulledBack`](@ref) rather than as a `density`.

# Why the pullback and not a density

A space's derivative tables are the **parameter** derivatives. On an unmapped domain those are
also the physical ones, which is why the box's keyword form is right and stays right. Here they
are not: the physical gradient is ``\nabla_x = J^{-T}\hat\nabla``, and a bracket that ignores
``J^{-T}`` is a bracket of the parametrisation. Handing `disk.pb` to the constructor fixes four
things from the one map — the measure, the frame the gradients and the perpendicular are read
in, the pairing the mass sandwich uses, and the points the mobility is sampled at.

**No structural check can tell the two apart.** Symmetry, positive semi-definiteness and the
degeneracy ``\mathbb{G}\,\partial H/\partial\hat{j} = 0`` hold for the framed bracket, for the
frameless one and for one whose frame is transposed. Passing the structural pass is necessary
and not sufficient on a mapped domain, so `scripts/verify_gradshafranov_disk.jl` measures the
covariance separately.

# The mobility is written in the physical radius

``M = Cr^2+D`` is sampled at ``F(\hat{x}_q)``, so `herrnegger_mobility(x[1])` needs no
composition with the map — a difference from the box, where the quadrature nodes already are
the physical points.

`mobility_derivative = 0` is not a shortcut: ``M`` depends on ``r`` alone, which is what makes
``\mathbb{G}`` quadratic in ``\hat{j}``, the dissipative residual an exact cubic polynomial and
its Jacobian analytic.

# The degeneracy

``\partial H/\partial\hat{j} = \mathbb{M}\Lambda\hat{j} = \mathbb{M}\hat\psi`` carries the
**physical** mass ``\mathbb{M} = \int\Psi_K\Psi_L\,dr\,dz``, and with a frame the sandwich pairs
with exactly that matrix rather than with the space's parameter-measure one. The two coincide on
an unmapped domain, which is why the box never needed the distinction.
"""
function gs_flow(disk::GradShafranovDisk)
    G = CollisionBracket(disk.space, disk.Λ, disk.pb;
        mobility = (x, u) -> herrnegger_mobility(x[1]),
        mobility_derivative = 0)
    MetriplecticFlow(disk.space, G, QuadraticHamiltonian(disk.MΛ),
        QuadraticHamiltonian(disk.W))
end

@doc raw"""
    gs_current(disk, ĵ)
    gs_ordinate(disk, ĵ)

The manuscript's ``u_h = r \, j_h`` and the scatter ordinate
``u_h/(Cr^2+D) = \sigma(r) j_h``, both on the quadrature grid.

``r`` is the **physical** radius, `nodes(disk.pb)`, not the radial parameter. On the box the two are
the same coordinate and the distinction does not arise; here it is the difference between
``[8,16]`` and ``[0,1]``.
"""
function gs_current(disk::GradShafranovDisk, ĵ::AbstractVector)
    [x[1] for x in nodes(disk.pb)] .* field(disk.space, ĵ, (0, 0))
end

function gs_ordinate(disk::GradShafranovDisk, ĵ::AbstractVector)
    [gs_entropy_weight(x[1]) for x in nodes(disk.pb)] .* field(disk.space, ĵ, (0, 0))
end

@doc raw"""
    gs_fit(disk, ĵ)

How close `ĵ` is to satisfying `eq:gs-ref`, ``u/(Cr^2+D) = \lambda\psi``, as
`(λ, residual, relative)` — the mapped-disk counterpart of
[`gs_fit`](@ref)`(::GradShafranovBox, ĵ)`, and read the same way.

The inner product is ``L^2(\mu)`` on the **physical** domain, so its weight is the parameter
quadrature's times the pulled-back measure. Both numbers are needed: ``\lambda`` is a
projection coefficient whose error is second order in `relative`, so it is already close for a
state that has barely moved.
"""
function gs_fit(disk::GradShafranovDisk, ĵ::AbstractVector)
    wμ = quadrature_weights(disk.space) .* measure(disk.pb)
    y = gs_ordinate(disk, ĵ)
    ψ = field(disk.space, disk.Λ * ĵ, (0, 0))
    λ = dot(wμ, y .* ψ) / dot(wμ, ψ .* ψ)
    res = sqrt(max(dot(wμ, (y .- λ .* ψ) .^ 2), 0.0))
    return (λ, res, res / sqrt(dot(wμ, y .^ 2)))
end

@doc raw"""
    gs_rayleigh(disk, ĵ)

The Rayleigh quotient ``S(\hat{j})/H(\hat{j})``, an upper bound on ``\lambda_h`` at every
admissible state and equal to it at an equilibrium. Identical in form to the box's, since it
reads the two Hamiltonian matrices and never touches a coordinate.
"""
function gs_rayleigh(disk::GradShafranovDisk, ĵ::AbstractVector)
    dot(ĵ, disk.W, ĵ) / dot(ĵ, disk.MΛ, ĵ)
end

@doc raw"""
    integrate(disk::GradShafranovDisk, ĵ)
    l2inner(disk::GradShafranovDisk, ĵ, v̂)
    l2norm(disk::GradShafranovDisk, ĵ)

``\int_\Omega j_h \, dx``, ``(j_h, v_h)_{L^2(dx)}`` and its norm — word for word the box's, and
for the same reasons, because ``\mathbb{M}`` and ``b`` already carry the map. Both are built
from `pbx`, the plain ``dr\,dz`` pullback, so nothing here needs to know the geometry a second
time.
"""
integrate(disk::GradShafranovDisk, ĵ::AbstractVector) = dot(disk.b, ĵ)

l2inner(disk::GradShafranovDisk, ĵ::AbstractVector, v̂::AbstractVector) = dot(ĵ, disk.M, v̂)

l2norm(disk::GradShafranovDisk, ĵ::AbstractVector) = sqrt(max(l2inner(disk, ĵ, ĵ), 0.0))

"""
    GradShafranovSolver

The two §5.5 discretisations, where a diagnostic reads only ``\\Lambda``, ``\\mathbb{W}``,
``\\mathbb{M}``, ``b`` and the space — all of which the box and the disk carry under the same
names. Used where widening a signature says the truth and a second method would only duplicate
it; where the two genuinely differ, as in [`gs_ordinate`](@ref), each keeps its own method.
"""
const GradShafranovSolver = Union{GradShafranovBox, GradShafranovDisk}
