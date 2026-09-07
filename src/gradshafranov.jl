#
# Section 5.5 of *Metriplectic relaxation to equilibria*: the metriplectic relaxation of the
# Grad-Shafranov problem `takeda.jl` states, on the rectangle C1.
#
# WHAT CHANGES FROM SECTION 5.4, and the first item is the whole file:
#
#   * the measure is `dmu = dr dz / r`, not `dx`, and it appears in THREE places -- the inner
#     moment quadrature of the bracket, the outer div_mu assembly, and Delta-star.  A stray
#     factor of r in any one of them is an O(1) error; `verify_gradshafranov.jl` measures each;
#   * THE STATE VARIABLE IS `j = u/r`, NOT `u`.  That is a choice, it is forced, and
#     `GradShafranovBox` explains it.  With it the whole mu-weighting collapses onto exactly the
#     Section 5.4 machinery: one plain mass matrix, one mu-weighted stiffness matrix, and a
#     bracket whose `density` carries the measure;
#   * the entropy is `s = y^2/2(Cr^2+D)`, still QUADRATIC in `y`, so unlike B3 there is no
#     positivity requirement and no floor.  In the `j` variable it is a quadratic form with the
#     x-dependent weight `sigma = r/(Cr^2+D)`;
#   * the mobility `M = Cr^2+D` does NOT depend on the state, so `mobility_derivative = 0`, the
#     dissipative residual is an exact cubic polynomial in the degrees of freedom, and the
#     Jacobian is analytic.  Finite differences have no business here.
#
# Everything else -- the collision bracket, implicit midpoint, the Newton solve, the dense
# Jacobian -- is what `euler.jl` already runs, unchanged.
#

@doc raw"""
    gs_entropy_weight(r)

The weight ``\sigma(r) = r / (C r^2 + D)`` of the Section 5.5 entropy in the state variable
``j = u/r``.

`eq:mm-entropy` gives ``\fun{S}(u) = \int_\Omega u^2 / 2(Cr^2+D) \, d\mu`` with
``d\mu = r^{-1} dr\,dz``. Substituting ``u = rj`` turns the measure and the Jacobian of the
substitution into a single power of ``r``:

```math
\fun{S} = \int_\Omega \frac{r^2 j^2}{2(Cr^2+D)} \frac{dr\,dz}{r}
        = \frac{1}{2} \int_\Omega \frac{r}{Cr^2+D} \, j^2 \, dr \, dz ,
```

so ``S`` is a quadratic form on the coefficients with an ``x``-dependent weight, and its
derivative in the **plain** ``L^2(dx)`` pairing is

```math
\frac{\delta \fun{S}}{\delta j} = \sigma(r) \, j = \frac{u}{C r^2 + D} ,
```

which is the ordinate of the manuscript's scatter plot and the left-hand side of `eq:gs-ref`.
See [`GradShafranovBox`](@ref) for why the plain pairing is the one that works.
"""
gs_entropy_weight(r) = r / (HERRNEGGER_C * r^2 + HERRNEGGER_D)

@doc raw"""
    GSSpec

One of the Section 5.5 runs: its geometry, its initial condition, and the time step and final
time this reproduction chose for it.

Everything about the initial condition is the manuscript's, `eq:sv-ic` with the **squared**
widths §5 prints — see [`gaussian_w2`](@ref). Neither ``\Delta t`` nor ``T`` is: §5 states no
time step, no final time, no stopping criterion, no nonlinear-solver description and no
tolerance. Both are recorded choices and [`SECTION55_RUNS`](@ref) says what they produced.

`cells` is the per-axis cell count, which has no counterpart in the manuscript either — the
spline discretisation is this reproduction's deviation, and the resolution is capped by the
bracket's ``O(N^3)`` Jacobian rather than chosen for accuracy. See [`GradShafranovBox`](@ref).
"""
struct GSSpec
    name::String
    section::String
    gaussian::Gaussian{Float64}
    cells::Tuple{Int, Int}
    degree::Int
    Δt::Float64
    T::Float64
end

"""
    initial_condition(spec::GSSpec)

The initial **current** ``u_0`` of run `spec` as a callable `(r, z)` — the Gaussian exactly as
printed. The state variable is ``j = u/r``; [`gs_state`](@ref) is where the division happens.
"""
initial_condition(spec::GSSpec) = spec.gaussian

@doc raw"""
The runs of Section 5.5, keyed by name.

| run | figs | domain | ``x_0`` | ``w_1^2`` | ``w_2^2`` | ``N`` |
|:--|:--|:--|:--|:--|:--|:--|
| `c1` | `gsr_*` | ``[1,7] \times [-9.5, 9.5]`` | ``(4, 0)`` | ``0.5`` | ``3.2`` | ``1`` |

The entropy and the profile are the same for every §5.5 run — ``s = y^2/2(Cr^2+D)`` with
``C = 0.6``, ``D = 0.2`` — so they are constants of `takeda.jl` rather than fields here.

**C2, the mapped disk, is not in this table.** It is deferred, and the reason is in
`CHANGELOG.md`: its domain is the image of the unit disk, the isogeometric route brings a polar
singularity at ``s = 0`` where the standard tensor-product basis is not ``C^1``, and a
tensor-product space on a box has no way to represent that without a polar-spline
construction. Nothing here is quietly regularised to make it fit.

# The recorded choices, and what they produced

| run | ``\Delta t`` | ``T`` | steps | cells | degree | ``N`` |
|:--|--:|--:|--:|--:|--:|--:|
| `c1` | 0.0625 | 25 | 400 | ``18 \times 21`` | 2 | 378 |

**``\Delta t`` is two orders of magnitude smaller than §5.4's, and §5.4's reasoning does not
carry over.** B1 runs at ``\Delta t = 1`` because its relaxation has no fast scale; C1's does.
The bracket's own magnitude is what sets it: ``m_0 = \int_\Omega M \, d\mu \approx 281`` on this
domain against ``1`` on the unit square, since ``\Omega`` is 114 units of area rather than one
and ``M = Cr^2+D \approx 10`` rather than ``1``. Measured at ``18 \times 21`` cells, the entropy
falls to within ``6 \times 10^{-5}`` of its floor by ``t = 5`` and to round-off by ``t = 18``.

**Crank-Nicolson makes too large a step look like a slower relaxation, which is the trap here.**
Its amplification factor for a mode of rate ``\lambda`` is ``(1-\lambda\Delta t/2) /
(1+\lambda\Delta t/2)``, which for ``\lambda \Delta t \gg 1`` is ``\approx -1 + 4/(\lambda\Delta
t)``: a stiff mode is damped by only ``4/(\lambda\Delta t)`` per step, so after ``n`` steps the
decay exponent is ``4t/(\lambda \Delta t^2)`` and the *apparent* relaxation time scales as
``\Delta t^2``. Measured, and this is the number that settles it: the excess entropy at
``t = 200`` is ``1.2 \times 10^{-5}`` at ``\Delta t = 2``, ``2 \times 10^{-13}`` at
``\Delta t = 1``, and at ``\Delta t = 0.25`` it is at round-off by ``t = 20``. A run at
``\Delta t = 2`` is not a slowly relaxing system; it is an unresolved one. Every claim it makes
about ``H``, monotonicity and the equilibrium still holds — the fixed point does not depend on
``\Delta t`` — but the trace is not the trajectory. That is what makes the step-size study below
load-bearing rather than a formality: nothing in the run's own output would have flagged it.

**The step at which the trajectory is resolved depends on the MESH, which is easy to miss.** The
stiffest mode's rate grows like ``h^{-2}``, so ``\lambda \Delta t \gg 1`` sets in at a smaller
``\Delta t`` on a finer space, and a step-size study run on a coarse mesh does not transfer.
Measured at ``t = 2.5`` on this space, the relative ``L^2`` difference between successive
halvings runs

| ``\Delta t`` | 0.5 → 0.25 | 0.25 → 0.125 | 0.125 → 0.0625 | 0.0625 → 0.03125 |
|:--|--:|--:|--:|--:|
| difference | ``1.33 \times 10^{-1}`` | ``1.35 \times 10^{-2}`` | ``9.8 \times 10^{-4}`` | ``2.7 \times 10^{-4}`` |
| ratio | | 9.8 | 13.8 | 3.6 |

and the **second-order regime begins at ``\Delta t = 0.0625``** — the ratio is 3.6 there and 10
to 14 above it, which is the stiff-mode error dying off rather than a convergence order. So
``\Delta t = 0.0625`` is the recorded step: it is the largest one whose difference from
``\Delta t/2`` is both small (``2.7 \times 10^{-4}``) and *asymptotic*. At ``\Delta t = 0.25``
the same difference is ``1.4 \times 10^{-2}``, which is a fourteenth of the state and not a
resolved trajectory.

``T = 25`` follows from the same stopping criterion §5.4 used, the excess entropy: ``S`` reaches
``\lambda_h H_0`` to twelve digits and the scatter residual reaches the projection-error floor
exactly. `run_c1.jl` re-measures both, and the step-size ratio, rather than trusting this
paragraph — and it does so **at ``t = 2.5``** rather than after five steps, because a comparison
taken in the first few steps of a fast transient measures the first step's error and not the
trajectory's.

**The entropy is monotone at every step size tested, and for this entropy that is a theorem
rather than a measurement.** ``S`` is *quadratic* in ``\hat{j}``, so
``S(\hat{j}^{n+1}) - S(\hat{j}^n) = \partial S/\partial \hat{j}(\bar{j}) \cdot
(\hat{j}^{n+1}-\hat{j}^n)`` holds exactly at the midpoint, and the increment is
``-\Delta t \, \mathbb{G}(\bar{j})\partial S/\partial\hat{j}(\bar{j})``, so the difference is
``-\Delta t \, g^T \mathbb{G} g \le 0`` by positive semi-definiteness alone. §5.4's B3 needed a
small step for monotonicity because ``y\log y`` is not quadratic; C1 does not, and its
``\Delta t`` is an accuracy choice throughout.

**The mesh is anisotropic, and matched to the initial condition rather than to the domain.**
``\Omega`` is three times longer in ``z`` than in ``r``, but the Gaussian is *also* wider in
``z``: ``w_1 = \sqrt{0.5} = 0.707`` against ``w_2 = \sqrt{3.2} = 1.789``. Cells per e-folding
width is what matters, so ``n_1/n_2 \approx (L_1/w_1)/(L_2/w_2) = 8.49/10.62 = 0.80``, and
``18 \times 21`` is that ratio at a degree-of-freedom count the ``O(N^3)`` Jacobian allows.
Measured on this machine at ``\Delta t = 2``, one step costs 0.67 s at ``N = 120``, 3.1 s at
``N = 224`` and 25 s at ``N = 378``; the manuscript's ``64 \times 64`` would be ``N = 4096`` and
days. Per-step cost falls sharply with the step size, because the Newton iteration is what
varies — at ``N = 378`` and ``\Delta t = 0.0625`` it is **2.1 s**, so the 400 steps of the
recorded run are some 14 minutes rather than the three hours the ``\Delta t = 2`` figure would
suggest. The cells are strongly anisotropic in *shape* — ``h_z/h_r = 2.7`` — which is the right
answer and not a defect.
"""
const SECTION55_RUNS = Dict(
    "c1" => GSSpec("c1", "5.5", gaussian_w2((4.0, 0.0), (0.5, 3.2), 1.0),
    (18, 21), 2, 0.0625, 25.0))

"The Section 5.5 runs, in the order the manuscript presents them."
const SECTION55_ORDER = ("c1",)

## The discretisation

@doc raw"""
    GradShafranovBox(cells, degree = 2; state = :dirichlet)

A tensor-product B-spline space on the Section 5.5 rectangle
``\Omega = [1,7] \times [-9.5,9.5]``, with the operators the relaxation needs assembled once:
``\Lambda`` the ``\Delta^*`` solution operator, ``\mathbb{M}`` the plain mass matrix,
``\mathbb{M}\Lambda`` the matrix of the energy, ``\mathbb{W}`` the matrix of the entropy, and
``b`` the basis integrals.

`cells` is a pair ``(n_r, n_z)`` or a single number for both.

# The state variable is ``j = u/r``, and that is what makes the discrete structure exact

The manuscript's state is ``u = (4\pi/c) r J_\varphi``, and its functional derivatives are
taken "with respect to the ``L^2``-product with the measure ``\mu``". This reproduction evolves
``j = u/r = (4\pi/c) J_\varphi`` — the toroidal current density itself — and takes its
derivatives in the **plain** ``L^2(dx)`` product. The two formulations are the same problem,
because

```math
\int_\Omega f \, \delta u \, d\mu = \int_\Omega f \, r \, \delta j \, \frac{dr\,dz}{r}
    = \int_\Omega f \, \delta j \, dx ,
```

so every functional derivative is **numerically the same function** in both:
``\delta \fun{H}/\delta j = \psi`` and ``\delta \fun{S}/\delta j = u/(Cr^2+D)``, exactly as
`eq:GradShafranov-S-H` and `eq:gs-ref` state them. The measure has not been dropped; it has
been moved out of the pairing and into the operators, where `eq:div-grad-ev` needs it anyway.

**This is not presentational.** [`CollisionBracket`](@ref) forms
``\mathbb{G} = \mathbb{M}^{-1}\mathbb{A}\mathbb{M}^{-1}`` with ``\mathbb{M}`` the space's own
mass matrix, and energy conservation is the statement
``\mathbb{G}\,\partial H/\partial \hat{u} = 0``, which holds exactly when
``\mathbb{M}^{-1}\partial H/\partial \hat{u}`` is the bracket's generating field ``\hat\psi``.
In the ``u`` formulation the elliptic solve is
``\Lambda = (\mathbb{K}^\mu)^{-1}\mathbb{M}^\mu`` with a ``\mu``-**weighted** mass matrix on the
right, so ``\partial H/\partial\hat{u} = \mathbb{M}^\mu\hat\psi`` and the degeneracy would need
``\mathbb{M}^\mu`` in the sandwich, which no tensor-product mass operator provides. In the ``j``
formulation the weak form of `eq:GradShafranov-psi` is

```math
\int_\Omega \nabla \psi_h \cdot \nabla v_h \, d\mu = \int_\Omega u_h v_h \, d\mu
    = \int_\Omega j_h v_h \, dx ,
```

i.e. ``\Lambda = (\mathbb{K}^\mu)^{-1}\mathbb{M}`` with the **plain** mass matrix on the right —
the ``1/r`` of the measure cancels against the ``r`` of the substitution. Then
``\Lambda^T \mathbb{K}^\mu = \mathbb{M}``, so

```math
\frac{\partial H}{\partial \hat{j}}
  = \Lambda^T \mathbb{K}^\mu \Lambda \hat{j} = \mathbb{M} \hat\psi ,
```

and the degeneracy is structural. The energy is still the poloidal magnetic energy,
``H = \tfrac12 \hat\psi^T \mathbb{K}^\mu \hat\psi = \tfrac12 \int_\Omega |\nabla\psi_h|^2 d\mu``
— nothing has been redefined to make the algebra work.

The price is stated rather than hidden: ``u_h = r j_h`` is a spline times ``r`` and not itself a
spline, so the trial space is ``r V`` where the manuscript's is ``V``. That is a different
discretisation of the same problem, of the same order, and it is the deviation this file
records.

# The two state spaces

``\psi`` is the homogeneous-Dirichlet ``\Delta^*`` solve either way; `state` says where ``j``
lives, and the answer matters for the same reason it did in §5.4 ([`SECTION5_RUNS`](@ref)).
`eq:gs-ref` is the ``\mu = 0``, ``c = 0`` member of the equilibrium family
``\delta S/\delta j = \lambda\psi + c\cdot x + \mu``, whose extra multipliers belong to the
bracket's mass and momentum Casimirs. Only a space that does not contain ``1``, ``x_1``,
``x_2`` forces them to vanish, so `:dirichlet` — ``j \in V_D`` — is what makes the manuscript's
reference exact, and the price is that ``\int_\Omega u \, d\mu = \int_\Omega j \, dx`` drifts.
`:free` is implemented, reaching the elliptic solve through the recombination matrix exactly as
`EulerSquare` does, and `verify_gradshafranov.jl` relaxes both spaces as a **control**: after
240 steps the free space's state misses `eq:gs-ref` by ``2.52 \times 10^{-1}`` against the
Dirichlet space's ``4.28 \times 10^{-4}`` — a factor 589 — and carries the three multipliers
``(\mu, c_1, c_2)`` at 2456 times the size, while its mass is conserved to
``4 \times 10^{-16}`` where the Dirichlet space's drifts 83 %. That is the Casimirs being
present or absent, seen directly. It is a control rather than an alternative formulation, so no
run uses it.

**The control has to be a relaxation and not an eigenvalue problem**, which is worth stating
because the cheaper version reports nothing: the pencil ``(\mathbb{M}\Lambda, \mathbb{W})``
carries no mass constraint, so its lowest eigenvector is the ``\mu = 0``, ``c = 0`` member in
*both* spaces — measured, their eigenvalues agree to ``5 \times 10^{-9}`` and both eigenvectors
fit `eq:gs-ref` to ``4 \times 10^{-4}``. What separates the spaces is that the **flow** in the
free space conserves the mass and the momenta and therefore cannot reach that member from an
initial state whose mass is nonzero.

# The mesh is capped by the Jacobian, not by accuracy

The manuscript uses ``64 \times 64``. This reproduction does not, for the reason `EulerSquare`
measures: the collision bracket is nonlocal, its operator is dense, and one Newton matrix is
``N`` dense assemblies of size ``N^2``. What resolution is needed for is the initial Gaussian's
narrower direction, ``w_1 = 0.707`` on a domain of width 6 — **not** the relaxed state, which is
the lowest ``\Delta^*`` eigenmode and is already resolved on any mesh these runs reach:
[`gs_eigenvalue`](@ref) is within ``4.9 \times 10^{-9}`` of the continuum eigenvalue at degree 3
on ``12 \times 14`` cells, and within ``6.4 \times 10^{-5}`` at degree 2 on ``10 \times 12``.
Measured, the ``L^2`` error of the projected Gaussian runs ``6.7\times10^{-2}``,
``1.5\times10^{-2}``, ``5.1\times10^{-3}``, ``2.3\times10^{-3}`` at ``10\times12``,
``14\times16``, ``18\times21``, ``22\times26``, falling at the space's own order.
See [`SECTION55_RUNS`](@ref) for the numbers chosen.
"""
struct GradShafranovBox{T, ST <: TensorSplineSpace{T, 2}, MT}
    space::ST
    Λ::Matrix{T}
    M::MT
    MΛ::Matrix{T}
    W::Matrix{T}
    b::Vector{T}
    state::Symbol
end

function GradShafranovBox(cells::Tuple{Int, Int}, degree::Int = 2;
        state::Symbol = :dirichlet)
    state in (:dirichlet, :free) || throw(ArgumentError(
        "the state space is :dirichlet or :free, but got :$(state)"))
    meshes = (UniformMesh(cells[1], GS_RADIAL), UniformMesh(cells[2], GS_AXIAL))
    bc = state === :dirichlet ? Dirichlet() : Free()
    s = TensorSplineSpace(meshes, degree, bc)

    M = sparse(mass_matrix(s))
    Kμ = sparse(gs_stiffness(s))
    Λ = if state === :dirichlet
        cholesky(Symmetric(Kμ)) \ Matrix(M)
    else
        # The embedding V_D ⊂ V, per axis, then tensored in REVERSE axis order because the flat
        # index runs the first axis fastest — the same constraint `EulerSquare` documents, and
        # here the two axes have different cell counts, so getting it wrong is a dimension
        # mismatch rather than a silent transpose.
        R = map(k -> recombination_matrix(BSplineBasis(meshes[k], degree, Dirichlet())), (
            1, 2))
        E = sparse(kron(R[2], R[1]))
        E * (cholesky(Symmetric(Matrix(E' * Kμ * E))) \ Matrix(E' * M))
    end

    # MΛ = M (K^μ)⁻¹ M is symmetric exactly; the factorisation makes it so only to round-off,
    # and `QuadraticHamiltonian` checks symmetry rather than assuming it.
    A = M * Λ
    σ = [gs_entropy_weight(x[1]) for x in quadrature_nodes(s)]
    W = Matrix(weighted_matrix(s, σ, (0, 0), (0, 0)))

    GradShafranovBox(s, Matrix(Λ), M, Matrix((A .+ A') ./ 2), (W .+ W') ./ 2,
        Vector{Float64}(basis_integrals(s)), state)
end

function GradShafranovBox(cells::Int, degree::Int = 2; kwargs...)
    GradShafranovBox((cells, cells), degree; kwargs...)
end

space(box::GradShafranovBox) = box.space
nbasis(box::GradShafranovBox) = nbasis(box.space)

function Base.show(io::IO, box::GradShafranovBox)
    print(io, "GradShafranovBox(cells=", ncells(box.space), ", p=", degree(box.space)[1],
        ", N=", nbasis(box.space), ", state=:", box.state, ")")
end

@doc raw"""
    integrate(box::GradShafranovBox, ĵ)

``\int_\Omega j_h \, dx = \int_\Omega u_h \, d\mu``, the mass Casimir of the collision
bracket.

The two integrals are equal because ``u = rj`` and ``d\mu = r^{-1}dr\,dz``, which is the same
cancellation the state variable of [`GradShafranovBox`](@ref) is chosen for. It is a Casimir
because its derivative ``\partial C/\partial\hat{j} = b`` satisfies
``\mathbb{M}^{-1} b = \mathbf{1}``, whose gradient vanishes, so ``Q_2`` annihilates it —
and it is conserved discretely only in the `:free` space, where the constant function is
present.
"""
integrate(box::GradShafranovBox, ĵ::AbstractVector) = dot(box.b, ĵ)

"``(j_h, v_h)_{L^2(dx)}``, through the plain mass matrix."
l2inner(box::GradShafranovBox, ĵ::AbstractVector, v̂::AbstractVector) = dot(ĵ, box.M, v̂)

l2norm(box::GradShafranovBox, ĵ::AbstractVector) = sqrt(max(l2inner(box, ĵ, ĵ), 0.0))

@doc raw"""
    gs_state(box, spec)

The initial state of run `spec`, as the ``L^2`` projection of ``j_0 = u_0/r`` onto the state
space of `box`.

The manuscript's initial condition is the Gaussian ``u_0``; the state variable is ``j = u/r``
([`GradShafranovBox`](@ref)), so the division happens here and exactly once.
"""
function gs_state(box::GradShafranovBox, spec::GSSpec)
    u₀ = initial_condition(spec)
    project(box.space, x -> u₀(x[1], x[2]) / x[1])
end

@doc raw"""
    gs_flow(box)

The [`MetriplecticFlow`](@ref) of the Section 5.5 relaxation: the collision-like bracket
generated by the ``\Delta^*`` energy, that energy, and the Herrnegger-Maschke entropy.

```math
H = \tfrac12 \hat{j}^T \mathbb{M}\Lambda \hat{j} , \qquad
S = \tfrac12 \hat{j}^T \mathbb{W} \hat{j} , \qquad
\mathbb{G} = \mathbb{M}^{-1}\mathbb{A}(\hat{j})\mathbb{M}^{-1} ,
```

with the bracket's `density` carrying ``d\mu = r^{-1}dr\,dz`` into both the inner moment
quadrature and the outer ``\operatorname{div}_\mu`` assembly, and its `mobility`
``M = Cr^2+D`` fixed by `eq:M-condition`.

`mobility_derivative = 0` is not a shortcut — ``M`` depends on ``r`` alone. That is what makes
``\mathbb{G}`` quadratic in ``\hat{j}``, ``\partial S/\partial\hat{j}`` linear, and hence the
whole dissipative residual an exact **cubic polynomial** in the degrees of freedom, with an
analytic Jacobian. [`CollisionBracket`](@ref) refuses to difference a mobility it was not given
the derivative of, for exactly this reason.

The Poisson half is absent, which is the case every §5 experiment runs. The energy and the
bracket are generated by the *same* ``\Lambda``, which is what makes the degeneracy
``\mathbb{G}\,\partial H/\partial\hat{j} = 0`` structural rather than accidental — see
[`GradShafranovBox`](@ref) for why that needs the ``j`` variable.
"""
function gs_flow(box::GradShafranovBox)
    G = CollisionBracket(box.space, box.Λ;
        mobility = (x, u) -> herrnegger_mobility(x[1]),
        mobility_derivative = 0,
        density = gs_density)
    MetriplecticFlow(box.space, G, QuadraticHamiltonian(box.MΛ),
        QuadraticHamiltonian(box.W))
end

## The closed-form reference

@doc raw"""
    gs_current(box, ĵ)

The manuscript's state variable ``u_h = r \, j_h``, sampled on the quadrature grid.
"""
function gs_current(box::GradShafranovBox, ĵ::AbstractVector)
    [x[1] for x in quadrature_nodes(box.space)] .* field(box.space, ĵ, (0, 0))
end

@doc raw"""
    gs_ordinate(box, ĵ)

The field ``u_h/(Cr^2+D) = \sigma(r) j_h`` on the quadrature grid — the ordinate of the
manuscript's scatter plot and the left-hand side of `eq:gs-ref`, and also
``\delta S/\delta j``.
"""
function gs_ordinate(box::GradShafranovBox, ĵ::AbstractVector)
    [gs_entropy_weight(x[1]) for x in quadrature_nodes(box.space)] .*
    field(box.space, ĵ, (0, 0))
end

@doc raw"""
    gs_fit(box, ĵ)

How close `ĵ` is to satisfying `eq:gs-ref`, ``u/(Cr^2+D) = \lambda\psi``, as
`(λ, residual, relative)`.

The fit is in ``L^2(\mu)``, the measure the variational principle is posed in:
``\lambda = \langle \sigma j, \psi\rangle_\mu / \langle \psi,\psi\rangle_\mu``, `residual` is
``\|\sigma j - \lambda\psi\|_{L^2(\mu)}`` and `relative` that over
``\|\sigma j\|_{L^2(\mu)}``.

The two numbers are a **two-sided** check and neither stands alone, exactly as for §5.4's
[`eigenmode_fit`](@ref). ``\lambda`` is a projection coefficient, so its error is *second*
order in `relative`: it is already close for states that are nothing like an equilibrium, and
asserting on it alone would pass a run that had barely moved. `relative` is the number that
says the scatter has collapsed onto a line.

[`gs_rayleigh`](@ref) is a third estimate of the same ``\lambda`` from the invariants rather
than from the fields, and the two agreeing is a statement about the state.
"""
function gs_fit(box::GradShafranovBox, ĵ::AbstractVector)
    s = box.space
    wμ = quadrature_weights(s) .* gs_density.(quadrature_nodes(s))
    y = gs_ordinate(box, ĵ)
    ψ = field(s, box.Λ * ĵ, (0, 0))
    λ = dot(wμ, y .* ψ) / dot(wμ, ψ .* ψ)
    res = sqrt(max(dot(wμ, (y .- λ .* ψ) .^ 2), 0.0))
    return (λ, res, res / sqrt(dot(wμ, y .^ 2)))
end

@doc raw"""
    gs_rayleigh(box, ĵ)

The Rayleigh quotient ``S(\hat{j}) / H(\hat{j})``, which is the Grad-Shafranov eigenvalue at
an equilibrium and an upper bound on it everywhere.

At a constrained critical point of `eq:VP-GS` the multiplier satisfies
``\delta S/\delta j = \lambda \, \delta H/\delta j``, hence

```math
S = \tfrac12 \int_\Omega j \, \frac{\delta S}{\delta j} \, dx
  = \tfrac{\lambda}{2} \int_\Omega j \, \psi \, dx = \lambda H ,
```

both functionals being quadratic. Minimising ``S/H`` reproduces
`eq:Grad-Shafranov-equation` — the stationarity condition is
``\sigma j = \lambda \psi``, and eliminating ``j`` gives
``-\Delta^*\psi = \lambda(Cr^2+D)\psi`` — so

```math
\frac{S(\hat{j})}{H(\hat{j})} \ge \lambda_{\min}
```

for **every** admissible state and not only along the flow. That makes it a precondition to
assert at ``t = 0`` rather than a result to check at ``t = T``: it is the Grad-Shafranov
counterpart of §5.4's Poincaré floor [`euler_entropy_floor`](@ref), read as a quotient so that
no second entropy-floor function is needed.

Its error at a nearly-relaxed state is second order in the residual of [`gs_fit`](@ref), which
makes it the sharpest of the three ``\lambda`` estimates and the one whose agreement with
[`takeda_iterate`](@ref) is the reproduction's actual claim.
"""
function gs_rayleigh(box::GradShafranovBox, ĵ::AbstractVector)
    dot(ĵ, box.W, ĵ) / dot(ĵ, box.MΛ, ĵ)
end

@doc raw"""
    gs_eigenvalue(box::GradShafranovBox)

The Grad-Shafranov eigenvalue of the space `box` is built on — the number a run on it can
actually reach, and what [`gs_rayleigh`](@ref) is bounded below by.

Read off the assembled invariants rather than from ``\mathbb{K}^\mu``, and **in the reciprocal
direction**, so that one method serves both state spaces. `eq:gs-ref` tested against the basis
is ``\mathbb{W}\hat{j} = \lambda \, \mathbb{M}\Lambda\hat{j}``, i.e. the stationary values of
[`gs_rayleigh`](@ref); but in the `:free` space ``\mathbb{M}\Lambda`` has rank
``\dim V_D < \dim V`` — the elliptic solve lands in ``V_D`` — so that pencil has a singular
right-hand matrix and no Cholesky. ``\mathbb{W}`` is positive definite in either space, since
``\sigma(r) = r/(Cr^2+D) > 0`` on ``\overline\Omega``, so the pencil is taken the other way
round and inverted:

```math
\lambda_h = 1 \big/ \lambda_{\max}(\mathbb{M}\Lambda, \mathbb{W}) .
```

[`gs_eigenvalue`](@ref)`(space)` computes the same number from ``\mathbb{K}^\mu`` and
``\mathbb{B}`` directly and shares no matrix with this; `verify_gradshafranov.jl` checks that
they agree, which is what says the energy and entropy forms assembled here are the ones the
eigenvalue problem is about.
"""
function gs_eigenvalue(box::GradShafranovBox)
    inv(maximum(real, eigvals(Symmetric(box.MΛ), Symmetric(box.W))))
end
