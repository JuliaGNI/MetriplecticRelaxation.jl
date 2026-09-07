#
# Section 5.4 of *Metriplectic relaxation to equilibria*: reduced Euler on the unit square.
#
# The problem definition and the discretisation together, unlike Section 4, where `torus.jl`
# holds the problem and `spectral.jl`/`spline.jl` two independent solvers of it. There is only
# one solver here -- the collision-like bracket has no spectral counterpart in this tree -- so a
# second file would separate nothing.
#
# WHAT CHANGES FROM SECTION 4, and none of it is cosmetic:
#
#   * the domain is [0,1]^2 with HOMOGENEOUS DIRICHLET boundary conditions rather than the
#     periodic torus, so the stiffness matrix is nonsingular, the bordered Poisson solve of
#     `spline.jl` is unnecessary, and CONSTANTS ARE NOT IN THE SPACE.  That last point is what
#     makes both of the manuscript's closed-form references exact -- see `SECTION5_RUNS`;
#   * the bracket is `CollisionBracket`, which is NONLOCAL: its operator is dense, and its
#     analytic Jacobian costs N dense assemblies.  That is what CAPS the mesh; what fixes it at
#     26 is resolving the initial condition's narrow direction.  Accuracy of the relaxed state
#     constrains neither -- see `EulerSquare`;
#   * the integrator is `ImplicitMidpoint`, which for this field IS Crank-Nicolson, so every
#     step is a Newton solve rather than four function evaluations;
#   * B3's entropy is `s = y log y` rather than `y^2/2`, which needs `omega > 0` pointwise and
#     is the case the manuscript warns needs small time steps.
#

@doc raw"""
The side length of the Section 5.4 domain ``\Omega = [0,1]^2``.
"""
const SQUARE_LENGTH = 1.0

@doc raw"""
The first Dirichlet eigenvalue of ``-\Delta`` on ``[0,1]^2``,
``\lambda_{1,1} = 2\pi^2 \approx 19.7392``, with eigenfunction
``\sin(\pi x_1) \sin(\pi x_2)``.

This is the reference B1 relaxes to: the constrained entropy minimiser of
`eq:entropy-2D` with ``s = y^2/2`` at fixed energy is ``\omega = \lambda_{1,1} \phi``, the
lowest eigenmode, and its entropy is ``S_\eta = \lambda_{1,1} H_0`` — see
[`euler_entropy_floor`](@ref).
"""
const DIRICHLET_EIGENVALUE = 2π^2

## The initial condition

@doc raw"""
    gaussian_w2(x₀, w², N)

The Section 5 Gaussian, built from the **squared** widths the manuscript prints.

```math
\omega_0(x) = \frac{1}{N} \exp\Big[
    - \frac{(x_1 - x_{0,1})^2}{w_1^2} - \frac{(x_2 - x_{0,2})^2}{w_2^2} \Big]
```

Section 4 states ``w``; **Section 5 states ``w^2``** — `w₁² = 0.01`, `w₂² = 0.07` for every B
run — and [`Gaussian`](@ref) takes ``w``. This function is the *only* place the square root is
taken, so that the conversion exists once and can be tested once: `verify_euler.jl` asserts that
the e-folding distance along each axis is ``\sqrt{w_k^2}`` and not ``w_k^2``, which is the
transcription error this whole indirection exists to make impossible.

`images = 0`, i.e. the formula exactly as printed. The periodic-image reading of
[`Gaussian`](@ref) is a Section 4 question about the torus and has no meaning on a bounded
square.
"""
gaussian_w2(x₀, w², N) = Gaussian(x₀, map(sqrt, w²), N)

@doc raw"""
    perturbation_b2(x₁, x₂)

B2's added mode ``\sin(6\pi x_1) \sin(4\pi x_2)``.

It vanishes on ``\partial\Omega`` and is an eigenfunction of ``-\Delta`` there, with eigenvalue
``(6\pi)^2 + (4\pi)^2 = 52\pi^2``. An eigenmode satisfies ``\omega = \lambda \phi`` exactly,
which by the equilibrium condition of [`euler_flow`](@ref) makes it a **stationary state** of
the relaxation — and, not being the lowest one, an unstable one.

That is what B2 is: an equilibrium plus a perturbation. Its Gaussian carries ``N = 100``, so
the bump is a 1 % disturbance of a mode of amplitude one, and the entropy therefore sits still
until the instability has grown — the plateau-then-decay signature the run checks for.
"""
perturbation_b2(x₁, x₂) = sin(6π * x₁) * sin(4π * x₂)

@doc raw"""
The positive background B3's initial condition is given, as an absolute vorticity offset.

**This is a departure from the printed initial condition, and it is forced.** ``s = y \log y``
is undefined at ``y \le 0`` and `eq:M-condition` gives it the mobility ``M = y``, which the same
equation requires to be positive — so B3's state must be **strictly** positive, everywhere, at
every step. The printed Gaussian is strictly positive as a *function*: with ``w_1^2 = 0.01`` it
decays to ``1.4 \times 10^{-11}`` of its peak at the far corner of the square, which is not
zero. It is, however, indistinguishable from zero for a Galerkin scheme, which is not
positivity-preserving: measured on the ``26^2`` space the runs use, the projected state starts
at ``\min\omega_0 = 6.0 \times 10^{-13}`` — the Gaussian's own value at that node — and the
**first** step at ``\Delta t = 0.02`` drives the iterate to ``-6.6 \times 10^{-5}``, where the
entropy raises a `DomainError`. No step completes: the margin is too thin to survive one.
`run_b3.jl` performs exactly that, as the control on this constant.

The continuum does not have this problem, and the reason is worth stating because it says what
the discretisation is losing: the mobility ``M = \omega`` *vanishes* where ``\omega`` does, so
the flux vanishes with it and the equation is degenerate-parabolic and positivity-preserving. A
Galerkin projection of it is not.

A floor is the smallest change that restores admissibility. At **1 % of the peak** it is four
orders above the undershoot it has to absorb, and on the homogeneous-Dirichlet space it leaves
``\min\omega_0 = 9.3 \times 10^{-3}`` — a real margin where the printed condition leaves
``6 \times 10^{-13}``. It costs 12 % of the mass, which §5.4 measures nothing against: the
reference ``\omega = e^{\lambda\phi-1}`` and its multiplier are both read off the *relaxed*
state, not off the initial data. The alternative — a positivity-preserving limiter, or evolving
``\log\omega`` — is a different scheme, and this reproduction does not write one.
"""
const B3_FLOOR = 0.1

## The runs

@doc raw"""
    EulerSpec

One of the three Section 5.4 runs: which entropy it dissipates, its initial condition, and the
time step and final time this reproduction chose for it.

`entropy` is `:quadratic` for ``s = y^2/2`` or `:gibbs` for ``s = y \log y``; it selects both
the entropy functional and, through `eq:M-condition` ``M \, \partial_y^2 s = 1``, the mobility
of the bracket. See [`euler_flow`](@ref).

`state` is which space the **vorticity** lives in, `:dirichlet` or `:free`. ``\phi`` is always
the homogeneous-Dirichlet solve; what differs is whether ``\omega`` is constrained to vanish on
``\partial\Omega`` too. It is not a free parameter and it is not the same for all three runs —
see [`EulerSquare`](@ref) and [`SECTION5_RUNS`](@ref).

Unlike Section 4's [`RunSpec`](@ref), **neither `Δt` nor `T` is the manuscript's**: §5 states no
time step, no final time and no stopping criterion. Both are recorded choices and the reasoning
is in [`SECTION5_RUNS`](@ref) and in `CHANGELOG.md`, beside the numbers they produced.
"""
struct EulerSpec{B}
    name::String
    section::String
    entropy::Symbol
    state::Symbol
    gaussian::Gaussian{Float64}
    background::B
    Δt::Float64
    T::Float64
end

"""
    initial_condition(spec::EulerSpec)

The initial condition of run `spec` as a callable `(x₁, x₂)`: the Gaussian, plus the
`background` mode where the run has one.
"""
function initial_condition(spec::EulerSpec)
    spec.background === nothing ? spec.gaussian :
    (x₁, x₂) -> spec.background(x₁, x₂) + spec.gaussian(x₁, x₂)
end

@doc raw"""
The three runs of Section 5.4, keyed by name.

| run | figs | ``s(y)`` | ``M`` | initial condition |
|:--|:--|:--|:--|:--|
| `b1` | `sv_*` | ``y^2/2`` | ``1`` | ``\omega_G``, ``x_0 = (\tfrac12,\tfrac12)``, ``w_1^2 = 0.01``, ``w_2^2 = 0.07``, ``N = 1`` |
| `b2` | `pe_*` | ``y^2/2`` | ``1`` | ``\sin(6\pi x_1)\sin(4\pi x_2) + \omega_G``, ``N = 100`` |
| `b3` | `ge_*` | ``y \log y`` | ``y`` | ``\omega_G``, ``1/N = 10``, **plus** [`B3_FLOOR`](@ref) |

Everything in that table is the manuscript's except B3's floor, which is forced and is explained
where it is defined. **Nothing in the two columns below is the manuscript's either**, because
§5 states no time step, no final time, no stopping criterion, no nonlinear-solver description
and no tolerance. Each is a choice of this reproduction, recorded here with what it produced.

| run | ``\Delta t`` | ``T`` | steps |
|:--|--:|--:|--:|
| `b1` | 1.0 | 200 | 200 |
| `b2` | 0.5 | 150 | 300 |
| `b3` | 0.02 | 3 | 150 |

**``\Delta t = 1`` is not a typo, and it is not a licence Crank-Nicolson's A-stability alone
would give.** It was measured: on B1 the whole entropy trace at ``\Delta t = 1`` agrees with the
trace at ``\Delta t = 0.1`` to five digits, and the energy is conserved at ``10^{-14}`` in both.
The relaxation has no fast scale to resolve — the flow is a nonlinear diffusion whose slowest
mode sets a decay rate near ``0.09`` per unit time — so accuracy does not constrain the step,
and the cost does: one step is one Newton solve, and one Newton matrix is ``N`` dense
assemblies of the nonlocal bracket. Each driver re-measures the step-halving difference rather
than trusting this paragraph.

**``T`` is set by the stopping criterion, which is also a choice: the excess entropy.** B1's
entropy approaches ``S_\eta = \lambda_{1,1} H_0`` at a rate near ``0.09``, so ``T = 200``
leaves it ``10^{-8}`` of the way above the floor and the state ``10^{-4}`` from
``\lambda_{1,1}\phi`` in ``L^2``; ``T = 30`` would leave 3 % and 17 %, which draws the figure
and cannot settle the reference.

B2 spends **half** the step on a **shorter** horizon, and both halves of that are the plateau's
doing. What its claim needs resolved is the transition: measured, its entropy production grows
by nearly three orders of magnitude between ``t = 0`` and ``t \approx 18``, and at
``\Delta t = 1`` that span is nine samples rather than eighteen. The horizon then has to come
back down to keep the step count affordable, and it can: B2's own late relaxation is about four
times slower than B1's, so ``T = 150`` leaves the state 2 % from ``\lambda_{1,1}\phi`` where
``T = 200`` would leave 1 % — a difference that changes no claim, against an hour of wall clock
that it does.

B3 is a different regime again and its ``\Delta t`` is the one number in this table that is
*measured* rather than chosen, because for ``s = y \log y`` Crank-Nicolson no longer guarantees
monotone dissipation and the manuscript says only that "sufficiently small time steps must be
used". The step-size study `run_b3.jl` performs locates the boundary, and it is **admissibility
rather than accuracy** that sets it: measured at 12 cells, ``\Delta t`` up to ``0.08`` is
monotone and ``0.16`` leaves the admissible set within ten steps; at 16 cells the boundary sits
between ``0.05`` and ``0.1``, and ``\Delta t = 0.05`` reproduces ``0.02``'s entropy at
``t = 0.5`` to six digits. It tightens with the mesh, so the run takes ``0.02`` — a factor of
four inside the coarse-mesh boundary and re-measured on the run's own mesh by the sweep. B3 also
relaxes two orders of magnitude faster than B1, its entropy having fallen 38 % by ``t = 0.5``,
so ``T = 3`` is six relaxation times.

**All three step counts are set by the same wall clock.** One step is one Newton solve and one
Newton matrix is ``N`` dense assemblies, so at ``N = 676`` a step costs **27 s** measured on
this machine at ``\Delta t = 1`` — the Newton iteration is what varies with the step size, and
a larger step buys fewer but more expensive steps rather than a proportional saving.

# The homogeneous-Dirichlet space is what makes both references exact

Not a boundary condition inherited from the geometry — a *modelling* choice, and both of §5.4's
closed-form limits depend on it. ``\phi`` is the Dirichlet solve either way; the question is
whether ``\omega`` is constrained too.

At equilibrium the flow satisfies ``\nabla(\delta S/\delta u) = \lambda \nabla(\delta H/\delta
u)``, hence

```math
\frac{\delta S}{\delta u} = \lambda \phi + c \cdot x + \mu
```

with ``\mu`` and ``c`` the multipliers of the bracket's other Casimirs, the mass
``\int_\Omega u`` and the momenta ``\int_\Omega x_k u`` — each a Casimir because its variational
derivative has vanishing gradient, so ``Q_2`` annihilates it identically. Both of the
manuscript's references are the ``\mu = 0``, ``c = 0`` case: ``\omega = \lambda_{1,1}\phi`` for
``s = y^2/2`` and ``\omega = e^{\lambda\phi - 1}`` for ``s = y\log y``, the second of which
gives ``\lambda = (M+S)/2H_0`` by one substitution.

| | ``\omega \in V_D`` (`:dirichlet`) | ``\omega \in V`` (`:free`) |
|:--|:--|:--|
| ``1``, ``x_1``, ``x_2`` in the space | no | yes |
| mass and momenta conserved | no | yes |
| ``\mu`` and ``c`` | forced to zero | fixed by ``M_0`` and ``P_0`` |
| the manuscript's references | **exact** | carry an extra ``\mu`` |

All three runs therefore use ``V_D``, and the price is that the mass drifts — B1's by 46 % over
its run — which the drivers report rather than assert on.

**The `:free` alternative is implemented, and B3 runs it as a control rather than as its
formulation.** Both spaces are viable for B3 once [`B3_FLOOR`](@ref) makes the state admissible,
and they disagree exactly as the table predicts: measured at 16 cells over the same 25 steps,
``V_D`` has ``\mu`` shrinking through ``-1.24, -0.85, -0.63, -0.47, -0.35`` toward zero while
``V`` holds it at ``-1.46``, and the residual against the manuscript's ``\mu = 0`` reference
falls to 0.147 in ``V_D`` and sticks at 0.41 in ``V``. That is the mass Casimir being present or
absent, seen directly, and it is why `run_b3.jl` measures both.
"""
const SECTION5_RUNS = Dict(
    "b1" => EulerSpec("b1", "5.4", :quadratic, :dirichlet,
        gaussian_w2((0.5, 0.5), (0.01, 0.07), 1.0), nothing, 1.0, 200.0),
    "b2" => EulerSpec("b2", "5.4", :quadratic, :dirichlet,
        gaussian_w2((0.5, 0.5), (0.01, 0.07), 100.0), perturbation_b2, 0.5, 150.0),
    "b3" => EulerSpec("b3", "5.4", :gibbs, :dirichlet,
        gaussian_w2((0.5, 0.5), (0.01, 0.07), 0.1), (x₁, x₂) -> B3_FLOOR, 0.02, 3.0))

"The three runs in the order Section 5.4 presents them."
const SECTION5_ORDER = ("b1", "b2", "b3")

## The Gibbs entropy

@doc raw"""
    GibbsEntropy()

The entropy ``S(\omega) = \int_\Omega \omega \log \omega \, dx`` of B3, `eq:entropy-2D` with
``s(y) = y \log y``.

```math
\frac{\partial S}{\partial \hat\omega_I} = \int_\Omega (\log \omega_h + 1) \, \Phi_I \, dx ,
\qquad
\frac{\partial^2 S}{\partial \hat\omega_I \partial \hat\omega_J} =
    \int_\Omega \frac{\Phi_I \Phi_J}{\omega_h} \, dx ,
```

both by quadrature on the space's own grid, both analytic. PoissonBrackets has no nonlinear
`DiscreteHamiltonian` — `QuadraticHamiltonian` and `MassCasimir` are the two it carries — and
this one is specific to §5.4, so it lives here. A candidate for the package once a second
caller turns up, as `LinearHamiltonian` and `EllipticEnergy` in `spline.jl` are.

!!! warning "``\omega > 0`` is a requirement of the problem, not of the discretisation"
    ``y \log y`` is undefined for ``y \le 0`` and its second derivative ``1/y`` — which is the
    mobility's reciprocal through `eq:M-condition` — is singular at zero. Nothing here guards
    against it: `log` of a negative number raises, which is the correct outcome, and a guard
    would convert a state the entropy does not admit into a silently wrong number.

    That makes the initial state a real constraint rather than a formality. On a
    homogeneous-Dirichlet space the ``L^2`` projection of B3's Gaussian oscillates next to the
    peak, and the undershoot is **negative** until the narrow direction is resolved: measured at
    degree 2, the minimum over the quadrature grid is ``-5.8 \times 10^{-4}`` at 16 cells,
    ``-3.2 \times 10^{-9}`` at 24, and positive from 26 on. That measurement, not the accuracy
    of the relaxed state, is what sets B3's mesh; `verify_euler.jl` asserts it and `run_b3.jl`
    re-checks positivity at every recorded sample.
"""
struct GibbsEntropy{T} <: DiscreteHamiltonian{T} end

GibbsEntropy() = GibbsEntropy{Float64}()

function hamiltonian(::GibbsEntropy, s::TensorSplineSpace, ω̂::AbstractVector)
    u = field(s, ω̂, (0, 0))
    dot(quadrature_weights(s), u .* log.(u))
end

function gradient(::GibbsEntropy, s::TensorSplineSpace, ω̂::AbstractVector)
    u = field(s, ω̂, (0, 0))
    basis_values(s, (0, 0)) * (quadrature_weights(s) .* (log.(u) .+ 1))
end

function hessian(::GibbsEntropy, s::TensorSplineSpace, ω̂::AbstractVector)
    u = field(s, ω̂, (0, 0))
    Φ = basis_values(s, (0, 0))
    Matrix(Φ * Diagonal(quadrature_weights(s) ./ u) * Φ')
end

## The discretisation

@doc raw"""
    EulerSquare(n, p = 2; state = :dirichlet)

A tensor-product B-spline space of `n` cells and degree `p` per direction on
``\Omega = [0,1]^2``, with the operators the Section 5.4 runs need assembled once.

`Λ` is the **assembled dense** solution operator of the Dirichlet Poisson problem
``-\Delta\phi = \omega``, ``\phi|_{\partial\Omega} = 0``; `M` the sparse mass matrix, `MΛ` the
matrix of the energy, and `b` the basis integrals, so that ``\int_\Omega \omega_h = b \cdot
\hat\omega``.

# The two state spaces

``\phi`` is the homogeneous-Dirichlet solve either way. `state` says where ``\omega`` lives, and
[`SECTION5_RUNS`](@ref) says why the answer is not the same for all three runs.

  - `:dirichlet` — ``\omega \in V_D``, the recombined homogeneous-Dirichlet basis, ``(n+p-2)^2``
    degrees of freedom, and ``\Lambda = \mathbb{K}^{-1}\mathbb{M}`` with ``\mathbb{K}``
    nonsingular because the constants are gone.
  - `:free` — ``\omega \in V``, the plain clamped basis, ``(n+p)^2`` degrees of freedom. The
    elliptic solve still happens in ``V_D \subset V``, reached through the **recombination
    matrix** ``\mathbb{E}``, which is exactly the matrix expressing each Dirichlet basis
    function in the clamped one, so

    ```math
    \Lambda = \mathbb{E} \, (\mathbb{E}^T \mathbb{K} \mathbb{E})^{-1} \mathbb{E}^T \mathbb{M} .
    ```

    ``\mathbb{M}\Lambda = \mathbb{M}\mathbb{E}\mathbb{K}_D^{-1}\mathbb{E}^T\mathbb{M}`` is
    symmetric positive semi-definite exactly, so the energy is still
    ``\tfrac12 \int|\nabla\phi|^2``, and ``\phi`` still vanishes on ``\partial\Omega`` — while
    ``\omega`` no longer has to.

# Why `Λ` is assembled, where `spline.jl` keeps a factorisation instead

Opposite choice, opposite reason. Section 4's [`PoissonMap`](@ref) is applied once per
Runge-Kutta stage and never differentiated, so a cached sparse factorisation beats a dense
matrix by the whole memory bandwidth. Here [`CollisionBracket`](@ref)'s analytic
[`metric_derivative`](@ref) reads ``\Lambda`` **column by column** — ``\partial\hat\phi /
\partial\hat\omega_m`` is `Λ[:, m]` — so the columns have to exist. Assembling costs one sparse
Cholesky of ``\mathbb{K}`` and ``N`` back-substitutions, once, before the time loop; storage is
``N^2``, which at the resolutions below is tens of megabytes.

That also settles the caching question the elliptic solve otherwise raises: ``\mathbb{K}`` is
fixed in time, and factorising it inside a residual would cost more than the entire moment
collapse of the bracket saves. It is factorised exactly once, here, and the factorisation is
discarded with the constructor.

The Dirichlet stiffness matrix is nonsingular, unlike the periodic case: constants are not in
``V_D``, so there is no kernel to border out as `spline.jl` has to.

# The mesh is capped by the Jacobian and fixed by the initial condition, not by accuracy

The manuscript uses ``64^2`` P2 elements. This reproduction does not, and the reason is
measurable rather than a matter of taste: the collision bracket is nonlocal, so its operator is
dense, and `metric_directional` is ``N`` dense assemblies of size ``N^2`` — one Newton matrix is
``O(N^3)``. Measured on this machine at degree 2, one implicit-midpoint step takes **1.0 s** at
``N = 256``, **7.7 s** at ``N = 576`` and **36 s** at ``N = 1024``; the manuscript's ``64^2``
would be ``N = 4096`` and some 40 minutes *per step*.

None of the §5.4 claims needs that resolution. The relaxed states are the lowest Dirichlet
eigenmode and a smooth exponential of it, and the discrete first eigenvalue is already within
``3.4 \times 10^{-5}`` of ``2\pi^2`` at 8 cells.

What *does* need resolution is the **narrow direction of the initial condition**, ``w_1 = 0.1``,
and that is what fixes the number at 26: it is the coarsest mesh on which the ``L^2`` projection
of the printed Gaussian stops oscillating below zero. Measured at degree 2, ``\min\omega_0`` runs
``-1.4 \times 10^{-2}``, ``-5.8 \times 10^{-4}``, ``-4.2 \times 10^{-7}``,
``+6.0 \times 10^{-13}`` at 12, 16, 20 and 26 cells — a resolution statement about the peak, not
about the boundary, and unchanged by refining the other axis. `verify_euler.jl` asserts it with
the three coarser rows as the control that 26 is a threshold rather than a preference.
"""
struct EulerSquare{T, ST <: TensorSplineSpace{T, 2}, MT}
    space::ST
    Λ::Matrix{T}
    M::MT
    MΛ::Matrix{T}
    b::Vector{T}
    state::Symbol
end

function EulerSquare(n::Int, p::Int = 2; state::Symbol = :dirichlet)
    state in (:dirichlet, :free) || throw(ArgumentError(
        "the state space is :dirichlet or :free, but got :$(state)"))
    s = TensorSplineSpace((n, n), p,
        state === :dirichlet ? Dirichlet() : Free(); L = SQUARE_LENGTH)
    M = sparse(mass_matrix(s))
    K = sparse(stiffness_matrix(s))
    Λ = if state === :dirichlet
        cholesky(Symmetric(K)) \ Matrix(M)
    else
        # The embedding V_D ⊂ V, per axis, then tensored. `basis_values` builds its tables as
        # kron over the axes in REVERSE order because the flat index runs the first axis
        # fastest, so the embedding has to be built the same way or the two disagree silently.
        R = recombination_matrix(BSplineBasis(UniformMesh(n, SQUARE_LENGTH), p, Dirichlet()))
        E = sparse(kron(R, R))
        E * (cholesky(Symmetric(Matrix(E' * K * E))) \ Matrix(E' * M))
    end
    # M K⁻¹ M is symmetric exactly; the factorisation makes it so only to round-off, and
    # `QuadraticHamiltonian` checks symmetry rather than assuming it.
    A = M * Λ
    EulerSquare(s, Matrix(Λ), M, Matrix((A .+ A') ./ 2),
        Vector{Float64}(basis_integrals(s)), state)
end

space(sq::EulerSquare) = sq.space
nbasis(sq::EulerSquare) = nbasis(sq.space)

function Base.show(io::IO, sq::EulerSquare)
    print(io, "EulerSquare(n=", ncells(sq.space)[1], ", p=", degree(sq.space)[1],
        ", N=", nbasis(sq.space), ", state=:", sq.state, ")")
end

"``\\int_\\Omega \\omega_h \\, dx``."
integrate(sq::EulerSquare, ω̂::AbstractVector) = dot(sq.b, ω̂)

"``(\\omega_h, v_h)_{L^2}``, through the mass matrix."
l2inner(sq::EulerSquare, ω̂::AbstractVector, v̂::AbstractVector) = dot(ω̂, sq.M, v̂)

l2norm(sq::EulerSquare, ω̂::AbstractVector) = sqrt(max(l2inner(sq, ω̂, ω̂), 0.0))
mean_value(sq::EulerSquare, ω̂::AbstractVector) = integrate(sq, ω̂) / SQUARE_LENGTH^2

@doc raw"""
    dirichlet_eigenvalue(sq)

The discrete approximation of ``\lambda_{1,1} = 2\pi^2`` on this space, as the reciprocal of the
largest eigenvalue of ``\Lambda``.

Read off ``\Lambda`` rather than from ``\mathbb{K}\hat{v} = \lambda\mathbb{M}\hat{v}`` so that
one method serves both state spaces. In the `:free` case the free stiffness matrix has the
constants in its kernel and its smallest generalised eigenvalue is zero, which is not the
Dirichlet eigenvalue and not an approximation of it; ``\Lambda`` is the Dirichlet solution
operator in either space, so ``1/\lambda_{\max}(\Lambda)`` is the right number by construction.

Reported alongside [`DIRICHLET_EIGENVALUE`](@ref) wherever a relaxed state is measured against
it, because that is the number the run can actually reach: a fitted ``\lambda`` that agreed with
``2\pi^2`` more closely than the space's own eigenvalue does would be a coincidence, not a
better result.
"""
dirichlet_eigenvalue(sq::EulerSquare) = 1 / maximum(real, eigvals(sq.Λ))

@doc raw"""
    euler_state(sq, spec)

The initial state of run `spec`, as the ``L^2`` projection of [`initial_condition`](@ref) onto
the state space of `sq`.

No mean is subtracted, unlike Section 4's [`spline_state`](@ref): there the state was
``\omega = u - u_\Omega`` because the periodic Poisson solve needs a mean-free right-hand side
and ``u_\Omega`` is a constant of the motion. Here the Dirichlet solve needs neither, the mean
is not conserved, and the constant function is not even in the space to subtract.
"""
function euler_state(sq::EulerSquare, spec::EulerSpec)
    ω₀ = initial_condition(spec)
    project(sq.space, x -> ω₀(x[1], x[2]))
end

@doc raw"""
    euler_flow(sq, spec)

The [`MetriplecticFlow`](@ref) of run `spec`: the collision-like bracket generated by the
elliptic energy, that energy, and the run's entropy.

| `spec.entropy` | ``s(y)`` | ``S`` | mobility ``M`` | ``\partial M/\partial u`` |
|:--|:--|:--|:--|:--|
| `:quadratic` | ``y^2/2`` | `QuadraticHamiltonian(M)` | ``1`` | ``0`` |
| `:gibbs` | ``y \log y`` | [`GibbsEntropy`](@ref) | ``y`` | ``1`` |

The mobility is not an independent knob: `eq:M-condition` fixes it as
``M = 1 / \partial_y^2 s``, so the two columns are one choice. Passing the derivative is
mandatory for a mobility given as a function — [`CollisionBracket`](@ref) refuses to difference
it, because the Jacobian of a metriplectic flow is analytic exactly when that derivative is
known.

The Poisson half is absent (`bracket = nothing`), which is the case every §5 experiment runs;
the energy is ``H = \tfrac12 \int_\Omega |\nabla\phi|^2 = \tfrac12 \hat\omega^T \mathbb{M}
\Lambda \hat\omega``, and the bracket is generated by the same ``\Lambda``, which is what makes
the degeneracy ``\mathbb{G} \, \partial H/\partial\hat\omega = 0`` structural rather than
accidental.
"""
function euler_flow(sq::EulerSquare, spec::EulerSpec)
    H = QuadraticHamiltonian(sq.MΛ)
    if spec.entropy === :quadratic
        G = CollisionBracket(sq.space, sq.Λ)
        return MetriplecticFlow(sq.space, G, H, QuadraticHamiltonian(Matrix(sq.M)))
    elseif spec.entropy === :gibbs
        G = CollisionBracket(sq.space, sq.Λ; mobility = (x, u) -> u,
            mobility_derivative = (x, u) -> one(u))
        return MetriplecticFlow(sq.space, G, H, GibbsEntropy())
    end
    throw(ArgumentError("unknown entropy $(spec.entropy) for run $(spec.name)"))
end

## The closed-form references

@doc raw"""
    euler_entropy_floor(H₀; λ = DIRICHLET_EIGENVALUE)

``S_\eta = \lambda_{1,1} H_0``, the constrained minimum of ``S = \tfrac12 \int_\Omega \omega^2``
at fixed energy ``H_0`` on the unit square with homogeneous Dirichlet conditions.

At the minimiser ``\omega = \lambda_{1,1}\phi``, so
``S = \tfrac12 \lambda_{1,1}^2 \|\phi\|^2`` and ``H_0 = \tfrac12 \lambda_{1,1} \|\phi\|^2``,
whence ``S = \lambda_{1,1} H_0``. Equivalently this is the Poincaré inequality
``\|\nabla\phi\|^2 \le \lambda_{1,1}^{-1} \|\Delta\phi\|^2``, so ``S \ge S_\eta`` holds for
*any* admissible state and not only along the flow — which is why it is asserted at ``t = 0``
as a precondition rather than only at the end.

Pass `λ` to measure against the space's own [`dirichlet_eigenvalue`](@ref) instead of the
continuum one.
"""
euler_entropy_floor(H₀; λ = DIRICHLET_EIGENVALUE) = λ * H₀

@doc raw"""
    eigenmode_fit(sq, ω̂)

How close `ω̂` is to a state of the form ``\omega = \lambda \phi``, as
`(λ, residual, relative)`.

The best ``\lambda`` in ``L^2`` is ``(\omega,\phi)/(\phi,\phi)``, and since ``\phi = \Lambda
\omega`` the numerator is ``2H``, so ``\lambda = 2H / \|\phi\|^2`` — no fitting is done and none
is needed. `residual` is ``\|\omega - \lambda\phi\|_{L^2}`` and `relative` that over
``\|\omega\|_{L^2}``.

For a relaxed B1 or B2 state ``\lambda`` must approach [`dirichlet_eigenvalue`](@ref) and
`relative` must approach zero **together**, and ``\lambda`` is the weaker of the two: its error
is **second** order in `relative`. Measured on B1 and B2 alike,

```math
\frac{|\lambda - \lambda_h|}{\lambda_h} = 0.2474 \, \mathrm{relative}^2 ,
```

to five digits and three orders of magnitude apart in `relative` — which is why ``\lambda`` is
already within a per cent of the eigenvalue for states that are nothing like the eigenmode, and
why the drivers put their tolerance on ``\lambda`` at `relative^2` rather than at a constant.
"""
function eigenmode_fit(sq::EulerSquare, ω̂::AbstractVector)
    φ̂ = sq.Λ * ω̂
    λ = l2inner(sq, ω̂, φ̂) / l2inner(sq, φ̂, φ̂)
    res = l2norm(sq, ω̂ .- λ .* φ̂)
    return (λ, res, res / l2norm(sq, ω̂))
end

@doc raw"""
    gibbs_lambda(mass, S, H₀)

``\lambda = (M(\omega) + S(\omega)) / 2 H_0``, the manuscript's multiplier for B3's relaxed
state ``\omega = e^{\lambda\phi - 1}``.

It is an identity rather than a fit, and the derivation is one line: substituting
``\log\omega = \lambda\phi - 1`` into ``S = \int_\Omega \omega\log\omega`` gives
``S = \lambda \int_\Omega \omega\phi - \int_\Omega \omega = 2\lambda H - M``. It therefore
holds **only** with the ``\mu = 0`` of [`SECTION5_RUNS`](@ref), i.e. only in a space where the
mass is not conserved, and `verify_euler.jl` checks that the identity closes on the run's own
final state rather than assuming it.
"""
gibbs_lambda(mass, S, H₀) = (mass + S) / (2H₀)

@doc raw"""
    interior_weights(sq; margin = 0)

The quadrature weights of `sq`, zeroed at every node within `margin` of ``\partial\Omega``.

Only B3 uses it, and only to separate the boundary layer from the interior. The manuscript's
reference ``e^{\lambda\phi - 1}`` does not vanish on ``\partial\Omega``:
``\phi|_{\partial\Omega} = 0`` leaves it at ``e^{-1} \approx 0.368``, while every
``\omega_h \in V_D`` is zero there. B3 runs in ``V_D`` all the same — that is what forces
``\mu = 0`` and makes the ``\lambda`` formula an identity ([`SECTION5_RUNS`](@ref)) — so the
outermost cell is where no state of the space can match the reference, and an unweighted fit of
``\log\omega = \lambda\phi + \mu - 1`` would be a fit of that cell. `margin = 2h` is what
[`gibbs_fit`](@ref) is called with for that reason.

`run_b3.jl` reports the residual at four margins on top of that, and asserts it **falls** with
the margin: a residual that collapses when the boundary is excluded localises the disagreement
to the boundary layer, which is what B3 claims, and one that does not would be a disagreement
in the interior.

B1's and B2's reference ``\omega = \lambda_{1,1}\phi`` vanishes on ``\partial\Omega`` together
with ``\phi``, so nothing needs excluding there at all.
"""
function interior_weights(sq::EulerSquare; margin = 0.0)
    w = copy(quadrature_weights(sq.space))
    margin <= 0 && return w
    L = SQUARE_LENGTH
    for (r, x) in enumerate(quadrature_nodes(sq.space))
        (minimum(x) < margin || maximum(x) > L - margin) && (w[r] = 0.0)
    end
    return w
end

@doc raw"""
    gibbs_fit(sq, ω̂; margin = 0)

The multipliers ``(\lambda, \mu)`` of ``\log\omega = \lambda\phi + \mu - 1`` fitted to `ω̂` by
weighted least squares over the quadrature grid, as `(λ, μ, residual)`.

This is the **control** on [`gibbs_lambda`](@ref), not a substitute for it. The manuscript's
reference sets ``\mu = 0``; a fit that returns ``\mu`` far from zero would say that the discrete
equilibrium carries a mass multiplier after all, which is a statement about the space rather
than about the run, and it would invalidate the ``\lambda`` formula rather than merely shift it.
Both numbers are therefore reported side by side.

`residual` is the weighted root-mean-square of ``\log\omega - \lambda\phi - \mu + 1`` relative
to the spread of ``\log\omega``, which is the quantity that says whether *any* member of the
family fits. `margin` excludes the boundary layer — see [`interior_weights`](@ref).
"""
function gibbs_fit(sq::EulerSquare, ω̂::AbstractVector; margin = 0.0)
    s = sq.space
    w = interior_weights(sq; margin = margin)
    u = field(s, ω̂, (0, 0))
    φ = field(s, sq.Λ * ω̂, (0, 0))
    y = log.(u) .+ 1
    # Two normal equations for (λ, μ) against the weighted basis (φ, 1).
    W = sum(w)
    a11, a12, a22 = dot(w, φ .* φ), dot(w, φ), W
    r1, r2 = dot(w, φ .* y), dot(w, y)
    det = a11 * a22 - a12^2
    λ = (r1 * a22 - r2 * a12) / det
    μ = (r2 * a11 - r1 * a12) / det
    r = y .- λ .* φ .- μ
    spread = sqrt(max(dot(w, (y .- r2 / W) .^ 2) / W, 0.0))
    return (λ, μ, sqrt(max(dot(w, r .^ 2) / W, 0.0)) / max(spread, 1e-300))
end

@doc raw"""
    gibbs_residual(sq, ω̂, λ, μ = 0; margin = 0)

The relative ``L^2`` distance between `ω̂` and ``e^{\lambda\phi + \mu - 1}``, measured on the
quadrature grid,

```math
\frac{\big( \int |\omega_h - e^{\lambda\phi_h + \mu - 1}|^2 \big)^{1/2}}
     {\big( \int |\omega_h|^2 \big)^{1/2}} ,
```

both integrals over the same node set, which `margin` restricts to the interior — see
[`interior_weights`](@ref).

``\mu = 0`` is the manuscript's own reference; passing the fitted ``\mu`` of
[`gibbs_fit`](@ref) instead measures the same state against the equilibrium family the discrete
flow actually has. The difference between the two numbers is what the mass multiplier costs, and
`run_b3.jl` reports both.

The reference is evaluated pointwise rather than projected, so this compares the run against
the closed form and not against the closed form's own discretisation.
"""
function gibbs_residual(sq::EulerSquare, ω̂::AbstractVector, λ, μ = 0.0; margin = 0.0)
    s = sq.space
    w = interior_weights(sq; margin = margin)
    u = field(s, ω̂, (0, 0))
    r = u .- exp.(λ .* field(s, sq.Λ * ω̂, (0, 0)) .+ μ .- 1)
    return sqrt(dot(w, r .^ 2)) / sqrt(dot(w, u .^ 2))
end

@doc raw"""
    state_extrema(sq, ω̂)

`(minimum, maximum)` of ``\omega_h`` over the quadrature grid.

The minimum is the quantity B3 depends on: ``y \log y`` and the mobility ``M = y`` are defined
only for ``\omega > 0``, and the grid is where they are evaluated. See [`GibbsEntropy`](@ref).
"""
state_extrema(sq::EulerSquare, ω̂::AbstractVector) = extrema(field(sq.space, ω̂, (0, 0)))
