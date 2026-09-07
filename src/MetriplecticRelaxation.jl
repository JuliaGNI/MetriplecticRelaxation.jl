module MetriplecticRelaxation

using CairoMakie
using FFTW
using LinearAlgebra
using Printf
using SparseArrays

# `mass_matrix`, `stiffness_matrix`, `evaluate`, `project` and friends are exported by both
# PoissonBrackets and SimpleSplines, which leaves every one of them ambiguous under a bare
# `using` of the two. PoissonBrackets is the entry point, and what is needed is named.
#
# `field` is the one name here that PoissonBrackets neither exports nor declares public
# (`Base.ispublic(PoissonBrackets, :field)` is `false`), so this reaches into its internals and
# a rename upstream breaks this package silently. It is non-public in SimpleSplines too, so
# there is nowhere public to take it from and the fix belongs upstream, not here — an `export`
# or a `public` declaration on whichever of the two owns it.
using PoissonBrackets
using PoissonBrackets: DiscreteSpace, DiscreteHamiltonian, TensorSplineSpace,
                       DoubleBracket, ProjectorBracket, CollisionBracket, MetriplecticFlow,
                       QuadraticHamiltonian,
                       mass_matrix, mass_factorization, stiffness_matrix,
                       tensor_weighted_matrix, basis_integrals, ncells, degree,
                       quadrature_nodes, quadrature_weights, basis_values, field,
                       project, evaluate, vectorfield

# The boundary conditions and the V_D ⊂ V embedding Section 5 needs. PoissonBrackets re-exports
# SimpleSplines' assembly interface but neither of these: `TensorSplineSpace(n, p, bc)`
# dispatches on the condition type, and `recombination_matrix` is what expresses a
# homogeneous-Dirichlet basis function in the clamped one — see `EulerSquare`.
using SimpleSplines: Dirichlet, Free, BSplineBasis, UniformMesh, recombination_matrix

# Extended, not shadowed. The spectral grid and the spline space are two more discretisations
# of the objects these generic functions already name, so they get methods rather than
# same-named functions of their own — otherwise `canonical_bracket` would mean one thing on a
# `TorusGrid` and an unrelated thing on a `SpectralTorus`, with no dispatch between them.
import PoissonBrackets: canonical_bracket, hamiltonian_field, integrate, space, nbasis,
                        hamiltonian, gradient, hessian, entropy

export DOMAIN_LENGTH, DOMAIN_AREA
export Gaussian, islands_h, ISLAND_CENTRES, CENTRAL_ISLANDS
export agm, contour_length, contour_length_quadrature, contour_average, relaxation_time,
       contour_samples, contour_deviation
export analytic_minimiser, analytic_entropy_minimum
export euler_minimiser, euler_entropy_minimum
export RunSpec, SECTION4_RUNS, SECTION4_ORDER, initial_condition, periodise

include("torus.jl")

export SpectralTorus, torus_field, poisson_periodic, canonical_bracket, hamiltonian_field,
       double_bracket_field, parallel_diffusion, projector_bracket_field,
       spectral_state, spectral_rhs, spectral_step

include("spectral.jl")

export SplineTorus, PoissonMap, LinearHamiltonian, EllipticEnergy,
       spline_state, spline_flow, spline_rhs, spline_step, spline_grid,
       fixed_double_operator

include("spline.jl")

export SQUARE_LENGTH, DIRICHLET_EIGENVALUE
export EulerSpec, SECTION5_RUNS, SECTION5_ORDER, gaussian_w2, perturbation_b2, B3_FLOOR
export EulerSquare, GibbsEntropy, euler_state, euler_flow, dirichlet_eigenvalue,
       euler_entropy_floor, eigenmode_fit, gibbs_lambda, gibbs_fit, gibbs_residual,
       interior_weights, state_extrema, euler_axis, euler_grid

include("euler.jl")

export HERRNEGGER_C, HERRNEGGER_D, GS_RADIAL, GS_AXIAL,
       GS_LAMBDA_RECTANGLE, GS_LAMBDA_CONTINUUM
export gs_density, herrnegger_mobility, herrnegger_profile
export gs_stiffness, gs_profile_matrix, gs_eigenvalue, separable_eigenvalue
export TakedaGrid, takeda_iterate, takeda_eigenvalue, takeda_field,
       takeda_order, takeda_extrapolate
export DISK_MAP, GS_LAMBDA_DISK, GS_LAMBDA_DISK_CONTINUUM
export disk_map, DiskTriangulation, disk_area, disk_matrices, disk_eigenvalue

include("takeda.jl")

export GSSpec, SECTION55_RUNS, SECTION55_ORDER, gs_entropy_weight
export GradShafranovBox, gs_state, gs_flow, gs_current, gs_ordinate, gs_fit, gs_rayleigh

include("gradshafranov.jl")

export Diagnostics, potential, energy, entropy, vorticity_mass, potential_norm²,
       Trace, record!, energy_error, entropy_monotone, entropy_plateau, best_fit_euler,
       cone_coordinates, cone_residual, fit_rate, scatter_data

include("diagnostics.jl")

export figure_fields, figure_traces, figure_scatter, figure_cone, figure_rates, figure_tau

include("figures.jl")

end
