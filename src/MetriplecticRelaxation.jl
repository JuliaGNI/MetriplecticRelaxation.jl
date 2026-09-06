module MetriplecticRelaxation

using CairoMakie
using FFTW
using LinearAlgebra
using Printf
using SparseArrays

# `mass_matrix`, `stiffness_matrix`, `evaluate`, `project` and friends are exported by both
# PoissonBrackets and SimpleSplines, which leaves every one of them ambiguous under a bare
# `using` of the two. PoissonBrackets is the entry point, and what is needed is named.
using PoissonBrackets
using PoissonBrackets: DiscreteSpace, DiscreteHamiltonian, TensorSplineSpace,
                       DoubleBracket, ProjectorBracket, MetriplecticFlow,
                       QuadraticHamiltonian,
                       mass_matrix, mass_factorization, stiffness_matrix,
                       tensor_weighted_matrix, basis_integrals, ncells, degree,
                       project, evaluate, vectorfield

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

export Diagnostics, potential, energy, entropy, vorticity_mass, potential_norm²,
       Trace, record!, energy_error, entropy_monotone, best_fit_euler,
       cone_coordinates, cone_residual, fit_rate, scatter_data

include("diagnostics.jl")

export figure_fields, figure_traces, figure_scatter, figure_cone, figure_rates, figure_tau

include("figures.jl")

end
