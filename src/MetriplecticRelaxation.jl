module MetriplecticRelaxation

using FFTW
using LinearAlgebra
using Printf

export DOMAIN_LENGTH, DOMAIN_AREA
export Gaussian, islands_h, ISLAND_CENTRES, CENTRAL_ISLANDS
export agm, contour_length, contour_length_quadrature, contour_average, relaxation_time
export analytic_minimiser, analytic_entropy_minimum
export euler_minimiser, euler_entropy_minimum
export RunSpec, SECTION4_RUNS

include("torus.jl")

export SpectralTorus, torus_field, poisson_periodic, canonical_bracket, hamiltonian_field,
       double_bracket_field, projector_bracket_field, spectral_state, spectral_step!

include("spectral.jl")

end
