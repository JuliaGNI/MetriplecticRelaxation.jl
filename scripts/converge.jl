#!/usr/bin/env julia
#
# The refinement study behind the deliberate deviation.
#
#     julia --project=scripts scripts/converge.jl [--results-dir DIR]
#
# The reproduction replaces the manuscript's Fourier spectral method with a B-spline Galerkin
# one, and the whole justification for reading the two as one result is that they agree. But
# "they agree at this resolution" is a weak statement: it is consistent with two methods being
# wrong in the same way, and it says nothing about which one to believe when they differ.
#
# What settles it is that the difference CONVERGES, at the order the spline space has. A degree-p
# B-spline approximation is O(h^{p+1}) in L2, so refining the mesh must divide the difference
# from the spectral run by 2^{p+1} -- 16 for the cubic spaces used throughout. A measured order
# near p+1 says the spectral run is the exact answer as far as the spline run can tell, which is
# exactly what a reference is.
#
# This runs to a SHORT final time on purpose. The question is whether the two discretisations
# solve the same equation, which is a statement about the operator and is visible in a few
# hundred steps; running to T would fold in the whole relaxation and measure something else.
#
# The control that can fail is the last section: the observed order must degrade for A4, whose
# initial condition as printed is discontinuous, since no method converges at 4th order to a
# field that is not in the space at any resolution.

using MetriplecticRelaxation
using MetriplecticRelaxation: SECTION4_RUNS, SpectralTorus, SplineTorus, RunSpec,
                              spectral_state, spectral_step!, spline_state, spline_rhs,
                              spline_step!, spline_grid, torus_field, periodise
using Printf

include(joinpath(@__DIR__, "check.jl"))
include(joinpath(@__DIR__, "runner.jl"))
using .Checks: header, check, summary
using .Runner

const opts = parse_options()

# The spectral reference. 192 points is well past the resolution any spline mesh below reaches,
# so its own error does not enter the comparison.
const NREF = 192
const CELLS = (16, 24, 32, 48, 64)
const NSTEPS = 200

"""
    reference(spec)

The spectral solution of `spec` after `NSTEPS` steps on the `NREF`-point grid.

Computed once per `spec` and reused across every mesh, since it does not depend on the spline
resolution. Recomputing it inside [`difference`](@ref) would repeat a 192² run of 200 steps for
each of the five meshes, four runs and three degrees — some forty times over, for an answer that
is the same every time.
"""
function reference(spec)
    g = SpectralTorus(NREF)
    ωs, _ = spectral_state(g, spec)
    ĥ = spec.h === nothing ? nothing : torus_field(g, spec.h)
    for _ in 1:NSTEPS
        ωs = spectral_step!(g, ωs, spec, ĥ)
    end
    return ωs
end

"""
    difference(spec, ref, cells, degree)

The relative `L²` difference between the spline solution of `spec` on `cells` cells of degree
`degree` and the spectral reference `ref`, sampled on the spectral grid.
"""
function difference(spec, ref, cells::Int, degree::Int)
    t = SplineTorus(cells, degree)
    ω̂, _ = spline_state(t, spec)
    rhs = spline_rhs(t, spec)
    for _ in 1:NSTEPS
        ω̂ = spline_step!(rhs, ω̂, spec.Δt)
    end
    d = spline_grid(t, ω̂, NREF) .- ref
    return sqrt(sum(abs2, d) / sum(abs2, ref))
end

"The least-squares order of convergence of `errs` against `cells`."
function observed_order(cells, errs)
    x = log.(collect(Float64, cells))
    y = log.(collect(Float64, errs))
    x̄, ȳ = sum(x) / length(x), sum(y) / length(y)
    return -sum((x .- x̄) .* (y .- ȳ)) / sum((x .- x̄) .^ 2)
end

println("Refinement of the spline discretisation against a $(NREF)² spectral reference")
@printf("    %d steps, degree 3 unless stated\n", NSTEPS)
flush(stdout)

const RESULTS = Dict{String, Any}()

for name in ("a1", "a2", "a3", "a4")
    spec = SECTION4_RUNS[name]
    header("$(name): L² difference from the spectral reference after $(NSTEPS) steps")
    ref = reference(spec)
    errs = Float64[]
    for n in CELLS
        e = difference(spec, ref, n, 3)
        push!(errs, e)
        r = length(errs) > 1 ? errs[end - 1] / errs[end] : NaN
        check(@sprintf("%s  %3d cells", name, n), true,
            @sprintf("rel L² = %.4e%s", e,
                isnan(r) ? "" :
                @sprintf("   ratio %.2f  (order %.2f)", r,
                    log2(r) / log2(CELLS[length(errs)] / CELLS[length(errs) - 1]) *
                    log2(2))))
    end
    p = observed_order(CELLS, errs)
    RESULTS[name] = (; cells = CELLS, errs = errs, order = p)

    # Cubic B-splines are 4th-order accurate in L², so the difference from the reference must
    # fall at that rate — for every run whose initial condition is actually in the space to
    # approximate. A4's is not: as printed it is discontinuous on T², so it is the one case
    # where the order must be WORSE, and if it were not, the discontinuity finding would be
    # wrong.
    if name == "a4"
        check("$(name): the order DEGRADES (its u₀ is discontinuous)", p < 2.5,
            @sprintf("observed order %.2f   (cubic splines give 4 for a smooth u₀)", p))
    else
        check("$(name): the order is that of cubic splines", p > 3.0,
            @sprintf("observed order %.2f   (expected ≈ 4)", p))
    end
end

# =============================================================================================
header("A4 periodised: the order must come back")

# The control for the control. If periodising A4's Gaussian restores 4th-order convergence, the
# degradation above is the initial condition's discontinuity and nothing else.
let spec = periodise(SECTION4_RUNS["a4"]), ref = reference(spec)
    errs = [difference(spec, ref, n, 3) for n in CELLS]
    p = observed_order(CELLS, errs)
    RESULTS["a4_periodic"] = (; cells = CELLS, errs = errs, order = p)
    for (n, e) in zip(CELLS, errs)
        check(@sprintf("a4-periodic  %3d cells", n), true, @sprintf("rel L² = %.4e", e))
    end
    check("periodising A4 RESTORES the spline order", p > 3.0,
        @sprintf("order %.2f (periodised) vs %.2f (as printed)", p, RESULTS["a4"].order))
end

# =============================================================================================
header("the degree is what sets the order")

# One more way for the claim to be wrong: if the observed order were an artefact of the time
# stepping or of the reference rather than of the spline space, it would not move with p.
let spec = SECTION4_RUNS["a3"], cells = (16, 24, 32, 48), ref = reference(spec)
    for degree in (2, 3, 4)
        errs = [difference(spec, ref, n, degree) for n in cells]
        p = observed_order(cells, errs)
        check(@sprintf("degree %d gives order ≈ %d", degree, degree + 1),
            abs(p - (degree + 1)) < 1.0,
            @sprintf("observed %.2f   errors %.2e -> %.2e", p, errs[1], errs[end]))
    end
end

# =============================================================================================
section("writing")

let lines = ["# Refinement of the spline discretisation", "",
        @sprintf("A %d² Fourier spectral run is the reference; %d steps at each run's own Δt.",
            NREF, NSTEPS),
        "The manuscript's method is the spectral one, so this measures whether the",
        "B-spline deviation converges TO it, and at what order.", "",
        "| run | " * join(string.(CELLS) .* " cells", " | ") * " | order |",
        "|:--|" * repeat("--:|", length(CELLS) + 1)]
    for name in ("a1", "a2", "a3", "a4", "a4_periodic")
        r = RESULTS[name]
        push!(lines,
            "| $(name) | " * join([@sprintf("%.3e", e) for e in r.errs], " | ") *
            @sprintf(" | %.2f |", r.order))
    end
    push!(lines, "",
        "Cubic B-splines are 4th-order accurate in L². A4 as printed is the exception and",
        "the reason is its initial condition: `eq:initial_gaussian` is written unmodified on",
        "T², and A4's Gaussian is centred π/2 from the boundary with w₂ = 1 and amplitude",
        "1.8, so it is discontinuous there by 8.5 % of its peak. Periodising it restores the",
        "order, which is what identifies the initial condition rather than the solver.")
    println("    report  -> ", report(opts, "converge", lines))
end

summary("converge.jl")
