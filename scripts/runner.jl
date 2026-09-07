#
# The shared driver behind `run_a1.jl` ... `run_a4.jl`.
#
#     include(joinpath(@__DIR__, "runner.jl"))
#     using .Runner
#
# Everything the four drivers have in common: option parsing, the time loop, and writing the
# run out. What differs between them is only the checks, which is what each `run_aN.jl` holds.
#
# Output paths come from `--runs-dir` and `--results-dir` and default to the repository's own
# `runs/` and `results/`, neither of which is tracked. They are resolved at call time and never
# baked into a constant, so the same script can write into a talk's figure directory without
# being copied.

module Runner

using MetriplecticRelaxation
using MetriplecticRelaxation: SpectralTorus, SplineTorus, Diagnostics, Trace,
                              spectral_state, spectral_rhs, spectral_step,
                              spline_state, spline_rhs, spline_step, record!,
                              EulerSquare, euler_state, euler_flow, state_extrema,
                              GradShafranovBox, gs_state, gs_flow, gs_fit
using PoissonBrackets: Integrator, ImplicitMidpoint, integrate_step!
using Printf
using Serialization

export Options, parse_options, run_spectral, run_spline, run_euler, run_gs, save_run,
       section, report

@doc raw"""
    Options

Where a driver writes, and at what resolution it runs.

`spectral` defaults to **256**, which is the manuscript's own grid for every Section 4 run.
`cells` is the number of spline cells per direction and has no counterpart in the
manuscript — the spline discretisation is this reproduction's deviation — so it is a recorded
choice: **64** cubic cells, i.e. 4096 degrees of freedom against the spectral run's 65 536.
That is chosen so the two runs cost comparably rather than so they have equal resolution; the
spline space is 4th-order accurate and the spectral one is not, so equal degrees of freedom
would not mean equal accuracy either.
"""
struct Options
    runs_dir::String
    results_dir::String
    spectral::Int
    cells::Int
    degree::Int
    samples::Int
    quiet::Bool
end

"""
    parse_options(args; kwargs...)

Read `--runs-dir`, `--results-dir`, `--spectral`, `--cells`, `--degree`, `--samples` and
`--quiet` from `args`, falling back to the repository's own directories.

An option that takes a value is rejected if the value is missing or is itself an option. Without
that, `--runs-dir --results-dir out` silently binds the flag name as a path and `mkpath` then
creates a directory called `--results-dir`, which is how one came to be sitting in the repository
root.
"""
function parse_options(args = ARGS; spectral = 256, cells = 64, degree = 3, samples = 400)
    root = dirname(@__DIR__)
    runs = joinpath(root, "runs")
    results = joinpath(root, "results")
    quiet = false

    function value(i)
        i < length(args) || throw(ArgumentError("$(args[i]) needs a value"))
        v = args[i + 1]
        startswith(v, "--") &&
            throw(ArgumentError("$(args[i]) needs a value, got the option $(v)"))
        return v
    end

    i = 1
    while i <= length(args)
        a = args[i]
        if a == "--runs-dir"
            runs = value(i)
            i += 1
        elseif a == "--results-dir"
            results = value(i)
            i += 1
        elseif a == "--spectral"
            spectral = parse(Int, value(i))
            i += 1
        elseif a == "--cells"
            cells = parse(Int, value(i))
            i += 1
        elseif a == "--degree"
            degree = parse(Int, value(i))
            i += 1
        elseif a == "--samples"
            samples = parse(Int, value(i))
            i += 1
        elseif a == "--quiet"
            quiet = true
        else
            throw(ArgumentError("unknown option $(a)"))
        end
        i += 1
    end
    mkpath(runs)
    mkpath(results)
    return Options(runs, results, spectral, cells, degree, samples, quiet)
end

"A section rule, matching the `Checks` harness's own headers."
function section(text)
    println()
    println(text)
    println("-"^length(text))
    flush(stdout)
end

@doc raw"""
    run_spectral(spec, opts)

Integrate `spec` on a [`SpectralTorus`](@ref) of `opts.spectral` points with classical RK4 at
the manuscript's own time step, returning `(grid, diagnostics, trace, u_Ω)`.

`opts.samples` samples are recorded, evenly spaced in step count. `stdout` is flushed at every
progress line: Julia buffers it when redirected, so an interrupted run would otherwise lose
all of its progress output.
"""
function run_spectral(spec, opts::Options; observer = nothing)
    g = SpectralTorus(opts.spectral)
    d = Diagnostics(g, spec)
    ω, uΩ = spectral_state(g, spec)
    # Built once, outside the loop, exactly as `spline_rhs` is below: for A1 the closure holds
    # a fixed `X_h`, and rebuilding it per step would recompute that 10⁵ times.
    rhs = spectral_rhs(g, spec)

    nsteps = round(Int, spec.T / spec.Δt)
    stride = max(1, nsteps ÷ opts.samples)
    tr = Trace(ω)
    record!(tr, d, 0.0, ω)
    observer === nothing || observer(g, 0.0, ω)

    t0 = time()
    for k in 1:nsteps
        ω = spectral_step(rhs, ω, spec.Δt)
        if k % stride == 0 || k == nsteps
            record!(tr, d, k * spec.Δt, ω)
            observer === nothing || observer(g, k * spec.Δt, ω)
            if !opts.quiet && (k % (10stride) == 0 || k == nsteps)
                @printf("    spectral %6.1f%%   t = %7.3f   H = %+.10e   S = %.10e   [%.0f s]\n",
                    100k / nsteps, k * spec.Δt, tr.H[end], tr.S[end], time() - t0)
                flush(stdout)
            end
        end
    end
    return (g, d, tr, uΩ)
end

@doc raw"""
    run_spline(spec, opts)

The same run on a [`SplineTorus`](@ref) of `opts.cells` cells and degree `opts.degree`,
returning `(torus, diagnostics, trace, u_Ω)`.

Same integrator and same time step as [`run_spectral`](@ref): the two runs differ only in the
spatial discretisation, which is what makes their agreement a statement about the equation.
"""
function run_spline(spec, opts::Options; observer = nothing)
    t = SplineTorus(opts.cells, opts.degree)
    d = Diagnostics(t, spec)
    ω̂, uΩ = spline_state(t, spec)
    rhs = spline_rhs(t, spec)

    nsteps = round(Int, spec.T / spec.Δt)
    stride = max(1, nsteps ÷ opts.samples)
    tr = Trace(ω̂)
    record!(tr, d, 0.0, ω̂)
    observer === nothing || observer(t, 0.0, ω̂)

    t0 = time()
    for k in 1:nsteps
        ω̂ = spline_step(rhs, ω̂, spec.Δt)
        if k % stride == 0 || k == nsteps
            record!(tr, d, k * spec.Δt, ω̂)
            observer === nothing || observer(t, k * spec.Δt, ω̂)
            if !opts.quiet && (k % (10stride) == 0 || k == nsteps)
                @printf("    spline   %6.1f%%   t = %7.3f   H = %+.10e   S = %.10e   [%.0f s]\n",
                    100k / nsteps, k * spec.Δt, tr.H[end], tr.S[end], time() - t0)
                flush(stdout)
            end
        end
    end
    return (t, d, tr, uΩ)
end

@doc raw"""
    run_euler(spec, opts; Δt = spec.Δt, T = spec.T, cells = opts.cells,
              degree = opts.degree, observer = nothing)

Integrate a Section 5.4 run on an [`EulerSquare`](@ref) with `ImplicitMidpoint` — which for
this field **is** Crank-Nicolson, the manuscript's own method — returning
`(square, diagnostics, flow, trace)`.

The nonlinear solve is [`Integrator`](@ref)'s existing `SimpleSolvers.NewtonSolver` with the
flow's analytic Jacobian and a **dense** LU. No solver code is written here and none should be:
``\phi = \delta H/\delta u`` comes from an elliptic solve, so the Jacobian carries a
structurally dense block, and a sparse solver handed that Jacobian would be handed one with the
``\phi`` coupling dropped — Newton would degrade to a slow quasi-Newton and report nothing.

`û₀` is passed to `Integrator`, which is what puts the residual tolerance at the round-off
floor of *this* field's amplitude rather than at that of a field of order one. Newton runs to
that tolerance and not to a fixed iteration count, because the claim being reproduced is energy
conservation to machine precision and a lagged Jacobian converges only linearly — see
`default_f_abstol`.

`Δt`, `T`, `cells` and `degree` are overridable because they are choices rather than data, and
the step-size study `run_b3.jl` performs varies the first of them.
"""
function run_euler(spec, opts::Options; Δt = spec.Δt, T = spec.T, cells::Int = opts.cells,
        degree::Int = opts.degree, observer = nothing)
    sq = EulerSquare(cells, degree; state = spec.state)
    d = Diagnostics(sq, spec)
    ω̂ = euler_state(sq, spec)
    f = euler_flow(sq, spec)
    integ = Integrator(f, ImplicitMidpoint(), Δt; û₀ = ω̂)

    nsteps = round(Int, T / Δt)
    stride = max(1, nsteps ÷ opts.samples)
    tr = Trace(ω̂)
    record!(tr, d, 0.0, ω̂)
    observer === nothing || observer(sq, f, 0.0, ω̂)

    t0 = time()
    for k in 1:nsteps
        integrate_step!(ω̂, integ)
        if k % stride == 0 || k == nsteps
            record!(tr, d, k * Δt, ω̂)
            observer === nothing || observer(sq, f, k * Δt, ω̂)
            if !opts.quiet && (k % (10stride) == 0 || k == nsteps)
                (lo, hi) = state_extrema(sq, ω̂)
                @printf("    euler %6.1f%%   t = %8.3f   H = %+.10e   S = %.10e   ω ∈ [%+.2e, %.3f]   [%.0f s]\n",
                    100k / nsteps, k * Δt, tr.H[end], tr.S[end], lo, hi, time() - t0)
                flush(stdout)
            end
        end
    end
    return (sq, d, f, tr)
end

@doc raw"""
    run_gs(spec, opts; Δt = spec.Δt, T = spec.T, cells = spec.cells,
           degree = spec.degree, state = :dirichlet, observer = nothing)

Integrate a Section 5.5 Grad-Shafranov run on a [`GradShafranovBox`](@ref) with
`ImplicitMidpoint`, returning `(box, diagnostics, flow, trace)`.

Same integrator, same nonlinear solve and same reasoning about the dense Jacobian as
[`run_euler`](@ref); what differs is the geometry, the measure and the state variable
``j = u/r``, all of which are the box's business rather than the loop's.

The resolution comes from `spec` rather than from `opts` — §5.5's mesh is **anisotropic**, so a
single `--cells` would have to pick one of the two axes, and the ratio is a recorded choice
([`SECTION55_RUNS`](@ref)) rather than a knob. `cells` overrides it as a pair for the
convergence and control studies the driver runs.
"""
function run_gs(spec, opts::Options; Δt = spec.Δt, T = spec.T, cells = spec.cells,
        degree::Int = spec.degree, state::Symbol = :dirichlet, observer = nothing)
    box = GradShafranovBox(cells, degree; state = state)
    d = Diagnostics(box, spec)
    ĵ = gs_state(box, spec)
    f = gs_flow(box)
    integ = Integrator(f, ImplicitMidpoint(), Δt; û₀ = ĵ)

    nsteps = round(Int, T / Δt)
    stride = max(1, nsteps ÷ opts.samples)
    tr = Trace(ĵ)
    record!(tr, d, 0.0, ĵ)
    observer === nothing || observer(box, f, 0.0, ĵ)

    t0 = time()
    for k in 1:nsteps
        integrate_step!(ĵ, integ)
        if k % stride == 0 || k == nsteps
            record!(tr, d, k * Δt, ĵ)
            observer === nothing || observer(box, f, k * Δt, ĵ)
            if !opts.quiet && (k % (10stride) == 0 || k == nsteps)
                (λ, _, rel) = gs_fit(box, ĵ)
                @printf("    gs %6.1f%%   t = %7.3f   H = %+.10e   S = %.10e   λ = %.8f   rel = %.3e   [%.0f s]\n",
                    100k / nsteps, k * Δt, tr.H[end], tr.S[end], λ, rel, time() - t0)
                flush(stdout)
            end
        end
    end
    return (box, d, f, tr)
end

@doc raw"""
    save_run(opts, name, payload)

Serialise `payload` to `<runs-dir>/<name>.jls`, and write the traces it contains to
`<results-dir>/<name>_trace.csv`.

`Serialization` rather than JLD2 or HDF5: `runs/` is regenerable and untracked, nothing
outside this repository reads it, and a stdlib does not need a `[compat]` bound. The CSV is
what a figure or a table is built from, and is the form a number quoted in prose can be traced
back to.
"""
function save_run(opts::Options, name::AbstractString, payload)
    path = joinpath(opts.runs_dir, name * ".jls")
    open(path, "w") do io
        serialize(io, payload)
    end

    csv = joinpath(opts.results_dir, name * "_trace.csv")
    open(csv, "w") do io
        println(io, "method,t,H,S,phi2,mass")
        for (label, tr) in payload.traces, i in eachindex(tr.t)

            @printf(io, "%s,%.10e,%.16e,%.16e,%.16e,%.10e\n",
                label, tr.t[i], tr.H[i], tr.S[i], tr.φ²[i], tr.M[i])
        end
    end
    return (path, csv)
end

"""
    report(opts, name, lines)

Write the run's own numbers to `<results-dir>/<name>.md`, as the record of what the scripts
printed. Prose that quotes a number quotes it from here.
"""
function report(opts::Options, name::AbstractString, lines)
    path = joinpath(opts.results_dir, name * ".md")
    open(path, "w") do io
        for l in lines
            println(io, l)
        end
    end
    return path
end

end # module Runner
