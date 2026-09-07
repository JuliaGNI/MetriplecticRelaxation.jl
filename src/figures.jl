#
# The figures of Section 4, regenerated from a saved run.
#
# Figures are NEVER committed: `results/` is gitignored and everything here is rebuilt from the
# `runs/` payload a driver wrote. The numbers quoted in prose come from the driver's own report,
# not from reading a plot.
#
# CairoMakie is loaded by the package rather than by an extension. An extension would be the
# right call if plotting were optional for a *library*; here the repository exists to produce
# these figures, `CairoMakie` is already a hard dependency in `Project.toml`, and a weak
# dependency would only add a loading path that is always taken.
#

@doc raw"""
    figure_fields(path, x₁, x₂, initial, final, contours; kwargs...)

The upper panels of Figs. 1, 2 and 4: the initial and final state as colour maps, with the
contours of the generating field over them.

`threshold` clips the colour map from below, which is the manuscript's own choice for Fig. 1 —
"for clarity, in the color maps we display the solution ``u`` only where ``u \geq 10^{-4}``".
Passing `nothing` disables it, which is what the reduced Euler figures want, their vorticity
being signed.
"""
function figure_fields(path, x₁, x₂, initial, final, contours;
        threshold = nothing, title_initial = "t = 0", title_final = "t = T",
        colormap = :viridis, label = "u")
    fig = Figure(size = (900, 400))
    lo, hi = extrema(vcat(vec(initial), vec(final)))
    for (k, (ttl, f)) in enumerate(((title_initial, initial), (title_final, final)))
        ax = Axis(fig[1, k]; title = ttl, xlabel = "x₁", ylabel = k == 1 ? "x₂" : "",
            aspect = DataAspect())
        z = threshold === nothing ? f : map(v -> v < threshold ? NaN : v, f)
        hm = heatmap!(ax, x₁, x₂, z; colormap = colormap, colorrange = (lo, hi))
        contour!(ax, x₁, x₂, contours; color = :black, linewidth = 0.7, levels = 8)
        k == 2 && Colorbar(fig[1, 3], hm; label = label)
    end
    save(path, fig)
    return path
end

@doc raw"""
    figure_traces(path, traces, Sη)

The lower panels of Figs. 2 and 4: the relative energy error and the entropy against time, with
``S_\eta`` marked.

Both discretisations are drawn on the same axes. That is the point of the figure in this
reproduction — the manuscript has one method and this has two, and the reader should be able to
see whether they agree without comparing two pages.

`Sη = nothing` omits the line. Section 5.4's B3 is the case that needs it: with
``s = \omega\log\omega`` the constrained entropy minimum has no closed form, and drawing the
run's own final entropy as if it were a reference would put a result on the axis where a
prediction belongs.
"""
function figure_traces(path, traces, Sη)
    fig = Figure(size = (900, 380))
    ax1 = Axis(fig[1, 1]; xlabel = "t", ylabel = "|H(t) - H₀| / |H₀|", yscale = log10,
        title = "energy conservation")
    ax2 = Axis(fig[1, 2]; xlabel = "t", ylabel = "S", title = "entropy")
    for (label, tr) in traces
        e = energy_error(tr)
        # A log axis cannot take the exact zeros a conserved quantity produces.
        lines!(ax1, tr.t, max.(e, 1e-18); label = label)
        lines!(ax2, tr.t, tr.S; label = label)
    end
    Sη === nothing ||
        hlines!(ax2, [Sη]; color = :black, linewidth = 2, linestyle = :dash, label = "S_η")
    axislegend(ax1; position = :rb)
    axislegend(ax2; position = :rt)
    save(path, fig)
    return path
end

@doc raw"""
    figure_scatter(path, initial, final; reference = nothing, xlabel = "φ")

The scatter plots of Figs. 1, 3 and 5: the points ``(\phi_{ij}, \omega_{ij})``, at the initial
and the final time.

A functional relation emerging from the initial cloud is the manuscript's evidence that the
relaxed state is an equilibrium. `reference` draws the closed-form comparison over it — the
contour averages for Fig. 1, the fitted member of ``\mathfrak{C}_\eta`` for Fig. 5.
"""
function figure_scatter(path, initial, final; reference = nothing, xlabel = "φ",
        ylabel = "ω")
    fig = Figure(size = (560, 420))
    ax = Axis(fig[1, 1]; xlabel = xlabel, ylabel = ylabel)
    scatter!(ax, initial[1], initial[2]; markersize = 2, color = (:black, 0.25),
        label = "t = 0")
    scatter!(
        ax, final[1], final[2]; markersize = 2, color = (:crimson, 0.5), label = "t = T")
    if reference !== nothing
        scatter!(ax, reference[1], reference[2]; markersize = 9, color = :dodgerblue,
            marker = :xcross, label = "reference")
    end
    axislegend(ax; position = :lt)
    save(path, fig)
    return path
end

@doc raw"""
    figure_cone(path, traces, H₀)

The left panel of Figs. 6 and 7: the trajectory in the plane
``(S(u), 1/\|\phi\|^2_{L^2})``, with the region excluded by `eq:theoretical-limits` shaded and
the shrunk cone `eq:shrunk-cone` with ``a = 1/2`` drawn as a dashed line.

The vertex of the cone is the minimum entropy state, at ``(S_\eta, 1/2H_0)``.
"""
function figure_cone(path, traces, H₀)
    fig = Figure(size = (560, 460))
    ax = Axis(fig[1, 1]; xlabel = "S(u)", ylabel = "1/‖φ‖²")
    Smax = maximum(maximum(tr.S) for (_, tr) in traces)
    Sr = range(H₀, Smax * 1.02; length = 200)
    lower = fill(1 / (2H₀), length(Sr))
    upper = 1 / (2H₀) .+ (Sr .- H₀) ./ (2H₀^2)
    band!(ax, Sr, lower, upper; color = (:grey, 0.18))
    lines!(ax, Sr, upper; color = :grey, linewidth = 1)
    lines!(ax, Sr, lower; color = :grey, linewidth = 1)
    # eq:shrunk-cone with a = 1/2, whose slope is half the upper boundary's.
    lines!(ax, Sr, 1 / (2H₀) .+ 0.5 .* (Sr .- H₀) ./ (2H₀^2);
        color = :black, linestyle = :dash, linewidth = 1, label = "shrunk cone, a = 1/2")
    for (label, tr) in traces
        lines!(ax, tr.S, 1 ./ tr.φ²; label = label, linewidth = 2)
    end
    scatter!(ax, [H₀], [1 / (2H₀)]; color = :black, markersize = 12, label = "vertex")
    # Top-left, not bottom-right. The trajectory runs along the cone's LOWER boundary, from the
    # top-right corner down to the vertex, so a legend at `:rb` sits on top of it and hides
    # where it starts. The top left is empty because the axis begins at the vertex and the
    # cone's upper boundary *rises* from there — everything above that boundary, and hence the
    # whole top-left corner, is outside the shaded band.
    axislegend(ax; position = :lt)
    save(path, fig)
    return path
end

@doc raw"""
    figure_rates(path, t, excess, td, dist)

The right panel of Figs. 6 and 7: the excess entropy ``S(u(t)) - S_\eta`` and the distance
``\|\omega(t) - \omega(T)\|``, on a semi-log scale, with reference slopes of 1 and 1/2.

The reference slopes are the **exact** linearised rates, not fitted ones — see
`scripts/verify_projector_rates.jl`.
"""
function figure_rates(path, t, excess, td, dist)
    fig = Figure(size = (560, 420))
    ax = Axis(fig[1, 1]; xlabel = "t", ylabel = "", yscale = log10)
    lines!(ax, t, max.(excess, 1e-18); label = "S(t) - S_η", linewidth = 2)
    lines!(ax, td, max.(dist, 1e-18); label = "‖ω(t) - ω(T)‖", linewidth = 2)
    # Reference slopes anchored at a quarter of the way through, where both are still large.
    i = max(fld(length(t), 4), 1)
    lines!(ax, t, excess[i] .* exp.(-(t .- t[i])); color = :black, linestyle = :dash,
        linewidth = 1, label = "rate 1 (exact)")
    j = max(fld(length(td), 4), 1)
    lines!(ax, td, dist[j] .* exp.(-0.5 .* (td .- td[j])); color = :grey,
        linestyle = :dot, linewidth = 1, label = "rate 1/2 (exact)")
    axislegend(ax; position = :rt)
    save(path, fig)
    return path
end

@doc raw"""
    figure_tau(path, hs, measured)

The lower panel of Fig. 1: the relaxation time ``\tau_h`` against ``h`` on the contours of the
two central islands.

Where the manuscript plots only `eq:relaxation-time` evaluated from the formula, this draws the
closed form of [`relaxation_time`](@ref) as a curve and the rate **measured from the run** as
points over it, which is what makes the panel a test rather than an illustration. The curve is
evaluated here over a refined ``h``, so only the contours `hs` and their `measured` times are
passed in.
"""
function figure_tau(path, hs, measured)
    fig = Figure(size = (560, 420))
    ax = Axis(fig[1, 1]; xlabel = "h", ylabel = "τ_h", yscale = log10)
    hh = range(max(minimum(hs) / 2, 1e-3), 1.0; length = 400)
    lines!(ax, hh, [relaxation_time(h) for h in hh]; color = :black, linewidth = 2,
        label = "closed form  1/(4h·agm(1,√h)²)")
    scatter!(ax, collect(hs), collect(measured); color = :crimson, markersize = 12,
        label = "measured from the run")
    axislegend(ax; position = :rt)
    save(path, fig)
    return path
end
