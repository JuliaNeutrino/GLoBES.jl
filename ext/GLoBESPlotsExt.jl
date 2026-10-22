# Plots for GLoBES.jl, loaded with `using Plots`: probabilities with bands,
# event spectra, sensitivity curves, Δχ² contours, NuFIT projections and
# multi-panel figures in the style of journal plots.
module GLoBESPlotsExt

using GLoBES
using Plots
import GLoBES:
    plot_probability,
    plot_events,
    plot_curves,
    plot_contour,
    plot_nufit,
    plot_grid,
    paper_theme!,
    savefig_

# Line styles of the curves, in order: SM first, then the new physics
const STYLES = (
    (:black, :solid),
    (RGB(0.0, 0.70, 0.15), :dashdot),
    (RGB(0.90, 0.10, 0.10), :dash),
    (RGB(0.10, 0.20, 0.90), :dot),
    (RGB(0.95, 0.55, 0.0), :dashdotdot),
    (RGB(0.55, 0.20, 0.70), :solid),
)
const BAND_COLORS = (RGB(1.0, 0.35, 1.0), RGB(0.45, 0.95, 0.45), RGB(0.55, 0.75, 1.0))

style(i) = STYLES[mod1(i, length(STYLES))]

function paper_theme!()
    default(
        framestyle = :box,
        grid = true,
        gridalpha = 0.12,
        gridstyle = :dot,
        minorticks = true,
        tick_direction = :in,
        guidefontsize = 11,
        tickfontsize = 9,
        legendfontsize = 8,
        titlefontsize = 10,
        linewidth = 2,
        foreground_color_legend = nothing,
        background_color_legend = RGBA(1, 1, 1, 0.85),
        size = (600, 450),
        dpi = 200,
        left_margin = 3Plots.mm,
        bottom_margin = 3Plots.mm,
    )
    return nothing
end

savefig_(plt, path) = savefig(plt, path)

# Labels are Unicode; pass LaTeXStrings as `xlabel`, `ylabel` for LaTeX
flavour_sub(f) = (f in (:μ, :mu, 2)) ? "μ" : (f in (:τ, :tau, 3)) ? "τ" : "e"

function probability_label(from, to, anti)
    a, b = flavour_sub(from), flavour_sub(to)
    nu = anti ? "ν̄" : "ν"
    "P($(nu)$(a) → $(nu)$(b))"
end

function plot_probability(
    curves::Pair...;
    from = :μ,
    to = :e,
    E,
    L = nothing,
    ρ = 0.0,
    antineutrino::Bool = false,
    bands = (),
    experiment = nothing,
    setup = current_setup(),
    xlabel = "Eν (GeV)",
    ylabel = probability_label(from, to, antineutrino),
    kwargs...,
)
    E = collect(Float64, E)
    plt = plot(; xlabel, ylabel, legend = :topright, kwargs...)
    prob(p) =
        experiment === nothing ? probability(p, from, to; E, L, ρ, antineutrino, setup) :
        probability(p, from, to, experiment; E, antineutrino)
    for (i, (label, (p, scan))) in enumerate(bands)
        lo, hi =
            experiment === nothing ?
            probability_band(p, from, to, scan; E, L, ρ, antineutrino, setup) :
            probability_band(p, from, to, experiment, scan; E, antineutrino)
        c = BAND_COLORS[mod1(i, length(BAND_COLORS))]
        plot!(
            plt,
            E,
            lo;
            fillrange = hi,
            fillcolor = c,
            fillalpha = 0.85,
            linealpha = 0,
            label = string(label),
        )
    end
    for (i, (label, p)) in enumerate(curves)
        c, ls = style(i)
        plot!(plt, E, prob(p); label = string(label), color = c, linestyle = ls, linewidth = 2)
    end
    return plt
end

# Bin edges of an experiment's (reconstructed) energy bins
function bin_edges(e)
    c, w = bin_centers(e), bin_widths(e)
    return [c .- w ./ 2; c[end] + w[end] / 2]
end

function plot_events(
    e::Experiment,
    curves::Pair...;
    rule = 1,
    xlabel = "Reconstructed Eν (GeV)",
    ylabel = "Events/bin",
    kwargs...,
)
    edges = bin_edges(e)
    plt = plot(; xlabel, ylabel, legend = :topright, kwargs...)
    for (i, (label, p)) in enumerate(curves)
        r = event_rates(e, p; rule).total
        c, ls = style(i)
        plot!(
            plt,
            edges,
            [r; r[end]];
            seriestype = :steppost,
            label = string(label),
            color = c,
            linestyle = ls,
            linewidth = 2,
        )
    end
    return plt
end

function plot_curves(
    x,
    curves::Pair...;
    xlabel = "",
    ylabel = "",
    markers = false,
    xscale = :identity,
    kwargs...,
)
    plt = plot(; xlabel, ylabel, xscale, legend = :best, kwargs...)
    for (i, (label, y)) in enumerate(curves)
        c, ls = style(i)
        m = markers === true || (markers isa AbstractVector && i in markers)
        plot!(
            plt,
            x,
            y;
            label = string(label),
            color = c,
            linestyle = ls,
            linewidth = 2,
            markershape = m ? :circle : :none,
            markersize = 4,
            markerstrokewidth = 0,
            markercolor = c,
        )
    end
    return plt
end

function plot_contour(
    x,
    y,
    χ²::AbstractMatrix;
    levels = (2.30, 6.18, 11.83),
    labels = ("1σ", "2σ", "3σ"),
    colors = (:black, RGB(0.90, 0.10, 0.10), RGB(0.10, 0.20, 0.90)),
    best = true,
    plt = nothing,
    kwargs...,
)
    plt = plt === nothing ? plot(; legend = :topright, kwargs...) : plt
    z = permutedims(χ² .- minimum(χ²))
    for (k, lvl) in enumerate(levels)
        c = colors[mod1(k, length(colors))]
        contour!(
            plt,
            x,
            y,
            z;
            levels = [lvl],
            color = c,
            linewidth = 2,
            colorbar = false,
            label = "",
        )
        plot!(plt, [NaN], [NaN]; color = c, linewidth = 2, label = labels[k])       # legend entry
    end
    if best
        i = argmin(χ²)
        scatter!(
            plt,
            [x[i[1]]],
            [y[i[2]]];
            marker = :star5,
            color = :black,
            markersize = 6,
            label = "",
        )
    end
    return plt
end

axis_label(n::Symbol) =
    n === :θ12 ? "sin²θ₁₂" :
    n === :θ13 ? "sin²θ₁₃" :
    n === :θ23 ? "sin²θ₂₃" :
    n === :δ ? "δCP (°)" : n === :Δm21 ? "Δm²₂₁ (eV²)" : n === :Δm31 ? "Δm²₃ₗ (eV²)" : string(n)

function plot_nufit(t::GLoBES.NuFitChi2, names::Symbol...; ordering = :both, kwargs...)
    orderings = ordering === :both ? (:normal, :inverted) : (ordering,)
    if length(names) == 1
        plt = plot(; xlabel = axis_label(names[1]), ylabel = "Δχ²", legend = :top, kwargs...)
        for (i, o) in enumerate(orderings)
            x, χ² = getindex(t, names[1]; ordering = o)
            c, ls = style(i)
            plot!(
                plt,
                x,
                χ²;
                label = "NuFIT $(t.version) $(o === :normal ? "NO" : "IO")",
                color = c,
                linestyle = ls,
            )
        end
        return plt
    end
    plt = plot(; xlabel = axis_label(names[1]), ylabel = axis_label(names[2]), kwargs...)
    for o in orderings
        x, y, χ² = getindex(t, names...; ordering = o)
        plot_contour(
            x,
            y,
            χ²;
            plt,
            best = o === orderings[1],
            labels = ("1σ $(o === :normal ? "NO" : "IO")", "2σ", "3σ"),
        )
    end
    return plt
end

function plot_grid(
    panels::AbstractMatrix;
    title = "",
    column_titles = (),
    row_titles = (),
    size = (
        330 * Base.size(panels, 2) + (isempty(row_titles) ? 0 : 60),
        270 * Base.size(panels, 1) + (isempty(title) ? 0 : 60),
    ),
    kwargs...,
)
    nr, nc = Base.size(panels)
    cells = Any[]
    for i = 1:nr
        if !isempty(row_titles)
            push!(
                cells,
                plot(;
                    framestyle = :none,
                    grid = false,
                    ticks = false,
                    annotations = (0.5, 0.5, Plots.text(row_titles[i], 11, rotation = 90)),
                ),
            )
        end
        for j = 1:nc
            p = panels[i, j]
            i == 1 && j <= length(column_titles) && plot!(p; title = column_titles[j])
            push!(cells, p)
        end
    end
    ncol = nc + (isempty(row_titles) ? 0 : 1)
    widths = isempty(row_titles) ? fill(1 / nc, nc) : [0.04; fill(0.96 / nc, nc)]
    body = plot(cells...; layout = grid(nr, ncol; widths), kwargs...)
    isempty(title) && return plot(body; size)
    header = plot(;
        framestyle = :none,
        grid = false,
        ticks = false,
        annotations = (0.5, 0.5, Plots.text(title, 12)),
    )
    return plot(header, body; layout = grid(2, 1; heights = [0.06, 0.94]), size)
end

end # module
