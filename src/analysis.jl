# Common analyses: probability bands, the sensitivities to the mass ordering
# and to CP violation, and Δχ² profiles of new-physics parameters.

"""
    probability_band(params, from, to, name => values; E, L, ρ=0, antineutrino=false) -> (lo, hi)
    probability_band(params, from, to, e::Experiment, name => values; E, antineutrino=false)

The range of P(ν_from → ν_to) at each energy when parameter `name` takes
all `values` (for example a phase over [0, 2π]): the envelope drawn as a
band in plots.

```julia
lo, hi = probability_band(nsi, :μ, :e, :δ_eτ => range(0, 2π, 37); E = range(0.5, 10, 200), L = 1300, ρ = 2.848)   # parameters, from, to, parameter => values; GeV, km, g/cm³
```
"""
function probability_band(params::OscParams, from, to, scan::Pair{Symbol}; kwargs...)
    curves =
        [probability(with_param(params, scan.first, v), from, to; kwargs...) for v in scan.second]
    return map(minimum, zip(curves...)) |> collect, map(maximum, zip(curves...)) |> collect
end

function probability_band(params::OscParams, from, to, e::Experiment, scan::Pair{Symbol}; kwargs...)
    curves = [
        probability(with_param(params, scan.first, v), from, to, e; kwargs...) for v in scan.second
    ]
    return map(minimum, zip(curves...)) |> collect, map(maximum, zip(curves...)) |> collect
end

# The best χ² over several starting points of a fit
function best_fit(starts; kwargs...)
    best = nothing
    for s in starts
        r = marginalize(s; kwargs...)
        (best === nothing || r.chi2 < best.chi2) && (best = r)
    end
    return best
end

# Starting points: both octants of θ23 and four values of δ
function fit_starts(p::OscParams; octants = true, phases = true)
    θ23s = octants ? unique([p.θ23, π / 2 - p.θ23]) : [p.θ23]
    δs = phases ? [p.δ, mod(p.δ + π / 2, 2π), mod(p.δ + π, 2π), mod(p.δ + 3π / 2, 2π)] : [p.δ]
    return [OscParams(p; θ23 = t, δ = d) for t in θ23s for d in δs]
end

"""
    mass_ordering_sensitivity(truth, δ_true; free=(:θ23, :δ, :Δm31), priors=(;), prior=nothing,
                              setup=current_setup()) -> Vector{Float64}

The significance √Δχ² with which the wrong mass ordering is excluded, for
each true δ_CP in `δ_true` (radians), the other true parameters as in
`truth`. The test ordering has Δm²31 → -Δm²31 + Δm²21 (the same |Δm²32|),
and χ² is minimised over `free` from several starting points (both θ23
octants, four δ). Load the experiments into `setup` first.

```julia
σ = mass_ordering_sensitivity(nufit_params("6.1"), range(-π, π, 21))   # true parameters, true δ_CP values: √Δχ² of the wrong ordering
σ = mass_ordering_sensitivity(truth, range(-π, π, 21); free = (:θ23, :δ, :Δm31), priors = nufit_priors(nufit("6.1")))   # fitted parameters, priors
```
"""
function mass_ordering_sensitivity(
    truth::OscParams,
    δ_true;
    free = (:θ23, :δ, :Δm31),
    priors = NamedTuple(),
    prior = nothing,
    setup::Setup = current_setup(),
)
    map(δ_true) do δ
        t = OscParams(truth; δ = δ)
        simulate!(t; setup)
        flipped = OscParams(t; Δm31 = -t.Δm31 + t.Δm21)
        r = best_fit(fit_starts(flipped); free, priors, prior, centers = flipped, setup)
        sqrt(max(r.chi2, 0.0))
    end
end

"""
    cp_violation_sensitivity(truth, δ_true; free=(:θ23, :Δm31), orderings=:both, priors=(;),
                             prior=nothing, setup=current_setup()) -> Vector{Float64}

The significance √Δχ² with which CP conservation (δ = 0 or π) is excluded,
for each true δ_CP in `δ_true`, minimised over `free` and, with
`orderings=:both`, over both mass orderings (`:same` for the true one only).

```julia
σ = cp_violation_sensitivity(truth, range(-π, π, 21))   # true parameters, true δ_CP values: √Δχ² of δ = 0 and π
σ = cp_violation_sensitivity(truth, range(-π, π, 21); orderings = :same)   # only the true mass ordering
```
"""
function cp_violation_sensitivity(
    truth::OscParams,
    δ_true;
    free = (:θ23, :Δm31),
    orderings = :both,
    priors = NamedTuple(),
    prior = nothing,
    setup::Setup = current_setup(),
)
    map(δ_true) do δ
        t = OscParams(truth; δ = δ)
        simulate!(t; setup)
        tests = OscParams[]
        for δ0 in (0.0, π)
            push!(tests, OscParams(t; δ = δ0))
            orderings === :both && push!(tests, OscParams(t; δ = δ0, Δm31 = -t.Δm31 + t.Δm21))
        end
        χ² = minimum(tests) do test
            best_fit(
                fit_starts(test; phases = false);
                free = Tuple(f for f in free if f !== :δ),
                priors,
                prior,
                centers = test,
                setup,
            ).chi2
        end
        sqrt(max(χ², 0.0))
    end
end

"""
    chi2_profile(truth, name => values; test=truth, free=(), priors=(;), prior=nothing,
                 setup=current_setup()) -> Vector{Float64}

Δχ² of a parameter: the data are simulated with `truth`, and the test
parameters (`test` with `name` set to each value) are fitted over `free`.
For example the sensitivity to a new-physics coupling, with the data
simulated without it.

```julia
χ² = chi2_profile(nsi, :η_eτ => range(-0.1, 0.1, 21); free = (:θ23, :δ, :Δm31))   # true parameters, scanned parameter => values, fitted parameters
```
"""
function chi2_profile(
    truth::OscParams,
    scan::Pair{Symbol};
    test::OscParams = truth,
    free = (),
    priors = NamedTuple(),
    prior = nothing,
    setup::Setup = current_setup(),
)
    simulate!(truth; setup)
    map(scan.second) do v
        p = with_param(test, scan.first, v)
        isempty(free) && prior === nothing && isempty(priors) && return chi2(p; setup)
        best_fit(
            fit_starts(p; phases = :δ in free, octants = :θ23 in free);
            free,
            priors,
            prior,
            centers = truth,
            setup,
        ).chi2
    end
end

# ---------------------------------------------------------------------------
#  Plots (implemented in the Plots.jl extension)
# ---------------------------------------------------------------------------

const PLOTS_MESSAGE = "plotting needs Plots.jl: run `using Plots` first"

"""
    plot_probability(curves...; from=:μ, to=:e, E, L, ρ=0, antineutrino=false, bands=(), kwargs...)

P(ν_from → ν_to) against energy. Each curve is `label => params`, each band
`label => (params, name => values)` (see [`probability_band`](@ref)). Other
keywords go to Plots. Needs `using Plots`.

```julia
plot_probability("SM" => sm, "η_eτ = 0.1" => nsi; E = range(0.5, 20, 400), L = 1300, ρ = 2.848)   # label => parameters; GeV, km, g/cm³
plot_probability("NSI" => nsi; E = range(0.5, 20, 400), L = 1300, ρ = 2.848, bands = ("δ_eτ ∈ [0, 2π]" => (nsi, :δ_eτ => range(0, 2π, 37)),))   # band: label => (parameters, parameter => values)
```
"""
function plot_probability end

"""
    plot_events(experiment, curves...; rule=1, kwargs...)

Event spectra of a rule as histograms, one for each `label => params`.
Needs `using Plots`.

```julia
plot_events(nova, "SM" => sm, "NSI" => nsi; rule = 1)   # experiment, label => parameters, rule
```
"""
function plot_events end

"""
    plot_curves(x, curves...; xlabel, ylabel, markers=false, xscale=:identity, kwargs...)

Lines for `label => y` (sensitivities against δ_CP, Δχ² profiles, ...), in
the paper style. Needs `using Plots`.

```julia
plot_curves(δs ./ π, "NH" => σ_nh, "IH" => σ_ih; xlabel = "δ_CP/π", ylabel = "σ", markers = true)   # x, label => y, axis labels
```
"""
function plot_curves end

"""
    plot_contour(x, y, χ²; levels=(2.30, 6.18, 11.83), labels=("1σ", "2σ", "3σ"), kwargs...)

Confidence regions of a two-parameter Δχ² (2 degrees of freedom by
default). Needs `using Plots`.

```julia
plot_contour(range(0.07, 0.11, 41), range(0, 2π, 61), grid)   # x, y, Δχ² matrix from chi2_scan: 1σ, 2σ and 3σ
plot_contour(range(0.07, 0.11, 41), range(0, 2π, 61), grid; levels = (4.61,), labels = ("90% CL",))   # other levels
```
"""
function plot_contour end

"""
    plot_nufit(t::NuFitChi2, names...; ordering=:both, kwargs...)

A NuFIT Δχ² projection: one parameter (a profile) or two (contours).
Needs `using Plots`.

```julia
plot_nufit(nufit_chi2("6.1"), :θ23)   # Δχ² tables, one parameter: a profile
plot_nufit(nufit_chi2("6.1"; full = true), :θ23, :δ)   # two parameters: contours
```
"""
function plot_nufit end

"""
    plot_grid(panels; title="", column_titles=(), row_titles=(), legend_panel=nothing, size=...)

Arrange plots in a grid with a title, column titles and row titles, like
the multi-panel figures of papers. `panels` is a matrix of plots. Needs
`using Plots`.

```julia
plot_grid([p1 p2; p3 p4]; title = "P(νμ → νe)", column_titles = ("NH", "IH"), row_titles = ("ν", "ν̄"))   # matrix of plots, titles
```
"""
function plot_grid end

"""
    paper_theme!()

Plots defaults for publication figures: framed axes with ticks inside and
on all sides, minor ticks, a faint grid and LaTeX-like fonts. Needs
`using Plots`.

```julia
paper_theme!()   # Plots defaults for journal figures, for all plots that follow
```
"""
function paper_theme! end
