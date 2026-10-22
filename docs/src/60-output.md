# LaTeX and plots

## LaTeX

[`latex`](@ref) gives LaTeX for models, parameters, NuFIT tables, fit results, spectra and matrices; [`TexReport`](@ref) assembles a document.

```julia
r = TexReport("Scalar NSI at DUNE"; author = "A. Physicist")
add!(r, load_model("examples/models/SNSI"); NuFit = false)
section!(r, "Results")
add!(r, latex(fit))
figure!(r, "pmue.pdf"; caption = "Appearance probability")
write_tex("report.tex", r; compile = true)
```

## Plots

With `using Plots`, GLoBES.jl draws figures directly, in the style of journal figures ([`paper_theme!`](@ref)):

```julia
using GLoBES, Plots
paper_theme!()
plot_probability("SM" => sm, "η_eτ = 0.1" => nsi; E = range(0.5, 10, 400), L = 1300, ρ = 2.848,
                 bands = ("δ_eτ ∈ [0, 2π]" => (nsi, :δ_eτ => range(0, 2π, 37)),))
plot_events(experiment, "SM" => sm, "NSI" => nsi; rule = 1)
plot_curves(δs ./ π, "NH" => σ_nh, "IH" => σ_ih; xlabel = "δ_CP/π", ylabel = "σ", markers = true)
plot_contour(x, y, χ²)                      # 1σ, 2σ, 3σ for 2 parameters
plot_nufit(nufit_chi2("6.1"), :θ23)
plot_grid([p1 p2; p3 p4]; column_titles = ("NH", "IH"), row_titles = ("ν", "ν̄"))
```

Sensitivities come from [`mass_ordering_sensitivity`](@ref), [`cp_violation_sensitivity`](@ref) and [`chi2_profile`](@ref), and bands from [`probability_band`](@ref).
