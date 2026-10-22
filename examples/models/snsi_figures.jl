# Probability figures of scalar NSI and decoherence at DUNE (L = 1300 km),
# drawn directly with Plots: run from the repository root with
#     julia --project=. examples/models/snsi_figures.jl
# (Plots must be in the active environment or the global one).
using GLoBES, Plots
paper_theme!()

m = load_model(joinpath(@__DIR__, "SNSI"))
use_model!(m)
E = range(0.5, 10, 400)
setting = (; L = 1300, ρ = 2.848)

panel(sm, nsi, deco, both; anti = false, from = :μ, to = :e) = plot_probability(
    "SM" => sm,
    "η_eτ = 0.1" => nsi,
    "Γ₀ = 10⁻²³ GeV" => deco,
    "η + Γ" => both;
    from,
    to,
    E,
    setting...,
    antineutrino = anti,
    bands = ("δ_eτ ∈ [0, 2π]" => (nsi, :δ_eτ => range(0, 2π, 37)),),
)

panels = map(Iterators.product((:normal, :inverted), (false, true))) do (ordering, anti)
    base = nufit_params(nufit("6.1"; ordering))
    sm = OscParams(m, base)
    nsi = OscParams(m, base; eta_etau = 0.1)
    deco = OscParams(m, base; Gamma0 = 1e-23 * 1e9)       # Γ in eV
    both = OscParams(m, base; eta_etau = 0.1, Gamma0 = 1e-23 * 1e9)
    panel(sm, nsi, deco, both; anti)
end
fig = plot_grid(
    permutedims(panels);
    title = "P(νμ → νe) at 1300 km",
    column_titles = ("Normal ordering", "Inverted ordering"),
    row_titles = ("ν", "ν̄"),
)
savefig(fig, "snsi_probabilities.pdf")
println("written snsi_probabilities.pdf")
