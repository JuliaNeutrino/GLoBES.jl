# A tour of GLoBES.jl with the high-level API: the example-tour.jl program
# (a port of example-tour.c) without the glb* functions.
#
# Run with `julia --project=<path to GLoBES.jl> tour.jl`.

using GLoBES
using Printf

cd(@__DIR__)                        # the AEDL file and its data files are here

# Oscillation parameters of the simulated data
truth = OscParams(θ12=asin(sqrt(0.8)) / 2, θ13=asin(sqrt(0.001)) / 2, θ23=π / 4,
                  δ=π / 2, Δm21=7e-5, Δm31=2e-3)

# A neutrino factory with a 3000 km baseline; simulate its data
nufact = Experiment("NFstandard.glb")
simulate!(truth)

# Oscillation probabilities ν_e → ν_e, ν_μ, ν_τ at 50 GeV
for f in (:e, :μ, :τ)
    @printf("P(νe → ν%s): %.6f in vacuum, %.6f in matter\n", f,
            probability(truth, :e, f; E=50, L=3000), probability(truth, :e, f, nufact; E=50))
end

# χ² of the hypothesis sin²2θ13 = 0.0015, minimised over systematics only
test = with_param(truth, :sin²2θ13, 0.0015)
@printf("\nχ² with systematics only:           %8.4f\n", chi2(test))
@printf("from the first rule only:           %8.4f\n", chi2(test; experiment=nufact, rule=1))
@printf("sum over the rules:                 %8.4f\n",
        sum(chi2(test; rule=r) for r in eachindex(rule_names(nufact))))
@printf("statistics only:                    %8.4f\n", chi2(test; systematics=false))

# Correlations: minimise over all parameters except θ13, with 10% priors on
# the solar parameters and 5% on the matter density
priors = (θ12=0.1truth.θ12, Δm21=0.1truth.Δm21, ρ=0.05)
fit = marginalize(test; fixed=(:θ13,), priors)
@printf("\nχ² with correlations:               %8.4f\n", fit.chi2)
println("at the minimum: ", fit.params)

@printf("with δ fixed as well:               %8.4f\n",
        marginalize(test; fixed=(:θ13, :δ), priors).chi2)
@printf("with only δ and ρ free:             %8.4f\n",
        marginalize(test; free=(:δ, :ρ), priors).chi2)

# The sign degeneracy: start the fit at Δm31 < 0 and minimise over everything,
# with an additional 1/3 prior on |Δm31|
mirror = OscParams(truth; Δm31=-truth.Δm31)
deg = marginalize(mirror; priors=(priors..., Δm31=truth.Δm31 / 3), centers=mirror)
@printf("\nχ² at the inverted-hierarchy minimum: %6.4f\n", deg.chi2)
println("at ", deg.params)

# Two experiments: a second neutrino factory at the "magic" baseline
magic = Experiment("NFstandard.glb")
density_profile!(magic, average_profile(7500)...)
L, ρ = density_profile(magic)
@printf("\nmagic baseline: %g km, average density %.4f g/cm³\n", L[1], ρ[1])
simulate!(truth)
@printf("χ² of both experiments:             %8.4f\n", chi2(test))
@printf("χ² at 7500 km:                      %8.4f\n", chi2(test; experiment=magic))
fit = marginalize(test; fixed=(:θ13,), priors)
@printf("with correlations:                  %8.4f\n", fit.chi2)
