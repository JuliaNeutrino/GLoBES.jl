# Tutorial

This tutorial follows [`examples/tour.jl`](https://github.com/JuliaNeutrino/GLoBES.jl/blob/main/examples/tour.jl), a translation of the tour of the C GLoBES manual.
Run it from the `examples` directory, where the AEDL file `NFstandard.glb` and its data files are.

## Parameters and simulated data

Oscillation parameters are [`OscParams`](@ref); they are immutable, and [`with_param`](@ref) returns a copy with one parameter changed.
Derived quantities such as `sin²2θ13` can be read and set by name.

```julia
using GLoBES

truth = OscParams(θ12 = asin(sqrt(0.8)) / 2, θ13 = asin(sqrt(0.001)) / 2, θ23 = π / 4,
                  δ = π / 2, Δm21 = 7e-5, Δm31 = 2e-3)
nufact = Experiment("NFstandard.glb")   # a neutrino factory with a 3000 km baseline
simulate!(truth)                        # compute the simulated data
```

## Oscillation probabilities

```julia
probability(truth, :e, :μ; E = 50, L = 3000)   # vacuum
probability(truth, :e, :μ, nufact; E = 50)     # along the experiment's matter profile
```

## Δχ² with systematics

[`chi2`](@ref) compares the simulated data with the prediction for a test hypothesis and minimises over the systematic nuisance parameters:

```julia
test = with_param(truth, :sin²2θ13, 0.0015)
chi2(test)                                  # 21.66
chi2(test; experiment = nufact, rule = 1)   # 20.45, first rule only
chi2(test; systematics = false)             # 35.18, statistics only
```

## Correlations and degeneracies

[`marginalize`](@ref) minimises the Δχ² over the oscillation parameters that are not fixed, with Gaussian priors:

```julia
priors = (θ12 = 0.1truth.θ12, Δm21 = 0.1truth.Δm21, ρ = 0.05)
fit = marginalize(test; fixed = (:θ13,), priors)
fit.chi2      # 2.03
fit.params    # the parameters at the minimum
```

The minimiser finds the minimum nearest to its starting point.
To find the sign degeneracy, start at inverted ordering:

```julia
mirror = OscParams(truth; Δm31 = -truth.Δm31)
marginalize(mirror; priors = (priors..., Δm31 = truth.Δm31 / 3), centers = mirror).chi2   # 4.56
```

## Several experiments

All experiments of a setup contribute to the Δχ².
A second neutrino factory at the magic baseline:

```julia
magic = Experiment("NFstandard.glb")
density_profile!(magic, average_profile(7500)...)
simulate!(truth)
chi2(test)                                     # 28.22, both experiments
marginalize(test; fixed = (:θ13,), priors).chi2   # 3.58
```

## Scans

[`chi2_scan`](@ref) evaluates the Δχ² (or its projection) on a grid:

```julia
chi2_scan(truth, :sin²2θ13 => exp10.(-4:0.04:-2), :δ => deg2rad.(0:4:200))
```
