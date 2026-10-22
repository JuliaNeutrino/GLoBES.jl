# NuFIT data

The best fits, 1σ and 3σ ranges of NuFIT 5.3, 6.0 and 6.1 are built in, for both orderings, with and without the Super-Kamiokande atmospheric data.

```julia
nf = nufit("6.1")                          # the best ordering
nf = nufit("6.1"; ordering = :inverted, sk = true)
nf.s23                                     # (best, plus, minus, lo, hi)
truth = nufit_params(nf)                   # the best fit as OscParams
priors = nufit_priors(nf)                  # 1σ Gaussian priors for marginalize
nufit_ranges(nf; sigma = 3)                # (lo, hi) of each parameter
```

`nufit_params` converts the published Δm²32 of the inverted ordering to Δm²31 = Δm²32 + Δm²21.

## Δχ² projections

[`nufit_chi2`](@ref) reads the Δχ² projections that NuFIT publishes; the one-dimensional projections of NuFIT 6.0 and 6.1 are bundled with the package, and `full = true` downloads the release files once (into a scratch space) for the two- and three-dimensional ones.

```julia
t = nufit_chi2("6.1")
t(truth)                                   # Δχ², the sum of the 1D projections
x, χ² = t[:θ23]                            # a projection
fit = marginalize(test; prior = nufit_prior(t; use = (:θ12, :θ13, :Δm21)))
```

Please cite NuFIT (Esteban et al., JHEP 12 (2024) 216, and `www.nu-fit.org`) when you use these data.
