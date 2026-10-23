# Particle models

[`load_particle_model`](@ref) reads a model directory of Mathematica-style files (the model file, `particles.m`, `parameters.m`, `SPheno.m`); [`BSMModel`](@ref) and `add_gauge!`, `add_fermion!`, `add_scalar!`, `add_lagrangian!` build one in Julia.

```julia
scoto = load_particle_model("Scotogenic/")
check_model(scoto)                # gauge and global invariance, mass dimension of each term
s = spectrum(scoto, values)       # tree-level scalar and Majorana masses, tadpoles solved
write_leshouches("LesHouches.in", scoto, values)
```

## Neutrino masses

[`ScotogenicLoop`](@ref) and [`TypeISeesaw`](@ref) give the light neutrino mass matrix of a model, and [`casas_ibarra`](@ref) the Yukawa couplings that reproduce oscillation data:

```julia
mech = ScotogenicLoop()           # parameter names as in the scotogenic model file
values[:Yn] = casas_ibarra(mech, values, nufit_params("6.1"); lightest = 1e-3)
neutrino_spectrum(mech, values)   # masses, angles, phases, mββ, Σm
br_lfv(mech, values, :μ, :e)      # BR(μ → eγ); MEG II: 1.5e-13
```
