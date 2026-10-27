```@meta
CurrentModule = GLoBES
```

# GLoBES.jl

A pure-Julia implementation of [GLoBES](https://www.mpi-hd.mpg.de/personalhomes/globes/), the General Long Baseline Experiment Simulator.
GLoBES.jl simulates long-baseline and reactor neutrino oscillation experiments described in the Abstract Experiment Definition Language (AEDL), computes event rates, and evaluates Δχ² functions with systematics, parameter correlations and degeneracies.

- The [Tutorial](@ref) introduces the high-level API.
- The [Reference](@ref reference) lists the docstrings of all functions, including the `glb*` functions of the C API.
- The [user manual (PDF)](assets/globes-jl-manual.pdf) explains the physics and formulas behind every computation, describes AEDL, and lists where GLoBES.jl deliberately differs from C GLoBES.

## Installation

```julia
using Pkg
Pkg.add(url = "https://github.com/JuliaNeutrino/GLoBES.jl")
```

## Conventions

Energies are in GeV, lengths in km, densities in g/cm³, mass squared differences in eV² and angles in radians.
Experiments, rules and channels are numbered from 1, both in the high-level API and in the `glb*` layer.
