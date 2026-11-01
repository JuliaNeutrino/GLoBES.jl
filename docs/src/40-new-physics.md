# New physics in propagation

[`NewPhysics`](@ref) combines blocks that change the Hamiltonian of propagation. The blocks come in packages of their own:

| Package | Block | Hamiltonian | Parameters |
| --- | --- | --- | --- |
| [NeutrinoScalarNSI.jl](https://github.com/JuliaNeutrino/NeutrinoScalarNSI.jl) | `ScalarNSI` | (M + δM)(M + δM)†/2E, δM = √\|Δm²31\| η | `η_αβ`, `δ_αβ`, `m_lightest` |
| [NeutrinoVectorNSI.jl](https://github.com/JuliaNeutrino/NeutrinoVectorNSI.jl) | `VectorNSI` | V (diag(1, 0, 0) + ε) | `ε_αβ`, `δε_αβ` |
| [NeutrinoDecoherence.jl](https://github.com/JuliaNeutrino/NeutrinoDecoherence.jl) | `Decoherence` | damping Γ (E/E0)ⁿ of the coherences | `Γ` or `Γ21`, `Γ31`, `Γ32` |

```julia
using GLoBES, NeutrinoScalarNSI, NeutrinoDecoherence
model = NewPhysics(ScalarNSI(:eτ), Decoherence(n = 0))
use_model!(model)                                        # the probability engine
p = OscParams(model, nufit_params("6.1"); η_eτ = 0.1, Γ = 1e-23)
probability(p, :μ, :e; E = 2.5, L = 1300, ρ = 2.848)
```

Blocks from different packages combine; rates, χ² and fits use their parameters like the standard ones.
For other Hamiltonians, write a model file (see [Oscillation models](@ref oscillation-models)).

## Writing a block

A package adds its own block by subtyping [`PropagationBlock`](@ref GLoBES.PropagationBlock) and defining the methods of the interface; it then combines with the other blocks in `NewPhysics`:

```julia
using GLoBES

struct ElectronShift <: GLoBES.PropagationBlock end
GLoBES.block_parameters(::ElectronShift) = [:Δε]
GLoBES.parameter_aliases(::ElectronShift) = [:deps => :Δε]
GLoBES.matter_correction!(ε, ::ElectronShift, p, x) = (ε[1, 1] += x[:Δε]; ε)

model = NewPhysics(ElectronShift())   # or with other blocks: NewPhysics(ElectronShift(), Decoherence())
```

Blocks change the Hamiltonian through `mass_correction!` (with `uses_absolute_masses`) and `matter_correction!`, and damp coherences through `damps` and `damping_rates`; see [`PropagationBlock`](@ref GLoBES.PropagationBlock) for the full list.
