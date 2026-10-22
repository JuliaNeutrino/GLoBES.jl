# New physics in propagation

[`NewPhysics`](@ref) combines blocks that change the Hamiltonian of propagation:

| Block | Hamiltonian | Parameters |
| --- | --- | --- |
| [`ScalarNSI`](@ref) | (M + δM)(M + δM)†/2E, δM = √\|Δm²31\| η | `η_αβ`, `δ_αβ`, `m_lightest` |
| [`MatterNSI`](@ref) | V (diag(1, 0, 0) + ε) | `ε_αβ`, `δε_αβ` |
| [`Decoherence`](@ref) | damping Γ (E/E0)ⁿ of the coherences | `Γ` or `Γ21`, `Γ31`, `Γ32` |

```julia
model = NewPhysics(ScalarNSI(:eτ), Decoherence(n = 0))
use_model!(model)                                        # the probability engine
p = OscParams(model, nufit_params("6.1"); η_eτ = 0.1, Γ = 1e-23)
probability(p, :μ, :e; E = 2.5, L = 1300, ρ = 2.848)
```

The scalar NSI and decoherence agree with the C engine of the SNSI-decoherence DUNE study to 10⁻¹² (see `test/reference`, built with CODATA constants).
