# [Oscillation models](@id oscillation-models)

Like SARAH for particle models, GLoBES.jl reads oscillation models from files.
The effective Hamiltonian is the vacuum Hamiltonian plus additional terms, and each ingredient is a block that takes its default when it is left out or left blank.

## Files

```text
MyModel/
  MyModel.m      the model: blocks for the mixing matrix, masses and additional terms
  parameters.m   parameter values (NuFIT by default), LaTeX names, units, ranges
  particles.m    optional: Flavours = {e, mu, tau, s}
  Run.m          what to calculate, plot, fit and write to LaTeX
```

[`load_model`](@ref) reads the model, [`run_model`](@ref) executes `Run.m`.
`examples/models` has the scalar NSI with decoherence (`SNSI`) and a 3+1 sterile neutrino (`Sterile`).

## Blocks

```mathematica
MixingMatrix[Vacuum]{ U = R23 . Ud . R13 . Conj[Ud] . R12 }      (* the default *)
MixingMatrix[R23 . R13 . R12]{ R13 = IdentityMatrix[] }          (* the definition in brackets *)
MassSquared[Vacuum]{ M2 = Diag[{0, dm21, dm31}] }                 (* the default *)
MassHierarchy[NH]{ m1 = mlightest, m2 = Sqrt[m1^2 + dm21], m3 = Sqrt[m1^2 + dm31] }
MassMatrix[Vacuum]{ M = U . Diag[{m1, m2, m3}] . Dagger[U] + Sqrt[Abs[dm31]] eta }
AdditionalTerms[Matter]{ H = A Diag[{1, 0, 0}] }                  (* the default *)
AdditionalTerms[LIV]{ H = c E Diag[{1, -1, 0}] }
Decoherence[Lindblad]{ Gamma = Gamma0 (E/GeV)^n }
```

The Hamiltonian is M M†/2E with a `MassMatrix`, or U M² U†/2E without one, plus all `AdditionalTerms`; a `Hamiltonian{ H = ... }` block replaces it.
Definitions inside a block are local to it and override the defaults (`R12`, `R13`, `R23`, `Ud = PhaseMatrix[3, deltaCP]`, `PMNS`); a blank block (`MixingMatrix[Vacuum]{ }`) takes the default.
`MassHierarchy[Name]` blocks define masses; `Hierarchy = Name` in the model selects one, otherwise NH or IH follows the sign of `dm31`.

Everything is in eV: `E` is the energy, `A` the matter potential, `L` the layer length; `GeV`, `MeV`, `km` and `Degree` are units.
For antineutrinos H → H* with A → -A.
`*` and `.` both multiply matrices.

## Calculations

```julia
m = load_model("examples/models/SNSI")
p = OscParams(m; eta_etau = 0.1)               # NuFIT and parameters.m by default
use_model!(m)                                  # the probability engine for experiments
hamiltonian(m, p; E = 2.5, ρ = 2.848)          # E in GeV, ρ in g/cm³
effective_mixing(m, p; E = 2.5, ρ = 2.848)     # eigenvectors of H
effective_masses(m, p; E = 2.5, ρ = 2.848)     # 2E × eigenvalues, eV²
effective_parameters(m, p; E = 2.5, ρ = 2.848) # θ̃, δ̃, Δm̃² in matter
model_probability(m, p, :mu, :s; E = 2.0, L = 1300)
```

Matter eigenstates are numbered like the vacuum states they evolve into as the density decreases.

## Run files

```mathematica
L = 1300 km;
rho = 2.848;

Calculate[EffectiveMixingMatrix]{ E = 2.5 GeV, Parameters = {eta_etau -> 0.1} }
Plot[Pmue]{ Type = Probability, E = Range[0.5, 10, 0.025],
            Curves = {"SM" -> {}, "eta_etau = 0.1" -> {eta_etau -> 0.1}},
            Band = {"delta_etau" -> {eta_etau -> 0.1, delta_etau -> Range[0, 2 Pi, Pi/18]}} }
Fit[Eta]{ Experiment = "DUNE.glb", Scan = {eta_etau -> Range[-0.1, 0.1, 0.01]},
          Free = {theta23, deltaCP, dm31}, Priors = NuFit, Output = "eta.dat" }
LaTeX["results.tex"]{ Title = "...", NuFit = Off Print["Text added to the document."] }
```

`Calculate` takes `EffectiveMixingMatrix`, `EffectiveMasses`, `EffectiveParameters`, `Hamiltonian`, `MixingMatrix` and `Probability`; `Fit` takes `Type = Scan`, `MassOrdering` or `CPViolation`; `Plot` needs `using Plots`.
The LaTeX file has a section per model block, the parameters and the results; `Model`, `Hamiltonian`, `Definitions`, `Parameters`, `NuFit`, the block kinds, `Results`, `Fits` and `Figures` switch parts `On` or `Off`.
