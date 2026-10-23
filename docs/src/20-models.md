# [Oscillation models](@id oscillation-models)

Like SARAH for particle models, GLoBES.jl reads oscillation models from files.
The effective Hamiltonian is the vacuum Hamiltonian plus additional terms; each ingredient is one assignment and takes its default when it is left out.

## Files

```text
MyModel/
  MyModel.m      the model: the mixing matrix, masses and additional terms
  parameters.m   parameter values (NuFIT by default), LaTeX names, units, ranges
  particles.m    optional: Flavours = {e, mu, tau, s}
  Run.m          what to calculate, plot, fit and write to LaTeX
```

[`load_model`](@ref) reads the model, [`run_model`](@ref) executes `Run.m`.
`examples/models` has the scalar NSI with decoherence (`SNSI`) and a 3+1 sterile neutrino (`Sterile`).

## Model files

Model files are written like SARAH model files: one assignment per ingredient, with a comment that explains each entry.

```mathematica
Model`Name = "MyModel";                                             (* name, NameLaTeX, Authors, Date *)
MixingMatrix = R23 . Ud . R13 . Conj[Ud] . R12;                     (* the default *)
MassSquared = Diag[{0, dm21, dm31}];                                (* the default *)
MassHierarchy[[1]] = {NH, {m1 -> mlightest, m2 -> Sqrt[m1^2 + dm21], m3 -> Sqrt[m1^2 + dm31]}};  (* {name, masses} *)
MassMatrix = U . Diag[{m1, m2, m3}] . Dagger[U] + Sqrt[Abs[dm31]] eta;
AdditionalTerms[[1]] = {Matter, A Diag[{1, 0, 0}]};                 (* {name, matrix}: the default *)
AdditionalTerms[[2]] = {LIV, c E Diag[{1, -1, 0}]};
Decoherence = Gamma0 (E/GeV)^n;
```

The Hamiltonian is M M†/2E with a `MassMatrix`, or U M² U†/2E without one, plus all `AdditionalTerms`; `Hamiltonian = ...;` replaces it.
Ingredients that are left out take the defaults above; the built-in `R12`, `R13`, `R23`, `Ud = PhaseMatrix[3, deltaCP]` and `PMNS` can be redefined (`R13 = IdentityMatrix[];`).
An entry with the name of a default replaces it, and `AdditionalTerms[[1]] = {Matter, Off};` removes the matter potential.
`Hierarchy = NH;` selects a `MassHierarchy`; otherwise NH or IH follows the sign of `dm31`.

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

`Run.m` is written like SARAH's `SPheno.m`: settings, then entries `{name, {options}}` that run in the order of the file.

```mathematica
L = 1300 km;
rho = 2.848;
Calculate[[1]] = {EffectiveMixingMatrix, {E -> 2.5 GeV, Parameters -> {eta_etau -> 0.1}}};
Plot[[1]] = {Pmue, {Type -> Probability, E -> Range[0.5, 10, 0.025], Curves -> {"SM" -> {}, "NSI" -> {eta_etau -> 0.1}}}};
Fit[[1]] = {Eta, {Experiment -> "DUNE.glb", Scan -> {eta_etau -> Range[-0.1, 0.1, 0.01]}, Free -> {theta23, deltaCP, dm31}, Priors -> NuFit}};
LaTeX[[1]] = {"results.tex", {Title -> "...", NuFit -> Off, Text -> {"Text added to the document."}}};
```

`Calculate` takes `EffectiveMixingMatrix`, `EffectiveMasses`, `EffectiveParameters`, `Hamiltonian`, `MixingMatrix` and `Probability`; `Fit` takes `Type -> Scan`, `MassOrdering` or `CPViolation`; `Plot` needs `using Plots` and takes `Band -> {label -> rules with one Range}`, `File` and `Caption`.
The LaTeX file has a section per model ingredient, the parameters and the results; `Model`, `Hamiltonian`, `Definitions`, `Parameters`, `NuFit`, the ingredients, `Results`, `Fits` and `Figures` switch parts `On` or `Off`, and `Text` holds paragraphs and `Section["title"]`.
