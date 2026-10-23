(* What to compute for the SNSI model. Run with
       using GLoBES, Plots
       run_model("examples/models/SNSI")
   Numbers are in GeV, km and g/cm^3 unless a unit is given. The entries run in the order of this file. *)

(* Settings used by every entry that does not give its own: baseline and matter density *)

L = 1300 km;
rho = 2.848;

(* Calculations. Each entry: {quantity, {options}}. The quantities are EffectiveMixingMatrix, EffectiveMasses,
   EffectiveParameters, Hamiltonian, MixingMatrix and Probability; the options E, rho, L, From, To, Antineutrino
   and Parameters (rules that change the defaults of parameters.m) *)

Calculate[[1]] = {EffectiveMixingMatrix, {E -> 2.5 GeV, Parameters -> {eta_etau -> 0.1}}};
Calculate[[2]] = {EffectiveMasses, {E -> 2.5 GeV, Parameters -> {eta_etau -> 0.1}}};
Calculate[[3]] = {EffectiveParameters, {E -> 2.5 GeV}};
Calculate[[4]] = {Probability, {From -> mu, To -> e, E -> {1, 2.5, 5}, Parameters -> {eta_etau -> 0.1, Gamma0 -> 1*^-23 GeV}}};

(* Plots (they need using Plots). Each entry: {name, {options}}. Curves: {label -> parameter rules},
   Band: {label -> rules with one parameter over a Range}; File and Caption go to the LaTeX report *)

Plot[[1]] = {Pmue, {
    Type -> Probability, From -> mu, To -> e, E -> Range[0.5, 10, 0.025],
    Curves -> {"SM" -> {}, "eta_etau = 0.1" -> {eta_etau -> 0.1}, "Gamma0 = 1e-23 GeV" -> {Gamma0 -> 1*^-23 GeV}},
    Band -> {"delta_etau in [0, 2 Pi]" -> {eta_etau -> 0.1, delta_etau -> Range[0, 2 Pi, Pi/18]}},
    File -> "pmue.pdf", Caption -> "Appearance probability at 1300 km"}};

(* LaTeX report. Entry: {file, {options}}. Parts of the model are switched On or Off by name;
   Text: paragraphs and Section["title"] added at the end *)

LaTeX[[1]] = {"SNSI.tex", {
    Title -> "Scalar NSI and decoherence at DUNE", Definitions -> On, NuFit -> Off,
    Text -> {"The Hamiltonian is diagonalised numerically in each layer of constant density."}}};
