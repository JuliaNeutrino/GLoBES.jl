(* What to compute for the SNSI model: run with
       using GLoBES, Plots
       run_model("examples/models/SNSI")
   Numbers are in GeV, km and g/cm^3 unless a unit is given. *)

L = 1300 km;
rho = 2.848;

Calculate[EffectiveMixingMatrix]{ E = 2.5 GeV, Parameters = {eta_etau -> 0.1} }
Calculate[EffectiveMasses]{ E = 2.5 GeV, Parameters = {eta_etau -> 0.1} }
Calculate[EffectiveParameters]{ E = 2.5 GeV }
Calculate[Probability]{ From = mu, To = e, E = {1, 2.5, 5}, Parameters = {eta_etau -> 0.1, Gamma0 -> 1*^-23 GeV} }

Plot[Pmue]{
    Type = Probability, From = mu, To = e, E = Range[0.5, 10, 0.025],
    Curves = {"SM" -> {}, "eta_etau = 0.1" -> {eta_etau -> 0.1}, "Gamma0 = 1e-23 GeV" -> {Gamma0 -> 1*^-23 GeV}},
    Band = {"delta_etau in [0, 2 Pi]" -> {eta_etau -> 0.1, delta_etau -> Range[0, 2 Pi, Pi/18]}},
    File = "pmue.pdf", Caption = "Appearance probability at 1300 km"
}

LaTeX["SNSI.tex"]{
    Title = "Scalar NSI and decoherence at DUNE",
    Definitions = On,
    NuFit = Off,
    Print["The Hamiltonian is diagonalised numerically in each layer of constant density."]
}
