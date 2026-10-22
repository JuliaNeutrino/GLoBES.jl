(* Scalar non-standard interactions and decoherence at long-baseline experiments.
   The effective Hamiltonian is built from the blocks below:
       H = M M^dagger / 2E + A Diag[{1, 0, 0}]
   with the scalar NSI correcting the mass matrix, M = U diag(m) U^dagger + Sqrt[|dm31|] eta,
   and decoherence damping the coherences of the propagation eigenstates. *)

Model`Name = "SNSI";
Model`NameLaTeX = "Scalar NSI with decoherence";
Model`Authors = "GLoBES.jl example";

(* The PDG mixing matrix R23 . Ud . R13 . Conj[Ud] . R12: a blank block takes the default *)
MixingMatrix[Vacuum]{ }

(* Absolute masses: the defaults, m1 = mlightest (NH) or m3 = mlightest (IH) *)
MassHierarchy[NH]{ }
MassHierarchy[IH]{ }

(* Scalar NSI: a Hermitian correction to the mass matrix, in units of Sqrt[|dm31|] *)
MassMatrix[Vacuum]{
    M = U . Diag[{m1, m2, m3}] . Dagger[U] + Sqrt[Abs[dm31]] eta
    eta = {{eta_ee, eta_emu Exp[I delta_emu], eta_etau Exp[I delta_etau]},
           {eta_emu Exp[-I delta_emu], eta_mumu, eta_mutau Exp[I delta_mutau]},
           {eta_etau Exp[-I delta_etau], eta_mutau Exp[-I delta_mutau], eta_tautau}}
}

(* The standard matter potential A Diag[{1, 0, 0}] is included by default *)

(* Decoherence: Gamma(E) = Gamma0 (E/GeV)^n for all pairs of eigenstates *)
Decoherence[Lindblad]{ Gamma = Gamma0 (E/GeV)^n }
