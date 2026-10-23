(* Scalar non-standard interactions (NSI) and decoherence at long-baseline experiments, in GLoBES.jl. *)
(* Main model file: the ingredients of the effective Hamiltonian
       H = M M^dagger / 2E + A Diag[{1, 0, 0}]
   and the damping of the coherences by decoherence. Every ingredient that is left out takes its default. *)

Model`Name = "SNSI";
Model`NameLaTeX = "Scalar NSI with decoherence";
Model`Authors = "GLoBES.jl example";
Model`Date = "2026-10-24";

(* Mixing matrix in vacuum. The default is the PDG parametrisation, written with the built-in rotations
   R12 = Rotation[1, 2, theta12], R13, R23 and Ud = PhaseMatrix[3, deltaCP]; it is repeated here to show the form. *)

MixingMatrix = R23 . Ud . R13 . Conj[Ud] . R12;

(* Absolute masses. Each MassHierarchy entry: {name, {masses as rules}}. The defaults are
   NH: m1 -> mlightest and IH: m3 -> mlightest, the others from dm21 and dm31; the hierarchy follows the sign of dm31.
   Both are left out here, so the defaults are used. *)

(* Scalar NSI: a Hermitian correction to the mass matrix, in units of Sqrt[|dm31|].
   eta_ab is the size of the coupling of flavours a and b, delta_ab its phase. *)

eta = {{eta_ee, eta_emu Exp[I delta_emu], eta_etau Exp[I delta_etau]},
       {eta_emu Exp[-I delta_emu], eta_mumu, eta_mutau Exp[I delta_mutau]},
       {eta_etau Exp[-I delta_etau], eta_mutau Exp[-I delta_mutau], eta_tautau}};

(* Mass matrix: with one, the vacuum Hamiltonian is M M^dagger / 2E instead of U Diag[{0, dm21, dm31}] U^dagger / 2E *)

MassMatrix = U . Diag[{m1, m2, m3}] . Dagger[U] + Sqrt[Abs[dm31]] eta;

(* Additional terms of the Hamiltonian. Each entry: {name, matrix in eV}; A is the matter potential.
   The standard matter potential is the default, so this line could be left out. *)

AdditionalTerms[[1]] = {Matter, A Diag[{1, 0, 0}]};

(* Decoherence: the damping Gamma(E) of the coherences of all pairs of propagation eigenstates, in eV *)

Decoherence = Gamma0 (E/GeV)^n;
