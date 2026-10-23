(* The 3+1 model in GLoBES.jl: a fourth, sterile neutrino mixing with the three active ones. *)
(* The flavours are listed in particles.m; with four of them the built-in R12, R13, R23 and PMNS are 4x4 matrices. *)

Model`Name = "Sterile";
Model`NameLaTeX = "3+1 sterile neutrino";

(* Rotations with the sterile state: Rotation[i, j, angle, phase] *)

R14 = Rotation[1, 4, theta14];
R24 = Rotation[2, 4, theta24, delta24];
R34 = Rotation[3, 4, theta34, delta34];

(* Mixing matrix in vacuum *)

MixingMatrix = R34 . R24 . R14 . PMNS;

(* Masses squared, relative to m1^2 *)

MassSquared = Diag[{0, dm21, dm31, dm41}];

(* Additional terms. Each entry: {name, matrix in eV}. Charged and neutral currents: relative to the sterile
   state, the active flavours feel A for nu_e and A Nn/(2 Ne) for all, with Nn = Ne *)

AdditionalTerms[[1]] = {Matter, A Diag[{1, 0, 0, 1/2}]};
