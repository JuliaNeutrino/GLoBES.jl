(* The 3+1 model: a fourth, sterile neutrino mixing with the active ones.
   U = R34 . R24 . R14 . PMNS, where the default R12, R13, R23 and PMNS are
   4x4 matrices here because the model has four flavours. *)

Model`Name = "Sterile";
Model`NameLaTeX = "3+1 sterile neutrino";

R14 = Rotation[1, 4, theta14];
R24 = Rotation[2, 4, theta24, delta24];
R34 = Rotation[3, 4, theta34, delta34];

MixingMatrix[Vacuum]{ U = R34 . R24 . R14 . PMNS }

MassSquared[Vacuum]{ M2 = Diag[{0, dm21, dm31, dm41}] }

(* Charged and neutral currents: relative to the sterile state, the active
   flavours feel A for nu_e and A Nn/(2 Ne) for all, with Nn = Ne *)
AdditionalTerms[Matter]{ H = A Diag[{1, 0, 0, 1/2}] }
