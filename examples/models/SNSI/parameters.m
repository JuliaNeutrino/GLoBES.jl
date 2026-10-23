(* Parameter values of the SNSI model. *)
(* The standard parameters (theta12, theta13, theta23, deltaCP, dm21, dm31) are taken from NuFIT unless they are given here. *)

(* NuFIT fit to take the standard parameters from: {version, mass ordering (Normal or Inverted), with SK atmospheric data} *)

NuFit = {Version -> "6.1", Ordering -> Normal, SK -> False};

(* Each entry: {name, {options}}. Value -> default value, LaTeX -> name in LaTeX output,
   Range -> {lower, upper} for scans, Unit -> unit of Value, Description -> text in the parameter table *)

ParameterDefinitions = {
    {eta_ee,      {Value -> 0, LaTeX -> "\\eta_{ee}", Range -> {-0.2, 0.2}}},
    {eta_emu,     {Value -> 0, LaTeX -> "\\eta_{e\\mu}", Range -> {-0.3, 0.3}}},
    {eta_etau,    {Value -> 0, LaTeX -> "\\eta_{e\\tau}", Range -> {-0.19, 0.13}}},
    {eta_mumu,    {Value -> 0, LaTeX -> "\\eta_{\\mu\\mu}"}},
    {eta_mutau,   {Value -> 0, LaTeX -> "\\eta_{\\mu\\tau}"}},
    {eta_tautau,  {Value -> 0, LaTeX -> "\\eta_{\\tau\\tau}"}},
    {delta_emu,   {Value -> 0, LaTeX -> "\\delta_{e\\mu}"}},
    {delta_etau,  {Value -> 0, LaTeX -> "\\delta_{e\\tau}"}},
    {delta_mutau, {Value -> 0, LaTeX -> "\\delta_{\\mu\\tau}"}},
    {Gamma0,      {Value -> 0 GeV, Unit -> GeV, LaTeX -> "\\Gamma_0", Description -> "decoherence at 1 GeV"}},
    {n,           {Value -> 0, LaTeX -> "n", Description -> "power of the energy dependence of Gamma"}}
};
