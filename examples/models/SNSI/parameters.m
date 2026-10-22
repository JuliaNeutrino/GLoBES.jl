(* Parameter values. The standard parameters (theta12, theta13, theta23, deltaCP,
   dm21, dm31) are taken from NuFIT unless they are given here. *)

NuFit = {Version -> "6.1", Ordering -> Normal, SK -> False};

ParameterDefinitions = {
    {eta_ee,     {Value -> 0, LaTeX -> "\\eta_{ee}", Range -> {-0.2, 0.2}}},
    {eta_emu,    {Value -> 0, LaTeX -> "\\eta_{e\\mu}", Range -> {-0.3, 0.3}}},
    {eta_etau,   {Value -> 0, LaTeX -> "\\eta_{e\\tau}", Range -> {-0.19, 0.13}}},
    {eta_mumu,   {Value -> 0, LaTeX -> "\\eta_{\\mu\\mu}"}},
    {eta_mutau,  {Value -> 0, LaTeX -> "\\eta_{\\mu\\tau}"}},
    {eta_tautau, {Value -> 0, LaTeX -> "\\eta_{\\tau\\tau}"}},
    {delta_emu,  {Value -> 0, LaTeX -> "\\delta_{e\\mu}"}},
    {delta_etau, {Value -> 0, LaTeX -> "\\delta_{e\\tau}"}},
    {delta_mutau, {Value -> 0, LaTeX -> "\\delta_{\\mu\\tau}"}},
    {Gamma0,     {Value -> 0 GeV, Unit -> GeV, LaTeX -> "\\Gamma_0", Description -> "decoherence at 1 GeV"}},
    {n,          {Value -> 0, LaTeX -> "n", Description -> "power of the energy dependence of Gamma"}}
};
