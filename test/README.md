# Tests

`runtests.jl` compares GLoBES.jl against numbers produced by the C library
(GLoBES 3.0.11, with the bug fixes of `reference/c-fixes.patch`) for the
experiments in `data/`: event rates, oscillation
probabilities, χ² with and without systematics, all built-in χ² functions,
the `glbChi*` projections, a user-defined χ² function (Double Chooz, from
`example5.c`) and AEDL variables.

Run with

    julia --project=. test/runtests.jl

(`Pkg.test()` works too where a package registry is available.)

## Regenerating `reference_values.jl`

Needs a C GLoBES 3.0.11 installation with `reference/c-fixes.patch` applied
(`patch -p1 < c-fixes.patch` in the source directory, then build and install
as usual; `globes-config` on the `PATH`, or pass the include and library paths
by hand):

    cd test/reference
    gcc make_reference.c $(globes-config --include) $(globes-config --libs) -lm -o make_reference
    cd ../data
    ../reference/make_reference > ../reference_values.jl

The setups in `make_reference.c` and `runtests.jl` must be kept in sync.

## Where GLoBES.jl deliberately differs from C

The C numbers are a cross-check, not the specification: where C GLoBES is
wrong, GLoBES.jl fixes it and the tests check the fix instead. Chapter
"Differences from C GLoBES" of the manual (`docs/manual`) explains each one.

Bugs with a small fix are patched in the C reference too
(`reference/c-fixes.patch`), so the comparison stays exact:

* **Physical constants.** ħc (C: 1.9747235e-10 eV km), the matter potential
  factor √2 G_F N_A (C: 7.5e-14 eV per g/cm³), π (C: 3.1415 in the fluxes),
  Avogadro's number and the muon and electron masses are CODATA 2018 values.
* **Flux tables.** C reads the 501-point flux tables with a step of
  range/501 instead of range/500, i.e. at 1.002 times the energy, and
  tabulates the built-in fluxes only up to 500/501 of their endpoint. The
  beta beam phase space factor uses log10 instead of log.
* **Background tilt prior.** `chiSpectrumTilt` and `chiTotalRatesTilt` put the
  prior of the background tilt on x[4]/(emax-emin) instead of x[4].

The others are not patched in C:

* **Empty bins.** C's Poisson χ² returns 1e100 whenever a fit rate is ≤ 0,
  even when the bin has no true events either. `NF_hR_lT` and `BBvar_WC` have
  such bins above the beam energy, so C gives χ² ≈ 1e100 for them. Where the
  reference value is that sentinel, the tests only require a finite χ²
  (`matches_c`); "Poisson likelihood" tests the new behaviour.
* **Energy scale shift.** C's `glbShiftEnergyScale` (used by
  `chiSpectrumCalib`) interpolates linearly between bin centres, so χ² has a
  kink at zero shift in every bin, it doubles the last bin when shifting past
  the upper end, and it assumes equal bins. GLoBES.jl shifts the bin edges
  and interpolates the cumulative spectrum instead. The C comparison
  therefore leaves `chiSpectrumCalib` out and runs the setups that use it
  (NOvA, the reactors) with `chiSpectrumTilt`; "energy scale shift" and
  "chiSpectrumCalib" test the new method against exact bin integrals.
  DChooz's own χ² function uses the shift too, so it is compared to 3%.
* **Minimizer.** C's Powell routine computes its direction-replacement test
  through a macro with undefined behaviour, and stops in long curved valleys
  well above the minimum. GLoBES.jl uses the intended formula and restarts
  Powell until a whole run no longer improves χ², so fits can end at lower χ²
  than in C ("minimizer restarts").
* Rule rates outside a rule's energy window are left uninitialised by C
  (`glbGetRuleRatePtr`), so they are compared only inside the window.
