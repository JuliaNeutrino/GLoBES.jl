# Tests

`runtests.jl` compares GLoBES.jl against numbers produced by the C library
(GLoBES 3.0.11) for the experiments in `data/`: event rates, oscillation
probabilities, χ² with and without systematics, all built-in χ² functions,
the `glbChi*` projections, a user-defined χ² function (Double Chooz, from
`example5.c`) and AEDL variables.

Run with

    julia --project=. test/runtests.jl

(`Pkg.test()` works too where a package registry is available.)

## Regenerating `reference_values.jl`

Needs a C GLoBES installation (`globes-config` on the `PATH`, or pass the
include and library paths by hand):

    cd test/reference
    gcc make_reference.c $(globes-config --include) $(globes-config --libs) -lm -o make_reference
    cd ../data
    ../reference/make_reference > ../reference_values.jl

The setups in `make_reference.c` and `runtests.jl` must be kept in sync.

## Where GLoBES.jl deliberately differs from C

The C numbers are a cross-check, not the specification: where C GLoBES is
wrong, GLoBES.jl fixes it and the tests check the fix instead.

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
