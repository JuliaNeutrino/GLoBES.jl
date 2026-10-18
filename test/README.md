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

## Known differences to C

* Rule rates outside a rule's energy window are left uninitialised by C
  (`glbGetRuleRatePtr`), so they are compared only inside the window.
* `chiSpectrumCalib` is only piecewise smooth in the energy calibration
  parameter. Its minimum can change by a few percent when the starting point
  changes in the last digits, in C as well, so the minimised value is compared
  loosely and the function itself exactly at fixed nuisance parameters.
* `NF_hR_lT` and `BBvar_WC` give χ² ≈ 1e100 (GLoBES' value for a fit rate
  ≤ 0) for some rules with the AEDL variables used here; the tests check that
  the Julia code reproduces this.
