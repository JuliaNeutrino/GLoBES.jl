# Examples

Julia ports of the example programs that come with C GLoBES 3.0.11. They use
the same experiment files (copied here) and write the same output files, so
they can be compared line by line with the C versions.

| File               | C original       | What it does                                                         | Output                  |
|--------------------|------------------|----------------------------------------------------------------------|-------------------------|
| `example-tour.jl`  | `example-tour.c` | Tour of the API: probabilities, χ², projections, degeneracy, two experiments | `gl-tour.dat`  |
| `example1.jl`      | `example1.c`     | χ² in the s22θ13–δCP plane (systematics only)                        | `test1.dat`             |
| `example2.jl`      | `example2.c`     | Projection onto the s22θ13 axis: δCP only vs. all parameters free    | `test2.dat`             |
| `example3.jl`      | `example3.c`     | Locate the sgn(Δm²31) degeneracy and map χ² around it                | `test3.dat`             |
| `example4.jl`      | `example4.c`     | θ13 sensitivity bar charts: statistics, systematics, correlations, degeneracies | `test4a.dat` … `test4d.dat` |
| `example5.jl`      | `example5.c`     | Double Chooz near+far sensitivity vs. running time with user-defined χ² functions | `test5a.dat` … `test5e.dat` |
| `example6.jl`      | `example6.c`     | User-defined probability engine (wave packet decoherence, 7th parameter σ_E) | `test6a.dat`, `test6b.dat` |

Run them from any directory; the output files are written to the current
directory:

    julia --project=/path/to/GLoBES.jl /path/to/GLoBES.jl/examples/example1.jl

`myio.jl` holds the small output helpers (`InitOutput`, `AddToOutput`,
`AddToOutput2`) that the C examples take from `myio.c`.

Differences to the C code:

* Experiment and rule numbers are 1-based, as everywhere in GLoBES.jl, so
  `glbChiSys(p, 0, 0)` becomes `glbChiSys(p, 1, 1)` and the non-standard
  parameter of example 6 has index 7.
* `example5.c` uses the GSL Brent root finder; `example5.jl` contains a port
  of it (`BrentSolver`) that reproduces GSL's iterations exactly.
* The C examples load experiment files from the current directory; the Julia
  ones load them from this directory, so they can be run from anywhere.

## Agreement with C GLoBES

All outputs were compared number by number with those of the C examples
(GLoBES 3.0.11, gcc -O2), run in the same directory. GLoBES.jl fixes some
problems of the C library (see `test/README.md`), so the results are not
identical:

| Example        | Largest relative difference | Why                                                  |
|----------------|-----------------------------|------------------------------------------------------|
| `example1`     | 8e-9                        | rounding                                             |
| `example2`     | 5e-5                        | fits stop at the minimizer tolerance                 |
| `example3`     | 2e-3                        | fits (the degeneracy is found at slightly lower χ²)  |
| `example4`     | 0 (4a, 4b), 3e-5 (4c), see below (4d) | fits                                       |
| `example5`     | 4e-2                        | energy scale shift (`glbShiftEnergyScale`)           |
| `example6`     | 7e-2                        | `Reactor1.glb` uses `chiSpectrumCalib`, i.e. the shift |
| `example-tour` | 2e-3 in χ²                  | fits; every χ² minimum is as low as or lower than in C |

In `test4d.dat`, the χ² for Δm²31 = 3e-3 differs by up to 13% at small θ13
and is up to 44% lower at s22θ13 > 3e-4. The degeneracy search converges at
θ13 ≈ 0 with the same χ² in both codes, but δCP is undetermined there and
the two codes end at different values (−0.17 and 3.0). The following fits
each start from the previous one, so C and GLoBES.jl follow two different
branches in δCP; the projected χ² is really the lower of the two.
