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
(GLoBES 3.0.11, gcc -O2), run in the same directory:

| Example        | Largest relative difference                                  |
|----------------|--------------------------------------------------------------|
| `example1`     | 8e-9 (χ² with systematics only)                              |
| `example2`     | 8e-5 (χ² with correlations)                                  |
| `example3`     | 2e-4 at the degeneracy, 3e-4 in the χ² map around it         |
| `example4`     | 0 (4a, 4b), 1e-5 (4c), 2e-6 (4d)                             |
| `example5`     | 0: all roots identical to printf precision                   |
| `example6`     | 0                                                            |
| `example-tour` | χ² values ≤ 2e-5; minimum positions ≤ 2e-4; iteration counts differ by a few % |

The fits with correlations stop where the minimizer's relative tolerance
(1e-5) is reached, and the last few digits of the χ² values differ at the
1e-11 level between the two codes (libm, summation order), so minimum
positions along flat directions and iteration counts are not bit-identical.
