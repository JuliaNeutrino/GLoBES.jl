# GLoBES.jl

[![Stable Documentation](https://img.shields.io/badge/docs-stable-blue.svg)](https://arnavk23.github.io/GLoBES.jl/stable)
[![Development documentation](https://img.shields.io/badge/docs-dev-blue.svg)](https://arnavk23.github.io/GLoBES.jl/dev)
[![User manual (PDF)](https://img.shields.io/badge/manual-PDF-orange.svg)](docs/manual/globes-jl-manual.pdf)
[![Test workflow status](https://github.com/arnavk23/GLoBES.jl/actions/workflows/Test.yml/badge.svg?branch=main)](https://github.com/arnavk23/GLoBES.jl/actions/workflows/Test.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/arnavk23/GLoBES.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/arnavk23/GLoBES.jl)
[![Lint workflow Status](https://github.com/arnavk23/GLoBES.jl/actions/workflows/Lint.yml/badge.svg?branch=main)](https://github.com/arnavk23/GLoBES.jl/actions/workflows/Lint.yml?query=branch%3Amain)
[![Docs workflow Status](https://github.com/arnavk23/GLoBES.jl/actions/workflows/Docs.yml/badge.svg?branch=main)](https://github.com/arnavk23/GLoBES.jl/actions/workflows/Docs.yml?query=branch%3Amain)

A pure-Julia implementation of [GLoBES](https://www.mpi-hd.mpg.de/personalhomes/globes/), the General Long Baseline Experiment Simulator.
GLoBES.jl simulates long-baseline and reactor neutrino oscillation experiments described in the Abstract Experiment Definition Language (AEDL), computes event rates, and evaluates Δχ² functions with systematics, parameter correlations and degeneracies.
It reads the same AEDL files as GLoBES 3.0, and the correlated systematics (`sys` groups, `chiMultiExp`) and tables of any length of GLoBES 3.2, so the [DUNE configuration](examples/DUNE) works.

GLoBES.jl offers two interfaces:

- a Julian high-level API (`OscParams`, `Experiment`, `simulate!`, `chi2`, `marginalize`, `chi2_scan`, ...), and
- the `glb*` functions of the C library, so C programs can be translated line by line (with 1-based numbering).

It also fixes several bugs of C GLoBES 3.0.11, among them outdated physical constants, flux-table and background-tilt errors, undefined behaviour in the minimiser and the energy-scale shift.
Each fix is described in the manual and in [`test/README.md`](test/README.md).

## Installation

GLoBES.jl is not registered yet:

```julia
using Pkg
Pkg.add(url = "https://github.com/arnavk23/GLoBES.jl")
```

## Example

```julia
using GLoBES

truth = OscParams(θ12 = asin(sqrt(0.8)) / 2, θ13 = asin(sqrt(0.001)) / 2, θ23 = π / 4,
                  δ = π / 2, Δm21 = 7e-5, Δm31 = 2e-3)

nufact = Experiment("examples/NFstandard.glb")   # a 3000 km neutrino factory
simulate!(truth)                                 # the simulated data

test = with_param(truth, :sin²2θ13, 0.0015)
chi2(test)                                       # Δχ² with systematics: 21.66
marginalize(test; fixed = (:θ13,),               # with correlations: 2.03
            priors = (θ12 = 0.1truth.θ12, Δm21 = 0.1truth.Δm21, ρ = 0.05)).chi2
```

[`examples/tour.jl`](examples/tour.jl) continues this example with projections, the sign degeneracy and a second detector at the magic baseline.
The [`examples`](examples) folder also contains Julia ports of all C GLoBES example programs.

## Notebooks

The Jupyter notebooks in [`notebooks`](notebooks) run on Google Colab. Choose **Runtime → Change runtime type → Julia** first. Their outputs are stored, so GitHub shows them without running anything.

| Notebook | Contents | Colab |
| --- | --- | --- |
| [Getting started](notebooks/01_getting_started.ipynb) | probabilities, DUNE event rates, χ², a model of new physics in a few lines | [![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/arnavk23/GLoBES.jl/blob/main/notebooks/01_getting_started.ipynb) |
| [Scalar NSI and decoherence at DUNE](notebooks/02_scalar_NSI_and_decoherence_at_DUNE.ipynb) | the SNSI model files, probabilities, mass-ordering and CP-violation sensitivities, the bound on η_eτ | [![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/arnavk23/GLoBES.jl/blob/main/notebooks/02_scalar_NSI_and_decoherence_at_DUNE.ipynb) |

## Documentation

- The [user manual](docs/manual/globes-jl-manual.pdf) (PDF, built from [`docs/manual`](docs/manual)) follows the structure of the GLoBES 3.0 manual: the physics and formulas behind every computation, the high-level API, AEDL, the `glb*` compatibility layer and the differences from C GLoBES. All of its numbers and figures are computed with GLoBES.jl.
- The [online documentation](https://arnavk23.github.io/GLoBES.jl/dev) contains a tutorial and the docstrings of all functions.

## Tests

```shell
julia --project=. test/runtests.jl
```

The tests compare event rates, probabilities and χ² values with C GLoBES 3.0.11 (with the bug fixes of [`test/reference/c-fixes.patch`](test/reference/c-fixes.patch) applied) for the experiments in [`test/data`](test/data).

## License

GLoBES.jl is a derivative of C GLoBES and, like it, is distributed under the GNU General Public License, version 2 or later; see [LICENSE](LICENSE).

## References

- P. Huber, M. Lindner and W. Winter, [Simulation of long-baseline neutrino oscillation experiments with GLoBES](https://doi.org/10.1016/j.cpc.2005.01.003), Comput. Phys. Commun. 167, 195 (2005).
- P. Huber, J. Kopp, M. Lindner, M. Rolinec and W. Winter, [New features in the simulation of neutrino oscillation experiments with GLoBES 3.0](https://doi.org/10.1016/j.cpc.2007.05.004), Comput. Phys. Commun. 177, 432 (2007).
- [GLoBES 3.0 manual](https://www.mpi-hd.mpg.de/personalhomes/globes/documentation/globes-manual-3.0.8.pdf).
