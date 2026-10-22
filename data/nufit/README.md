# NuFIT Δχ² projections

One-dimensional Δχ² projections of the NuFIT global fits, cut from the release
files on <http://www.nu-fit.org>:

| file                      | release | ordering | SK atmospheric χ² table |
|---------------------------|---------|----------|-------------------------|
| `v61-TBoff-NO-1d.txt`     | 6.1     | normal   | no (IC23)               |
| `v61-TBoff-IO-1d.txt`     | 6.1     | inverted | no (IC23)               |
| `v61-TByes-NO-1d.txt`     | 6.1     | normal   | yes (IC24)              |
| `v61-TByes-IO-1d.txt`     | 6.1     | inverted | yes (IC24)              |
| `v60-*-1d.txt`            | 6.0     | as above | as above (IC19 / IC24)  |

Each file holds the `T13`, `T12`, `T23`, `DCP`, `DMS` and `DMA` projections
exactly as published. Δχ² is relative to the global minimum of both orderings.
`DMA` is Δm²31 for normal and Δm²32 for inverted ordering.
`nufit_chi2(...; full=true)` downloads the complete files with the two- and
three-dimensional projections.

If you use these tables, cite NuFIT: I. Esteban, M. C. Gonzalez-Garcia,
M. Maltoni, I. Martinez-Soler, J. P. Pinheiro and T. Schwetz, *NuFit-6.0:
updated global analysis of three-flavor neutrino oscillations*, JHEP 12 (2024)
216, arXiv:2410.05380, and NuFIT 6.1, `www.nu-fit.org`.
