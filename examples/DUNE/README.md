# DUNE

The GLoBES configuration of DUNE from the DUNE Technical Design Report
(DUNE collaboration, [arXiv:2103.04797](https://arxiv.org/abs/2103.04797)), as used for the study of scalar NSI and decoherence at DUNE.
It is written for GLoBES 3.2: correlated systematics (`sys(#name)` groups and `chiMultiExp`), flux and cross section tables of any length and `$emin = 0`.

Only the AEDL files are in git. The data tables (`eff/*.dat`, `flux/*.txt`, `smr/*.txt`, `xsec/*.dat`) are the ancillary files of arXiv:2103.04797; `dune_configuration()` downloads them on first use (checking each against its SHA-256) and returns the path of a ready-to-load copy:

```julia
dune = Experiment(dune_configuration())
```

The DUNE tests (`test/multiexp.jl`) do the same, so they run anywhere with network access.
