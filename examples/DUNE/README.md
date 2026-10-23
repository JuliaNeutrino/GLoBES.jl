# DUNE

The GLoBES configuration of DUNE from the DUNE Technical Design Report
(DUNE collaboration, [arXiv:2103.04797](https://arxiv.org/abs/2103.04797)), as used for the study of scalar NSI and decoherence at DUNE.
It is written for GLoBES 3.2: correlated systematics (`sys(#name)` groups and `chiMultiExp`), flux and cross section tables of any length and `$emin = 0`.

```julia
dune = Experiment("examples/DUNE/DUNE_GLoBES.glb")   # 40 kt, 1284.9 km, 6.5 + 6.5 years
```
