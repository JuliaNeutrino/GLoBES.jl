# NuFIT global fits: tables of best fits and ranges, and Δχ² projections

@testset "NuFIT" begin
    nf = nufit("6.1"; ordering = :normal)
    @test nf.s12.best == 0.3088 && nf.dcp.best == 207.0
    p = nufit_params(nf)
    @test sin(p.θ12)^2 ≈ 0.3088 && sin(p.θ13)^2 ≈ 0.02249 && sin(p.θ23)^2 ≈ 0.470
    @test p.δ ≈ deg2rad(207) && p.Δm21 ≈ 7.537e-5 && p.Δm31 ≈ 2.521e-3
    # Inverted: NuFIT gives Δm²32; Δm²31 = Δm²32 + Δm²21
    pio = nufit_params("6.1"; ordering = :inverted)
    @test pio.Δm31 ≈ -2.500e-3 + 7.537e-5
    @test nufit("6.1").ordering === :normal               # the best ordering
    @test nufit("5.3"; ordering = :normal).s23.best ≈ 0.572 rtol = 0.05
    @test_throws ArgumentError nufit("1.0")

    pr = nufit_priors(nf)
    @test pr.θ12 ≈ deg2rad((0.42 + 0.41) / 2)
    r = nufit_ranges(nf)
    @test r.θ12[1] ≈ deg2rad(32.54) && r.θ12[2] ≈ deg2rad(35.03)

    # Δχ² from the published one-dimensional projections
    t = nufit_chi2("6.1")
    @test t(p) < 0.1
    @test t(pio) ≈ 1.5 atol = 0.1
    @test t(OscParams(p; θ23 = asin(sqrt(0.6)))) > 4
    x, χ² = t[:θ23]
    @test minimum(χ²) < 0.05
    prior = nufit_prior(t; use = (:θ23,))
    @test prior(p) ≈ t(p; use = (:θ23,)) - min(t.offset[:normal], t.offset[:inverted])
end
