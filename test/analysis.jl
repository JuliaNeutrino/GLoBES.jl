# Bands and sensitivities

@testset "analysis" begin
    p = nufit_params("6.1")
    E = [1.0, 2.0, 3.0]
    lo, hi = probability_band(p, :μ, :e, :δ => range(0, 2π, 25); E, L = 1300, ρ = 2.848)
    P = probability(p, :μ, :e; E, L = 1300, ρ = 2.848)
    @test all(lo .<= P .<= hi) && all(hi .- lo .> 0.01)

    setup = Setup()
    with_setup(setup) do
        cd(() -> Experiment("NOvA.glb"), DATADIR)
    end
    χ² = chi2_profile(p, :θ23 => [p.θ23 - 0.05, p.θ23, p.θ23 + 0.05]; setup)
    @test χ²[2] ≈ 0 atol = 1e-10
    @test χ²[1] > 0 && χ²[3] > 0
    σ = mass_ordering_sensitivity(p, [-π / 2]; free = (:Δm31,), setup)
    @test only(σ) > 1
    σcp = cp_violation_sensitivity(p, [-π / 2, 0.0]; free = (), orderings = :same, setup)
    @test σcp[1] > 0.5 && σcp[2] < 1e-6
end
