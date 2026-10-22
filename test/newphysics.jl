# Propagation with new physics: scalar NSI, matter NSI and decoherence,
# compared with the C engine of SNSI-Decoherence-DUNE (reference/)

# P(μ → e), P(μ → μ) from snsi_reference_driver.c: config, sign of CP, E, P_μe, P_μμ
const SNSI_REFERENCE = [
    parse.(Float64, split(l)) for
    l in eachline(joinpath(@__DIR__, "reference", "snsi_c_reference.txt")) if !startswith(l, "#")
]

snsi_base(config) =
    config < 3 ?
    OscParams(
        θ12 = asin(sqrt(0.3088)),
        θ13 = asin(sqrt(0.02249)),
        θ23 = asin(sqrt(0.470)),
        δ = deg2rad(207),
        Δm21 = 7.537e-5,
        Δm31 = 2.521e-3,
    ) :
    OscParams(
        θ12 = asin(sqrt(0.3088)),
        θ13 = asin(sqrt(0.02261)),
        θ23 = asin(sqrt(0.555)),
        δ = deg2rad(283),
        Δm21 = 7.537e-5,
        Δm31 = -2.5e-3 + 7.537e-5,
    )

@testset "new physics engine" begin
    np = NewPhysics(ScalarNSI(:eμ, :eτ), Decoherence())
    @test :η_eτ in np.names && :δ_eμ in np.names && :Γ in np.names && :m_lightest in np.names
    E = [0.5, 1.0, 2.5, 5.0]

    # No new physics: the standard engine
    p = nufit_params("6.1")
    standard =
        [probability(p, :μ, :e; E, L = 1300, ρ = 2.848, antineutrino = a) for a in (false, true)]
    use_model!(np)
    q = OscParams(np, p)
    @test probability(q, :μ, :e; E, L = 1300, ρ = 2.848) ≈ standard[1] rtol = 1e-12
    @test probability(q, :μ, :e; E, L = 1300, ρ = 2.848, antineutrino = true) ≈ standard[2] rtol =
        1e-12

    # The C engine
    couplings = Dict(
        1 => (η_eτ = 0.1, Γ = 1e-23),
        2 => (η_eμ = 0.05, δ_eμ = 1.0),
        3 => (η_eτ = -0.05, δ_eτ = π / 3, Γ = 5e-24),
    )
    for (config, cp, E, pμe, pμμ) in SNSI_REFERENCE
        q = OscParams(np, snsi_base(Int(config)); couplings[Int(config)]...)
        P = probability(q, :μ, :e; E = [E], L = 1300, ρ = 2.848, antineutrino = cp < 0)[1]
        @test P ≈ pμe rtol = 1e-10
        @test probability(q, :μ, :μ; E = [E], L = 1300, ρ = 2.848, antineutrino = cp < 0)[1] ≈ pμμ rtol =
            1e-10
    end

    # Matter NSI: ε_ee = x is the standard potential times 1 + x
    mn = NewPhysics(MatterNSI(:ee, :μτ))
    use_model!(mn)
    @test probability(OscParams(mn, p), :μ, :e; E, L = 1300, ρ = 2.848) ≈ standard[1] rtol = 1e-12
    nsi = probability(OscParams(mn, p; ε_ee = 0.5), :μ, :e; E, L = 1300, ρ = 2.848)
    standard_engine!()
    @test nsi ≈ probability(p, :μ, :e; E, L = 1300, ρ = 1.5 * 2.848) rtol = 1e-10
    use_model!(mn)

    # Decoherence damps the oscillations to the average
    dec = NewPhysics(Decoherence())
    use_model!(dec)
    d = OscParams(dec, p; Γ = 1e-18)
    U = GLoBES.pmns_matrix(p)
    average = sum(abs2(U[2, i]) * abs2(U[1, i]) for i = 1:3)
    @test probability(d, :μ, :e; E = [1.0], L = 1300, ρ = 0.0)[1] ≈ average rtol = 1e-6
    standard_engine!()
end
