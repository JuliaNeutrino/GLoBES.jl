# Neutrino mass mechanisms: the scotogenic loop and the type-I seesaw

@testset "neutrino mass mechanisms" begin
    mech = ScotogenicLoop()
    vals = Dict{Symbol, Any}(
        :v => 246.22,
        :lambda3 => 0.2,
        :lambda4 => -0.1,
        :lambda5 => 1e-9,
        :mEt2 => 1e6,
        :Mn => Diagonal([1e3, 2e3, 3e3]),
    )
    s = scotogenic_masses(mech, vals)
    @test s.mηp ≈ sqrt(1e6 + 0.2 * 246.22^2 / 2)

    # PMNS and masses
    p = nufit_params("6.1")
    U = pmns_matrix(p; majorana_phases = (0.3, 0.7))
    @test U' * U ≈ I
    @test abs2(U[1, 3]) ≈ sin(p.θ13)^2
    m = light_masses(p; lightest = 0.01)
    @test m[2]^2 - m[1]^2 ≈ p.Δm21 && m[3]^2 - m[1]^2 ≈ p.Δm31
    R = orthogonal_matrix(0.3 + 0.1im, 0.2, -0.4)
    @test transpose(R) * R ≈ I

    # Casas-Ibarra: the Yukawas reproduce the oscillation data
    vals[:Yn] = casas_ibarra(mech, vals, p; lightest = 1e-3)
    ν = neutrino_spectrum(mech, vals)
    @test ν.ordering === :normal
    @test ν.params.values ≈ p.values atol = 1e-10
    @test ν.m ≈ light_masses(p; lightest = 1e-3)
    @test ν.Σm ≈ sum(ν.m)

    pio = nufit_params("6.1"; ordering = :inverted)
    vals[:Yn] = casas_ibarra(mech, vals, pio; lightest = 0.01, majorana_phases = (1.0, 2.0), R)
    ν = neutrino_spectrum(mech, vals)
    @test ν.ordering === :inverted
    @test ν.params.values ≈ pio.values atol = 1e-10
    @test all(ν.majorana_phases .≈ (1.0, 2.0))

    # μ → eγ grows as Y⁴: halving λ5 doubles Y² at fixed light masses
    vals[:Yn] = casas_ibarra(mech, vals, p; lightest = 1e-3)
    br = br_lfv(mech, vals)
    @test 1e-14 < br < 1e-12
    vals[:lambda5] = 0.5e-9
    vals[:Yn] = casas_ibarra(mech, vals, p; lightest = 1e-3)
    @test br_lfv(mech, vals) ≈ 4br rtol = 1e-4

    # Type-I seesaw: Mν = -v²/2 Yᵀ M⁻¹ Y
    seesaw = TypeISeesaw()
    sv = Dict{Symbol, Any}(:v => 246.22, :Mn => Diagonal([1e14, 2e14, 3e14]))
    sv[:Yv] = casas_ibarra(seesaw, sv, p; lightest = 0.002)
    @test neutrino_spectrum(seesaw, sv).params.values ≈ p.values atol = 1e-10
end
