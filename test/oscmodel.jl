# Oscillation models defined in model files (load_model)

const MODELS = joinpath(@__DIR__, "..", "examples", "models")

@testset "oscillation models" begin
    # The standard model: every block from the defaults
    std = load_model()
    @test GLoBES.nflavours(std) == 3
    @test [p.name for p in model_parameters(std)][1:6] == [:theta12, :theta13, :theta23, :deltaCP, :dm21, :dm31]
    p = OscParams(std)
    @test p.values ≈ [nufit_params("6.1").values; 0.0]
    E = [0.5, 1.0, 2.5, 5.0]
    use_model!(std)
    a = [probability(p, :μ, :e; E, L = 1300, ρ = 2.848, antineutrino = x) for x in (false, true)]
    standard_engine!()
    b = [
        probability(nufit_params("6.1"), :μ, :e; E, L = 1300, ρ = 2.848, antineutrino = x) for
        x in (false, true)
    ]
    @test a ≈ b rtol = 1e-12

    # Hamiltonian, effective mixing and masses
    H = hamiltonian(std, p; E = 2.5, ρ = 2.848)
    @test H ≈ H'
    U = GLoBES.pmns_matrix(nufit_params("6.1"))
    A = 2.848 * GLoBES.GLB_V_FACTOR * GLoBES.GLB_Ne_MANTLE
    @test H ≈ U * Diagonal([0, p.Δm21, p.Δm31]) * U' / 5e9 + A * Diagonal([1, 0, 0]) rtol = 1e-12
    @test hamiltonian(std, p; E = 2.5, ρ = 2.848, antineutrino = true) ≈
          conj(U * Diagonal([0, p.Δm21, p.Δm31]) * U' / 5e9 - A * Diagonal([1, 0, 0])) rtol = 1e-12
    q = effective_parameters(std, p; E = 2.5, ρ = 0.0)          # vacuum
    @test q.values[1:6] ≈ p.values[1:6] rtol = 1e-8
    @test effective_masses(std, p; E = 2.5) ≈ [0, p.Δm21, p.Δm31] atol = 1e-12
    q = effective_parameters(std, p; E = 2.5, ρ = 2.848)
    @test sin(q.θ12)^2 > 0.9 && sin(q.θ13) > sin(p.θ13)        # above the solar resonance, θ13 enhanced
    qbar = effective_parameters(std, p; E = 2.5, ρ = 2.848, antineutrino = true)
    @test sin(qbar.θ13) < sin(p.θ13)

    # The SNSI example and the C engine
    snsi = load_model(joinpath(MODELS, "SNSI"))
    @test snsi.name == "SNSI" && isempty(snsi.warnings)
    couplings = Dict(
        1 => (eta_etau = 0.1, Gamma0 = 1e-14),
        2 => (eta_emu = 0.05, delta_emu = 1.0),
        3 => (eta_etau = -0.05, delta_etau = π / 3, Gamma0 = 5e-15),     # Γ in eV
    )
    use_model!(snsi)
    for (config, cp, E, pμe, pμμ) in SNSI_REFERENCE
        b = snsi_base(Int(config))
        q = OscParams(snsi, b; couplings[Int(config)]...)
        @test model_probability(snsi, q, :mu, :e; E, L = 1300, ρ = 2.848, antineutrino = cp < 0) ≈
              pμe rtol = 1e-10
        @test probability(q, :μ, :μ; E = [E], L = 1300, ρ = 2.848, antineutrino = cp < 0)[1] ≈ pμμ rtol =
            1e-10
    end
    standard_engine!()

    # A model string: a mixing matrix from redefined rotations, a named
    # hierarchy, additional terms with energy dependence
    m = model_from_string(
        """
        R12 = Rotation[1, 2, th];
        R13 = IdentityMatrix[];
        R23 = IdentityMatrix[];
        MixingMatrix = R23 . R13 . R12;
        Hierarchy = Mine;
        MassHierarchy[[1]] = {Mine, {m1 -> 0.01, m2 -> Sqrt[m1^2 + dm21], m3 -> Sqrt[m1^2 + dm31]}};
        MassMatrix = U . Diag[{m1, m2, m3}] . Dagger[U];
        AdditionalTerms[[2]] = {LIV, c E Diag[{1, -1, 0}]};
        """;
        parameters = "ParameterDefinitions = {{th, {Value -> 0.5}}, {c, {Value -> 0}}};",
    )
    @test Set(b.name for b in m.blocks if b.kind === :AdditionalTerms) == Set([:Matter, :LIV])
    p = OscParams(m)
    H = hamiltonian(m, p; E = 1.0)
    R = [cos(0.5) sin(0.5) 0; -sin(0.5) cos(0.5) 0; 0 0 1]
    m2 = [0.01, sqrt(0.01^2 + p.Δm21), sqrt(0.01^2 + p.Δm31)] .^ 2
    @test H ≈ R * Diagonal(m2) * R' / 2e9 rtol = 1e-12
    H2 = hamiltonian(m, OscParams(m; c = 1e-22); E = 2.0)
    @test H2 - hamiltonian(m, p; E = 2.0) ≈ 1e-22 * 2e9 * Diagonal([1, -1, 0]) rtol = 1e-10

    # Four flavours
    st = load_model(joinpath(MODELS, "Sterile"))
    @test GLoBES.nflavours(st) == 4
    p = OscParams(st)
    @test sum(model_probability(st, p, :mu, j; E = 2.0, L = 1300, ρ = 2.8) for j = 1:4) ≈ 1 rtol =
        1e-12
    @test model_probability(st, p, :mu, :s; E = 2.0, L = 1300) > 0
    p0 = OscParams(st; theta14 = 0, theta24 = 0)
    use_model!(st)
    a = probability(p0, :μ, :e; E, L = 1300, ρ = 2.8)
    standard_engine!()
    @test a ≈ probability(nufit_params("6.1"), :μ, :e; E, L = 1300, ρ = 2.8) rtol = 1e-10

    # Errors
    @test_throws ArgumentError model_from_string("Mixing[[1]] = {X, 1};")
    @test_throws ArgumentError model_from_string("AdditionalTerms[[1]] = A;")
    @test_throws ArgumentError model_from_string("MassHierarchy[[1]] = {NH, {m1 = 0}};")
    # A later entry of the same name replaces an earlier one; Off removes the default
    m = model_from_string(
        "AdditionalTerms[[1]] = {Matter, 2 A};\nAdditionalTerms[[1]] = {Matter, Off};",
    )
    @test hamiltonian(m, OscParams(m); E = 1.0, ρ = 3.0) ≈ hamiltonian(m, OscParams(m); E = 1.0)
    @test_throws ArgumentError OscParams(std; nonsense = 1.0)
    @test_throws DimensionMismatch hamiltonian(
        model_from_string("MixingMatrix = Diag[{1, 1}];"),
        OscParams(std);
        E = 1.0,
    )
    w = model_from_string("AdditionalTerms[[1]] = {X, g A Diag[{1, 0, 0}]};")
    @test any(occursin("g", x) for x in w.warnings)
end
