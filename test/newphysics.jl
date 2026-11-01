# The new-physics engine through the PropagationBlock interface, with blocks
# defined here (scalar NSI, vector NSI and decoherence are packages of their
# own, tested there against the C engine of SNSI-Decoherence-DUNE, as is the
# SNSI model file in oscmodel.jl)

# Matter term Δε_ee: the standard potential times 1 + Δε
struct TestMatterShift <: GLoBES.PropagationBlock end
GLoBES.block_parameters(::TestMatterShift) = [:Δε_test]
GLoBES.parameter_aliases(::TestMatterShift) = [:deps_test => :Δε_test]
GLoBES.matter_correction!(ε, ::TestMatterShift, p, x) = (ε[1, 1] += x[:Δε_test]; ε)

# The same Γ for all pairs, times E/GeV
struct TestDamping <: GLoBES.PropagationBlock end
GLoBES.block_parameters(::TestDamping) = [:Γ_test]
GLoBES.damps(::TestDamping) = true
GLoBES.damping_rates(::TestDamping, p, x, E) = ntuple(_ -> E * x[:Γ_test], 3)

# A mass-matrix correction proportional to the identity, δM = c 𝟙
struct TestMassShift <: GLoBES.PropagationBlock end
GLoBES.block_parameters(::TestMassShift) = [:c_test]
GLoBES.uses_absolute_masses(::TestMassShift) = true
GLoBES.mass_correction!(δM, ::TestMassShift, p, x) = (δM[diagind(δM)] .+= x[:c_test]; δM)

struct TestIncomplete <: GLoBES.PropagationBlock end
struct TestNoRates <: GLoBES.PropagationBlock end
GLoBES.block_parameters(::TestNoRates) = [:Γ_none]
GLoBES.damps(::TestNoRates) = true

@testset "new physics engine" begin
    p = nufit_params("6.1")
    E = [0.5, 1.0, 2.5, 5.0]
    dune = (; L = 1300, ρ = 2.848)
    standard = [probability(p, :μ, :e; E, dune..., antineutrino = a) for a in (false, true)]

    model = NewPhysics(TestMassShift(), TestMatterShift(), TestDamping())
    @test model.names == [:m_lightest, :c_test, :Δε_test, :Γ_test]
    use_model!(model)

    # No new physics: the standard engine, whatever the lightest mass
    for mℓ in (0.0, 0.03), anti in (false, true)
        q = OscParams(model, p; m_lightest = mℓ)
        @test probability(q, :μ, :e; E, dune..., antineutrino = anti) ≈ standard[1 + anti] rtol =
            1e-10
    end

    # δM = c 𝟙 shifts all masses by c: with m_lightest = 0.01 and c = 0.02,
    # the masses of m_lightest = 0.03 (normal ordering)
    q = OscParams(model, p; m_lightest = 0.01, c_test = 0.02)
    shifted = GLoBES.light_masses(p; lightest = 0.01) .+ 0.02
    p_shifted = OscParams(p; Δm21 = shifted[2]^2 - shifted[1]^2, Δm31 = shifted[3]^2 - shifted[1]^2)
    P = probability(q, :μ, :e; E, dune...)
    standard_engine!()
    @test P ≈ probability(p_shifted, :μ, :e; E, dune...) rtol = 1e-10

    # Δε_ee = x is the standard potential times 1 + x (ASCII alias accepted)
    use_model!(model)
    P = probability(OscParams(model, p; deps_test = 0.5), :μ, :e; E, dune..., antineutrino = true)
    standard_engine!()
    @test P ≈ probability(p, :μ, :e; E, L = 1300, ρ = 1.5 * 2.848, antineutrino = true) rtol = 1e-10

    # Damping: large Γ gives the average of the oscillations, and probability is conserved
    use_model!(model)
    U = pmns_matrix(p)
    average = sum(abs2(U[2, i]) * abs2(U[1, i]) for i = 1:3)
    d = OscParams(model, p; Γ_test = 1e-18)
    @test probability(d, :μ, :e; E = [1.0], L = 1300, ρ = 0.0)[1] ≈ average rtol = 1e-6
    d = OscParams(model, p; Γ_test = 1e-23, Δε_test = 0.2, c_test = 0.01)
    @test sum(probability(d, :μ, β; E, dune...) for β in (:e, :μ, :τ)) ≈ ones(4) rtol = 1e-12
    standard_engine!()

    # Interface errors
    @test_throws ArgumentError NewPhysics(TestMatterShift(), TestMatterShift())
    @test_throws ArgumentError NewPhysics(TestIncomplete())
    use_model!(NewPhysics(TestNoRates()))
    @test_throws ArgumentError probability(
        OscParams(NewPhysics(TestNoRates()), p),
        :μ,
        :e;
        E,
        dune...,
    )
    standard_engine!()

    @test GLoBES.flavour_pairs((:eτ, :mumu)) == [(1, 3), (2, 2)]
    @test GLoBES.coupling_parameters([(1, 1), (1, 3)], :ε, :δε) == [:ε_ee, :ε_eτ, :δε_eτ]
    @test GLoBES.coupling_aliases([(1, 3)], :ε, :δε, :eps, :deltaeps) ==
          [:eps_etau => :ε_eτ, :deltaeps_etau => :δε_eτ]
end
