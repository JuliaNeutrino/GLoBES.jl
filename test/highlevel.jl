# Tests of the high-level API. They check it against the glb* functions,
# which are tested against the C library in runtests.jl.

const TRUTH = OscParams(θ12=asin(sqrt(0.3)), θ13=0.5asin(sqrt(0.1)), θ23=π/4, δ=3π/2,
                        Δm21=7.9e-5, Δm31=2.6e-3)

to_glb_params(p::OscParams, n_exp) =
    (g = glbAllocParams(); glbDefineParams(g, p.values[1:6]...); glbSetDensityParams(g, p.ρ, GLB_ALL); g)

@testset "high-level API" begin
    @testset "OscParams" begin
        p = OscParams(theta12=0.6, theta13=0.15, theta23=0.8, dm21=7e-5, dm31=2.5e-3)
        @test p.θ12 == 0.6 && p.δ == 0.0 && p.ρ == 1.0
        @test p.sin²2θ13 ≈ sin(0.3)^2
        @test p.sin2theta23 ≈ sin(0.8)^2
        q = OscParams(p; δ=1.0, ρ=[1.0, 1.1])
        @test q.δ == 1.0 && q.ρ == [1.0, 1.1] && p.δ == 0.0
        @test with_param(p, :sin²2θ13, 0.08).sin²2θ13 ≈ 0.08
        @test with_param(p, :sin2theta12, 0.31).sin²θ12 ≈ 0.31
        @test_throws ArgumentError OscParams(θ12=0.6)
        @test_throws ArgumentError OscParams(p; θ14=0.1)
        @test OscParams(p; δ=0.0) == p
        @test occursin("θ13", sprint(show, p))
    end

    s = Setup()
    with_setup(s) do
        nova = Experiment("NOvA.glb")
        @test experiments() == [nova]
        @test rule_names(nova)[1] == "Nu_E_Appearance"
        @test_throws GLoBESError chi2(TRUTH)               # nothing simulated yet
        simulate!(TRUTH)
        @test simulated_params() == TRUTH

        @testset "rates" begin
            r = event_rates(nova; rule="Nu_E_Appearance")
            @test r.total == glbGetRuleRatePtr(1, 1)
            @test sum(r.total) ≈ glbTotalRuleRate(1, 1, GLB_ALL, GLB_W_EFF, GLB_W_BG, GLB_W_COEFF, GLB_SIG) +
                                 glbTotalRuleRate(1, 1, GLB_ALL, GLB_W_EFF, GLB_W_BG, GLB_W_COEFF, GLB_BG)
            @test event_rates(nova, TRUTH; rule=2).total ≈ event_rates(nova; rule=2).total
            other = OscParams(TRUTH; Δm31=2.4e-3)
            @test total_events(nova, other; rule=2) != total_events(nova; rule=2)
            @test event_rates(nova; rule=2).total == glbGetRuleRatePtr(1, 2)  # truth unchanged
            @test channel_rates(nova, 4) == glbGetChannelRates(1, 4, GLB_POST)
            @test bin_centers(nova) == glbGetBinCentersListPtr(1)
        end

        @testset "χ² against glb*" begin
            test = OscParams(TRUTH; δ=0.0, θ23=0.7)
            g = to_glb_params(test, 1)
            @test chi2(test) == glbChiSys(g, GLB_ALL, GLB_ALL)
            @test chi2(test; experiment=nova, rule=2) == glbChiSys(g, 1, 2)
            @test chi2(TRUTH) ≈ 0 atol = 1e-9
            nosys = chi2(test; systematics=false)
            @test nosys > chi2(test)
            @test glbGetSysOnOffState(1, 1) == GLB_ON          # restored
            @test chi2_scan(TRUTH, :δ => [0.0, π]) == [chi2(OscParams(TRUTH; δ=x)) for x in (0.0, π)]
            grid = chi2_scan(TRUTH, :sin²2θ13 => [0.09, 0.1], :δ => [0.0, 3π/2])
            @test size(grid) == (2, 2)
            @test grid[2, 2] ≈ 0 atol = 1e-9
        end

        @testset "marginalize" begin
            start = OscParams(TRUTH; δ=0.0)
            r = marginalize(start; fixed=(:θ12, :θ13, :δ, :Δm21), priors=(ρ=0.05,))
            @test Set(r.free) == Set([:θ23, :Δm31, :ρ])
            @test r.params.δ == 0.0 && r.params.θ12 == TRUTH.θ12
            @test r.chi2 < chi2(start)
            @test chi2(r.params) + ((r.params.ρ - 1) / 0.05)^2 ≈ r.chi2 rtol = 1e-6

            # the same fit with glb*
            proj = glbAllocProjection()
            glbDefineProjection(proj, GLB_FIXED, GLB_FIXED, GLB_FREE, GLB_FIXED, GLB_FIXED, GLB_FREE)
            glbSetDensityProjectionFlag(proj, GLB_FREE, GLB_ALL)
            glbSetProjection(proj)
            err = glbAllocParams(); glbDefineParams(err, 0, 0, 0, 0, 0, 0); glbSetDensityParams(err, 0.05, GLB_ALL)
            glbSetInputErrors(err); glbSetCentralValues(to_glb_params(TRUTH, 1))
            @test glbChiNP(to_glb_params(start, 1), nothing, GLB_ALL) ≈ r.chi2 rtol = 1e-12

            # marginalize leaves the glb* projection and priors as they were
            marginalize(start; free=(:θ23,))
            after = glbAllocProjection(); glbGetProjection(after)
            @test after.osc == [GLB_FIXED, GLB_FIXED, GLB_FREE, GLB_FIXED, GLB_FIXED, GLB_FREE]
            e2 = glbAllocParams(); glbGetInputErrors(e2)
            @test glbGetDensityParams(e2, 1) == 0.05

            # a user-defined prior adds to χ²
            plain = marginalize(start; free=(:θ23,))
            with_prior = marginalize(start; free=(:θ23,), prior=p -> 100 * (p.θ23 - π/4)^2)
            @test with_prior.chi2 > plain.chi2
            @test abs(with_prior.params.θ23 - π/4) < abs(plain.params.θ23 - π/4)
            @test_throws ArgumentError marginalize(start; free=(:θ23,), fixed=(:δ,))
            @test_throws ArgumentError marginalize(start; free=(:nope,))
        end

        @testset "experiment properties" begin
            L = baseline(nova)
            baseline!(nova, 1000.0)
            @test baseline(nova) == 1000.0 == glbGetBaselineInExperiment(1)
            baseline!(nova, L)
            m = target_mass(nova)
            target_mass!(nova, 2m)
            @test glbGetTargetMass(1) == 2m
            target_mass!(nova, m)
            @test energy_window_bins(nova) == 1:20
            lengths, ρs = prem_profile(3000.0, 10)
            @test sum(lengths) ≈ 3000.0 && length(ρs) == 10
            @test_throws GLoBESError chi2(TRUTH; experiment=nova, rule="nope")
        end

        @testset "probabilities" begin
            @test probability(TRUTH, :μ, :e; E=2.0, L=810, ρ=2.8) ≈
                  (glbSetOscillationParameters(to_glb_params(TRUTH, 1));
                   glbConstantDensityProbability(2, 1, +1, 2.0, 810.0, 2.8))
            @test probability(TRUTH, :mu, :e; E=2.0, L=810) == glbVacuumProbability(2, 1, +1, 2.0, 810.0)
            Ps = probability(TRUTH, :μ, :e, nova; E=[1.5, 2.0], antineutrino=true)
            @test Ps == [glbProfileProbability(1, 2, 1, -1, E) for E in (1.5, 2.0)]
            @test sum(probability(TRUTH, :μ, f; E=1.0, L=1300, ρ=3.0) for f in (:e, :μ, :τ)) ≈ 1
            @test_throws ArgumentError probability(TRUTH, :x, :e; E=1.0, L=1.0)
        end

        @testset "user-defined χ²" begin
            # Reproduce chiNoSysSpectrum
            define_chi2("myNoSys", 0) do e, rule, x, errors
                data = true_rates(e; rule)
                s, b = fit_rates(e; rule)
                bgc = e.data.rules[rule].bg_centers[1]
                sum(poisson_chi2(data[i], s[i] + bgc * b[i]) for i in energy_window_bins(e; rule))
            end
            test = OscParams(TRUTH; θ23=0.7)
            ref = chi2(test; systematics=false)
            chi2_function!(nova, "myNoSys")
            @test chi2_function(nova; rule=3) == "myNoSys"
            @test chi2(test) ≈ ref rtol = 1e-12
            chi2_function!(nova, "chiSpectrumTilt"; errors=[0.05, 0.0001, 0.05, 0.0001], rule=1)
            @test systematic_errors(nova; rule=1) == [0.05, 0.0001, 0.05, 0.0001]
            @test_throws GLoBESError chi2_function!(nova, "noSuchFunction")
        end

        @testset "exceptions" begin
            @test_throws GLoBESError Experiment("does-not-exist.glb")
            @test length(experiments()) == 1
            @test_throws GLoBESError OscParams(TRUTH; ρ=[1.0, 1.0, 1.0]) |> chi2
        end
    end

    @testset "setups are independent" begin
        a = Setup(); b = Setup()
        with_setup(() -> Experiment("NOvA.glb"; BASELINE=812), a)
        with_setup(() -> Experiment("T2K.glb"), b)
        simulate!(TRUTH; setup=a)
        simulate!(OscParams(TRUTH; δ=0.0); setup=b)
        test = OscParams(TRUTH; δ=π/2)
        ca = chi2(test; setup=a)
        cb = chi2(test; setup=b)
        @test ca != cb
        @test chi2(TRUTH; setup=a) ≈ 0 atol = 1e-9
        @test chi2(TRUTH; setup=b) > 0
        @test length(experiments(setup=a)) == 1 && name(experiments(setup=b)[1]) == "T2K"
        @test current_setup() !== a && current_setup() !== b
        clear_experiments!(setup=a)
        @test isempty(experiments(setup=a))
    end

    @testset "probability engine" begin
        s = Setup()
        with_setup(s) do
            nova = Experiment("NOvA.glb")
            simulate!(TRUTH)
            standard = chi2(OscParams(TRUTH; δ=0.0))
            Pstd = probability(TRUTH, :μ, :e; E=2.0, L=812, ρ=2.8)
            # The standard probabilities, damped by an extra parameter
            engine = GLoBES.StandardEngine()
            probability_engine!((:damping,)) do P, p, E, anti, lengths, densities
                GLoBES.set_oscillation_parameters!(engine, p.values[1:6])
                GLoBES.probability_matrix!(P, engine, anti ? -1 : +1, E, length(lengths),
                                           lengths, densities, -1.0)
                P .= (1 - p.damping) .* P
            end
            t = OscParams(TRUTH; extra=(damping=0.0,))
            @test_throws GLoBESError chi2(t)                  # simulate! again for the new engine
            simulate!(t)
            @test chi2(OscParams(t; δ=0.0)) ≈ standard rtol = 1e-12
            @test probability(t, :μ, :e; E=2.0, L=812, ρ=2.8) ≈ Pstd rtol = 1e-12
            @test probability(OscParams(t; damping=0.5), :μ, :e; E=2.0, L=812, ρ=2.8) ≈ 0.5Pstd rtol = 1e-12
            r = marginalize(OscParams(t; damping=0.1); free=(:damping,))
            @test abs(r.params.damping) < 1e-3
            @test_throws GLoBESError chi2(TRUTH)              # 6 parameters for a 7-parameter engine
            standard_engine!()
            simulate!(TRUTH)
            @test chi2(OscParams(TRUTH; δ=0.0)) ≈ standard rtol = 1e-12
        end
    end
end
