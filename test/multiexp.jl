# GLoBES 3.2 features with the DUNE configuration: correlated systematics
# (sys groups and chiMultiExp), tables of any length and emin = 0. The
# numbers come from C GLoBES 3.2.18, see reference/dune_reference.c.

# The data tables of the DUNE configuration are downloaded from arXiv by
# dune_configuration(); without network access the test is skipped.
dune_file = try
    dune_configuration()
catch e
    @info "Could not download the DUNE data tables: skipping the DUNE tests" exception = e
    nothing
end
if dune_file !== nothing
    @testset "DUNE (GLoBES 3.2 systematics)" begin
        s = Setup()
        with_setup(s) do
            dune = Experiment(dune_file)
            @test rule_names(dune) == ["nue_app", "nuebar_app", "numu_dis", "numubar_dis"]
            @test dune.data.sys_errors["#err_nutau_bg"] == 0.2
            @test dune.data.rules[1].multiex_bg[5] == ["#err_nutau_bg"]
            truth = OscParams(
                θ12 = asin(sqrt(0.307)),
                θ13 = asin(sqrt(0.02203)),
                θ23 = asin(sqrt(0.572)),
                δ = deg2rad(198),
                Δm21 = 7.41e-5,
                Δm31 = 2.511e-3,
            )
            simulate!(truth)
            signal = [2330.191549, 511.372874, 20180.94105, 11810.09404]
            background = [898.4129854, 595.6927677, 451.2026604, 273.2343884]
            for r = 1:4
                @test sum(event_rates(dune; rule = r).signal) ≈ signal[r] rtol = 1e-6
                @test sum(event_rates(dune; rule = r).background) ≈ background[r] rtol = 1e-6
            end
            # Δχ² of δ = 0 and 90°: all rules share the systematics, so their sum is
            # not the sum of the rules minimised one by one
            for (δ, all, rules, nosys) in (
                (0, 105.7785831, [74.05701347, 17.20134517, 12.14434871, 2.066132805], 111.7673843),
                (
                    90,
                    150.9455862,
                    [124.7575348, 21.26299381, 2.865233035, 0.4885778645],
                    205.3334198,
                ),
            )
                test = OscParams(truth; δ = deg2rad(δ))
                @test chi2(test) ≈ all rtol = 1e-5
                @test [chi2(test; rule = r) for r = 1:4] ≈ rules rtol = 1e-5
                @test chi2(test; systematics = false) ≈ nosys rtol = 1e-5
            end
            @test chi2(OscParams(truth; Δm31 = -2.45e-3)) ≈ 980.6670283 rtol = 1e-5
            minimizer!(:hybrid; setup = s)
            @test_throws GLoBESError chi2(truth)
        end
    end
end
