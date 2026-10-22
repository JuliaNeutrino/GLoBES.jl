using GLoBES
using Test
using LinearAlgebra
using SpecialFunctions: erf

# Reference numbers from the C library, see reference/make_reference.c.
# The setups below must match the ones there.
include("reference_values.jl")

const DATADIR = joinpath(@__DIR__, "data")

const SETUPS = [
    ("NOvA", ["NOvA.glb"], Pair{String, Float64}[]),
    ("NOvA_prepost", ["NOvA_pre_post_smearing.glb"], Pair{String, Float64}[]),
    ("T2K", ["T2K.glb"], Pair{String, Float64}[]),
    ("T2HK", ["T2HK.glb"], Pair{String, Float64}[]),
    ("NFstandard", ["NFstandard.glb"], Pair{String, Float64}[]),
    ("NFvar", ["NFvar.glb"], ["BASELINE" => 3000.0, "emax" => 50.0]),
    ("NF_hR_lT", ["NF_hR_lT.glb"], ["BASELINE" => 4000.0, "emax" => 50.0]),
    ("BB_350", ["BB_350.glb"], Pair{String, Float64}[]),
    (
        "BBvar_WC",
        ["BBvar_WC.glb"],
        ["gammafactor" => 350.0, "baselinefactor" => 730.0 / 350.0, "EXP_FACTOR" => 1.0],
    ),
    ("Reactor1", ["Reactor1.glb"], Pair{String, Float64}[]),
    ("Reactor2", ["Reactor2.glb"], Pair{String, Float64}[]),
    ("DChooz", ["D-Chooz_far.glb", "D-Chooz_near.glb"], Pair{String, Float64}[]),
]

# chiSpectrumCalib is not compared with C: GLoBES.jl uses a different energy
# scale shift (see "energy scale shift" below), and setups that use it are
# switched to chiSpectrumTilt here and in make_reference.c.
const ALT_CHI = [
    "chiSpectrumTilt" => [0.025, 0.01, 0.05, 0.01],
    "chiNoSysSpectrum" => Float64[],
    "chiSpectrumOnly" => [0.0],
    "chiTotalRatesTilt" => [0.025, 0.01, 0.05, 0.01],
    "chiNoSysTotalRates" => Float64[],
    "chiZero" => Float64[],
]

# Tolerances. Rates and χ² at fixed nuisance parameters must agree to rounding.
# glbChiSys still minimizes over the systematics (Powell, TOLSYS = 1e-8), and
# the glbChi* projections also over oscillation parameters (TOLOSC = 1e-5).
# Fitted parameter values are only loosely determined along flat directions.
const RTOL_EXACT = 1e-9
const RTOL_SYS = 1e-6
const ATOL_SYS = 1e-8
const RTOL_FIT = 1e-3
const ATOL_FIT = 1e-6
const RTOL_FIT_PARAMS = 2e-2

# DChooz's own χ² function (chiDCNorm, from example5.c) shifts the energy
# scale with glbShiftEnergyScale, which GLoBES.jl does differently from C, so
# its χ² values agree with C only to about a percent.
const RTOL_DCHOOZ = 3e-2

# C GLoBES returns χ² = 1e100 whenever a fit rate is ≤ 0, even for bins with
# no true and no fit events (NF_hR_lT and BBvar_WC have such bins above the
# beam energy). GLoBES.jl fixes this, so where C gives the sentinel we only
# require a finite χ².
const C_SENTINEL = 1e99
matches_c(x, key; rtol, atol = 0.0) =
    R[key] >= C_SENTINEL ? (isfinite(x) && -atol <= x < C_SENTINEL) :
    isapprox(x, R[key]; rtol, atol)

function define_params(th12, th13, th23, delta, dm21, dm31, rho)
    p = glbAllocParams()
    glbDefineParams(p, th12, th13, th23, delta, dm21, dm31)
    glbSetDensityParams(p, rho, GLB_ALL)
    return p
end
true_point() = define_params(asin(sqrt(0.307)), 0.149, asin(sqrt(0.55)), 1.2π, 7.4e-5, 2.5e-3, 1.0)
test_point() = define_params(asin(sqrt(0.307)), 0.149, asin(sqrt(0.45)), 0.5π, 7.4e-5, 2.45e-3, 1.0)

params_vector(p) = [
    [glbGetOscParams(p, i) for i = 1:6];
    [glbGetDensityParams(p, i) for i = 1:glbNumOfExps()]
]

# example5.c of the C distribution: Double Chooz far/near with flux, fiducial
# mass and energy calibration errors
const EXP_FAR = 1
const EXP_NEAR = 2
function chiDCNorm(exp, rule, n_params, x, errors, user_data)
    n_bins = glbGetNumberOfBins(EXP_FAR)
    true_N = glbGetRuleRatePtr(EXP_NEAR, 1)
    true_F = glbGetRuleRatePtr(EXP_FAR, 1)
    emin, emax = glbGetEminEmax(exp)
    ew_low, ew_high = glbGetEnergyWindowBins(exp, rule)
    fit_F = glbShiftEnergyScale(x[4], glbGetSignalFitRatePtr(EXP_FAR, 1), n_bins, emin, emax)
    fit_N = glbShiftEnergyScale(x[5], glbGetSignalFitRatePtr(EXP_NEAR, 1), n_bins, emin, emax)
    nF = 1.0 + x[1] + x[2]
    nN = 1.0 + x[1] + x[3]
    chi2 = 0.0
    for i = ew_low:ew_high
        true_F[i] > 0 && (chi2 += (true_F[i] - nF * fit_F[i])^2 / true_F[i])
        true_N[i] > 0 && (chi2 += (true_N[i] - nN * fit_N[i])^2 / true_N[i])
    end
    for i = 1:n_params
        chi2 += (x[i] / errors[i])^2
    end
    return chi2
end

function load_setup(files, vars)
    glbClearAEDLVariables()
    for (name, value) in vars
        glbDefineAEDLVariable(name, value)
    end
    glbClearExperimentList()
    for f in files
        glbInitExperiment(f)
    end
    replace_calib()
end

function replace_calib()
    for e = 1:glbNumOfExps(), r = 1:glbGetNumberOfRules(e)
        if glbGetChiFunction(e, r, GLB_ON) == "chiSpectrumCalib"
            glbSetChiFunction(
                e,
                r,
                GLB_ON,
                "chiSpectrumTilt",
                copy(glbGetSysErrorsListPtr(e, r, GLB_ON)),
            )
        end
    end
end

function run_setup(name, files, vars)
    rtol_sys = name == "DChooz" ? RTOL_DCHOOZ : RTOL_SYS
    rtol_fit = name == "DChooz" ? RTOL_DCHOOZ : RTOL_FIT
    load_setup(files, vars)
    @test glbNumOfExps() == length(files)

    p0, pt, fit = true_point(), test_point(), glbAllocParams()
    err = define_params(
        0.1 * asin(sqrt(0.307)),
        0,
        0.1 * asin(sqrt(0.55)),
        0,
        0.05 * 7.4e-5,
        0.05 * 2.5e-3,
        0.05,
    )
    glbSetOscillationParameters(p0)
    glbSetRates()
    glbSetCentralValues(p0)
    glbSetInputErrors(err)
    glbSetStartingValues(p0)

    @testset "rates" begin
        for e = 1:glbNumOfExps()
            @test glbGetBinCentersListPtr(e) ≈ R["$name/$e/bincenters"] rtol = RTOL_EXACT
            for r = 1:glbGetNumberOfRules(e)
                k = "$name/$e/rule$r"
                # C leaves the rule rate outside the energy window uninitialised
                lo, hi = glbGetEnergyWindowBins(e, r)
                @test glbGetRuleRatePtr(e, r)[lo:hi] ≈ R["$k/rate"][lo:hi] rtol = RTOL_EXACT
                @test glbGetSignalRatePtr(e, r) ≈ R["$k/signal"] rtol = RTOL_EXACT
                @test glbGetBGRatePtr(e, r) ≈ R["$k/bg"] rtol = RTOL_EXACT
                @test glbTotalRuleRate(e, r, GLB_ALL, GLB_W_EFF, GLB_W_BG, GLB_W_COEFF, GLB_SIG) ≈
                      R["$k/total"] rtol = RTOL_EXACT
            end
            for c = 1:glbGetNumberOfChannels(e)
                @test glbGetChannelRatePtr(e, c, GLB_POST) ≈ R["$name/$e/channel$c/post"] rtol =
                    RTOL_EXACT
            end
        end
    end

    @testset "chi2 without fit" begin
        @test matches_c(
            glbChiSys(pt, GLB_ALL, GLB_ALL),
            "$name/chisys";
            rtol = rtol_sys,
            atol = ATOL_SYS,
        )
        for e = 1:glbNumOfExps(), r = 1:glbGetNumberOfRules(e)
            @test matches_c(
                glbChiSys(pt, e, r),
                "$name/$e/rule$r/chisys";
                rtol = rtol_sys,
                atol = ATOL_SYS,
            )
        end
        glbSwitchSystematics(GLB_ALL, GLB_ALL, GLB_OFF)
        @test matches_c(
            glbChiSys(pt, GLB_ALL, GLB_ALL),
            "$name/chisys_off";
            rtol = rtol_sys,
            atol = ATOL_SYS,
        )
        glbSwitchSystematics(GLB_ALL, GLB_ALL, GLB_ON)
    end

    @testset "chi2 with fit" begin
        @test matches_c(
            glbChiTheta13(pt, fit, GLB_ALL),
            "$name/chitheta13";
            rtol = rtol_fit,
            atol = ATOL_FIT,
        )
        # where C's χ² is the 1e100 plateau, its fitted parameters are arbitrary
        chidelta = glbChiDelta(pt, fit, GLB_ALL)
        @test matches_c(chidelta, "$name/chidelta"; rtol = rtol_fit, atol = ATOL_FIT)
        R["$name/chidelta"] < C_SENTINEL &&
            @test params_vector(fit) ≈ R["$name/chidelta_fit"] rtol = RTOL_FIT_PARAMS
        @test matches_c(
            glbChiTheta13Delta(pt, nothing, GLB_ALL),
            "$name/chitheta13delta";
            rtol = rtol_fit,
            atol = ATOL_FIT,
        )

        proj = glbAllocProjection()
        glbDefineProjection(proj, GLB_FIXED, GLB_FIXED, GLB_FREE, GLB_FIXED, GLB_FIXED, GLB_FREE)
        glbSetDensityProjectionFlag(proj, GLB_FREE, GLB_ALL)
        glbSetProjection(proj)
        @test matches_c(
            glbChiNP(pt, nothing, GLB_ALL),
            "$name/chinp";
            rtol = rtol_fit,
            atol = ATOL_FIT,
        )
    end

    if name != "DChooz"
        @testset "built-in chi2 functions" begin
            for (chi, errs) in ALT_CHI
                glbSetChiFunction(GLB_ALL, GLB_ALL, GLB_ON, chi, errs)
                @test matches_c(
                    glbChiSys(pt, GLB_ALL, GLB_ALL),
                    "$name/alt/$chi";
                    rtol = RTOL_SYS,
                    atol = ATOL_SYS,
                )
            end
            # At fixed nuisance parameters, with the fit rates of the last glbChiSys
            x = [0.01, -0.005, 0.02, 0.003]
            for e = 1:glbNumOfExps(), r = 1:glbGetNumberOfRules(e), (chi, errs) in ALT_CHI
                f = glbGetChiFunctionPtr(chi)
                @test matches_c(
                    f(e, r, length(errs), x, errs, nothing),
                    "$name/$e/rule$r/fixed/$chi";
                    rtol = RTOL_EXACT,
                    atol = ATOL_SYS,
                )
            end
        end
    end

    @testset "exposure" begin
        glbSetTargetMass(1, 2.0 * glbGetTargetMass(1))
        glbSetRunningTime(1, 1, 0.5 * glbGetRunningTime(1, 1))
        glbSetRates()
        for r = 1:glbGetNumberOfRules(1)
            @test glbTotalRuleRate(1, r, GLB_ALL, GLB_W_EFF, GLB_W_BG, GLB_W_COEFF, GLB_SIG) ≈
                  R["$name/1/rule$r/total_scaled"] rtol = RTOL_EXACT
        end
    end
end

cd(DATADIR) do
    glbSetVerbosityLevel(0)
    glbDefineChiFunction(chiDCNorm, 5, "chiDCNorm")

    @testset "GLoBES.jl" begin
        @testset "probabilities" begin
            energies = R["prob/energies"]
            load_setup(["NFstandard.glb"], [])
            glbSetOscillationParameters(true_point())
            probs(f) = [f(a, b, cp, E) for cp in (+1, -1) for E in energies for a = 1:3 for b = 1:3]
            @test probs((a, b, cp, E) -> glbVacuumProbability(a, b, cp, E, 1300.0)) ≈
                  R["prob/vacuum"] rtol = RTOL_EXACT atol = 1e-14
            @test probs((a, b, cp, E) -> glbConstantDensityProbability(a, b, cp, E, 1300.0, 2.8)) ≈
                  R["prob/constant"] rtol = RTOL_EXACT atol = 1e-14
            @test probs((a, b, cp, E) -> glbProfileProbability(1, a, b, cp, E)) ≈ R["prob/profile"] rtol =
                RTOL_EXACT atol = 1e-14

            # Unitarity
            for E in energies, cp in (+1, -1), a = 1:3
                @test sum(glbConstantDensityProbability(a, b, cp, E, 1300.0, 2.8) for b = 1:3) ≈ 1 atol =
                    1e-12
            end
        end

        @testset "energy scale shift" begin
            # Spectra with known bin integrals: G is the cumulative number of events
            # (zero outside the binned range), so a shift by b moves bin i to
            # G((1+b) hi) - G((1+b) lo)
            binned(G, edges) = [G(edges[i + 1]) - G(edges[i]) for i = 1:(length(edges) - 1)]
            centers(edges) = (edges[1:(end - 1)] .+ edges[2:end]) ./ 2
            equal = collect(range(0.5, 10.5, length = 21))
            unequal = [
                0.5,
                1.5,
                2.5,
                3.0,
                3.5,
                4.0,
                4.25,
                4.5,
                4.75,
                5.0,
                5.25,
                5.5,
                6.0,
                6.5,
                7.5,
                8.5,
                10.5,
            ]
            spectra = [
                ("peak", E -> 0.5 * erf((E - 5) / 1.5)),
                ("rising to the upper end", E -> E^3 / 300),
                ("falling from the lower end", E -> -exp(-E / 2)),
            ]
            for edges in (equal, unequal), (label, cdf) in spectra
                G = E -> cdf(clamp(E, edges[1], edges[end]))
                rates = binned(G, edges)
                shift = b -> GLoBES.shift_energy_scale(b, rates, centers(edges), diff(edges))
                @test shift(0.0) ≈ rates rtol = 1e-14
                for b in (-0.05, -0.01, 0.003, 0.01, 0.05)
                    exact = binned(E -> G((1 + b) * E), edges)
                    out = shift(b)
                    @test all(out .>= 0)
                    # events are conserved, up to the resolution of the last bin
                    @test sum(out) ≈ sum(exact) rtol = (edges === equal ? 1e-3 : 5e-3)
                    # C GLoBES' linear interpolation is off by up to 44% here
                    @test maximum(abs.(out .- exact)) <=
                          (edges === equal ? 2e-3 : 5e-2) * maximum(exact)
                end
            end

            # Smooth in b at b = 0 for a spectrum that vanishes at both ends
            # (C's version has a kink there in every bin)
            peak = binned(E -> 0.5 * erf((clamp(E, 0.5, 10.5) - 5) / 1.0), equal)
            shifted = b -> glbShiftEnergyScale(b, peak, 20, 0.5, 10.5)
            h = 1e-6
            @test (shifted(h) .- shifted(0.0)) ./ h ≈ (shifted(0.0) .- shifted(-h)) ./ h rtol = 1e-3 atol =
                1e-6

            # Shifting past the upper end must not create events (C doubles the last bin)
            flat = ones(20)
            @test glbShiftEnergyScale(0.01, flat, 20, 0.5, 10.5)[end] < 1
        end

        @testset "Poisson likelihood" begin
            L = glb_likelihood
            @test L(0.0, 0.0) == 0                       # empty bin (C: 2e100)
            @test L(0.0, 1.5) ≈ 3.0
            @test L(4.0, 4.0) == 0
            @test L(4.0, 2.0) ≈ 2 * (2.0 - 4.0 + 4.0 * log(2.0))
            # with observed events, a vanishing or negative fit rate gives a large,
            # finite, continuous and increasing penalty instead of a jump to 1e100
            f0 = GLoBES.FIT_RATE_FLOOR
            @test L(4.0, f0 * (1 + 1e-9)) ≈ L(4.0, f0 * (1 - 1e-9)) rtol = 1e-9
            @test L(4.0, 1e-3) < L(4.0, 0.0) < L(4.0, -1.0) < Inf
        end

        @testset "chiSpectrumCalib" begin
            load_setup(["NOvA.glb"], [])
            glbSetChiFunction(
                GLB_ALL,
                GLB_ALL,
                GLB_ON,
                "chiSpectrumCalib",
                [0.025, 0.01, 0.05, 0.01],
            )
            glbSetOscillationParameters(true_point())
            glbSetRates()
            pt = test_point()
            glbChiSys(pt, 1, 1)                          # computes the fit rates at pt
            # without energy scale shifts, the χ² is the plain normalisation fit
            x = [0.01, 0.0, 0.02, 0.0]
            errs = [0.025, 0.01, 0.05, 0.01]
            true_rates = glbGetRuleRatePtr(1, 1)
            sig, bg = glbGetSignalFitRatePtr(1, 1), glbGetBGFitRatePtr(1, 1)
            bg_center = glbGetBGCenters(1, 1)[1]
            lo, hi = glbGetEnergyWindowBins(1, 1)
            expected =
                sum(
                    glb_likelihood(true_rates[i], 1.01 * sig[i] + 1.02 * bg_center * bg[i]) for
                    i = lo:hi
                ) +
                (0.01 / 0.025)^2 +
                (0.02 * bg_center / 0.05)^2
            @test glbGetChiFunctionPtr("chiSpectrumCalib")(1, 1, 4, x, errs, nothing) ≈ expected rtol =
                1e-12
            # its minimum no longer depends on the last digits of the starting point
            chi(eps) = (
                p = test_point();
                glbSetOscParams(p, glbGetOscParams(p, GLB_THETA_23) * (1 + eps), GLB_THETA_23);
                glbChiSys(p, 1, 1)
            )
            @test chi(1e-12) ≈ chi(0.0) rtol = 1e-9
            @test chi(-1e-12) ≈ chi(0.0) rtol = 1e-9
        end

        @testset "minimizer restarts" begin
            # The sgn(Δm²31) degeneracy fit of example-tour.c, a long curved
            # valley where a single Powell run can stop early. C GLoBES (with
            # test/reference/c-fixes.patch) reaches 4.555672592; the restarted
            # minimizer must get there too, within its tolerance TOLOSC.
            load_setup(["NFstandard.glb"], [])
            th12, sdm, ldm = asin(sqrt(0.8)) / 2, 7e-5, 2e-3
            truth = define_params(th12, asin(sqrt(0.001)) / 2, π / 4, π / 2, sdm, ldm, 1.0)
            glbSetOscillationParameters(truth)
            glbSetRates()
            start = define_params(th12, asin(sqrt(0.001)) / 2, π / 4, π / 2, sdm, -ldm, 1.0)
            glbSetCentralValues(start)
            glbSetInputErrors(define_params(0.1 * th12, 0, 0, 0, 0.1 * sdm, ldm / 3, 0.05))
            @test glbChiAll(start, nothing, GLB_ALL) <= 4.555672592 * (1 + 1e-5)
        end

        for (name, files, vars) in SETUPS
            @testset "$name" run_setup(name, files, vars)
        end

        @testset "user-defined probability engine" begin
            # Wrapping the default engine must give identical rates
            load_setup(["NOvA.glb"], [])
            p0 = true_point()
            glbSetOscillationParameters(p0)
            glbSetRates()
            reference = copy(glbGetRuleRatePtr(1, 1))
            calls = Ref(0)
            prob =
                (P, cp, E, psteps, L, rho, fs, ud) ->
                    (calls[] += 1; GLoBES.default_prob_func(P, cp, E, psteps, L, rho, fs, ud))
            glbRegisterProbabilityEngine(6, prob, nothing, nothing)
            try
                glbSetOscillationParameters(p0)
                glbSetRates()
                @test calls[] > 0
                @test glbGetRuleRatePtr(1, 1) ≈ reference rtol = 1e-14
            finally
                glbRegisterProbabilityEngine(6, nothing, nothing, nothing)
            end
        end

        @testset "user-defined prior" begin
            load_setup(["NOvA.glb"], [])
            p0, pt = true_point(), test_point()
            glbSetOscillationParameters(p0)
            glbSetRates()
            base = glbChiSys(pt, GLB_ALL, GLB_ALL)
            glbRegisterPriorFunction((p, ud) -> 42.0, nothing, nothing)
            try
                # glbChiSys does not fit, but glbChiDelta adds the prior at the minimum
                @test glbChiDelta(pt, nothing, GLB_ALL) ≥ 42.0
            finally
                glbRegisterPriorFunction(nothing, nothing, nothing)
            end
            @test glbChiSys(pt, GLB_ALL, GLB_ALL) == base
        end

        @testset "parameter vectors" begin
            glbClearExperimentList()
            glbInitExperiment("NOvA.glb")
            p = true_point()
            q = glbAllocParams()
            glbCopyParams(p, q)
            @test params_vector(q) == params_vector(p)
            glbSetOscParams(q, 0.3, GLB_DELTA_CP)
            @test glbGetOscParams(q, GLB_DELTA_CP) == 0.3
            @test glbGetOscParams(p, GLB_DELTA_CP) == 1.2π
            glbSetOscillationParameters(p)
            r = glbAllocParams()
            glbGetOscillationParameters(r)
            @test params_vector(r)[1:6] ≈ params_vector(p)[1:6]
        end

        @testset "errors" begin
            glbClearExperimentList()
            @test glbInitExperiment("does-not-exist.glb") != 0
            @test glbNumOfExps() == 0
            @test glbGetNumberOfBins(1) == -1
        end

        include("highlevel.jl")
    end
end

# Beyond GLoBES: models, NuFIT, new physics, LaTeX and analyses
include("mathematica.jl")
include("sarah.jl")
include("mechanisms.jl")
include("nufit.jl")
include("newphysics.jl")
include("oscmodel.jl")
include("latex.jl")
include("analysis.jl")
