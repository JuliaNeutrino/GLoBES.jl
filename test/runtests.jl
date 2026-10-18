using GLoBES
using Test

# Reference numbers from the C library, see reference/make_reference.c.
# The setups below must match the ones there.
include("reference_values.jl")

const DATADIR = joinpath(@__DIR__, "data")

const SETUPS = [
    ("NOvA",         ["NOvA.glb"],                   Pair{String,Float64}[]),
    ("NOvA_prepost", ["NOvA_pre_post_smearing.glb"], Pair{String,Float64}[]),
    ("T2K",          ["T2K.glb"],                    Pair{String,Float64}[]),
    ("T2HK",         ["T2HK.glb"],                   Pair{String,Float64}[]),
    ("NFstandard",   ["NFstandard.glb"],             Pair{String,Float64}[]),
    ("NFvar",        ["NFvar.glb"],    ["BASELINE" => 3000.0, "emax" => 50.0]),
    ("NF_hR_lT",     ["NF_hR_lT.glb"], ["BASELINE" => 4000.0, "emax" => 50.0]),
    ("BB_350",       ["BB_350.glb"],                 Pair{String,Float64}[]),
    ("BBvar_WC",     ["BBvar_WC.glb"],
                     ["gammafactor" => 350.0, "baselinefactor" => 730.0 / 350.0, "EXP_FACTOR" => 1.0]),
    ("Reactor1",     ["Reactor1.glb"],               Pair{String,Float64}[]),
    ("Reactor2",     ["Reactor2.glb"],               Pair{String,Float64}[]),
    ("DChooz",       ["D-Chooz_far.glb", "D-Chooz_near.glb"], Pair{String,Float64}[]),
]

const ALT_CHI = ["chiSpectrumTilt" => [0.025, 0.01, 0.05, 0.01], "chiNoSysSpectrum" => Float64[],
                 "chiSpectrumOnly" => [0.0], "chiTotalRatesTilt" => [0.025, 0.01, 0.05, 0.01],
                 "chiNoSysTotalRates" => Float64[], "chiSpectrumCalib" => [0.025, 0.01, 0.05, 0.01],
                 "chiZero" => Float64[]]

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

# chiSpectrumCalib shifts the energy scale by linear interpolation between
# bins, so χ² is only piecewise smooth in the calibration parameter. Its
# minimum moves by up to a few percent when the starting point changes in the
# 12th digit, in C as much as here. The function itself is checked exactly in
# "built-in chi2 functions" via the "fixed" values.
const RTOL_CALIB = 5e-2

function define_params(th12, th13, th23, delta, dm21, dm31, rho)
    p = glbAllocParams()
    glbDefineParams(p, th12, th13, th23, delta, dm21, dm31)
    glbSetDensityParams(p, rho, GLB_ALL)
    return p
end
true_point() = define_params(asin(sqrt(0.307)), 0.149, asin(sqrt(0.55)), 1.2π, 7.4e-5, 2.5e-3, 1.0)
test_point() = define_params(asin(sqrt(0.307)), 0.149, asin(sqrt(0.45)), 0.5π, 7.4e-5, 2.45e-3, 1.0)

params_vector(p) = [[glbGetOscParams(p, i) for i in 1:6];
                    [glbGetDensityParams(p, i) for i in 1:glbNumOfExps()]]

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
    for i in ew_low:ew_high
        true_F[i] > 0 && (chi2 += (true_F[i] - nF * fit_F[i])^2 / true_F[i])
        true_N[i] > 0 && (chi2 += (true_N[i] - nN * fit_N[i])^2 / true_N[i])
    end
    for i in 1:n_params
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
end

function run_setup(name, files, vars)
    load_setup(files, vars)
    @test glbNumOfExps() == length(files)

    p0, pt, fit = true_point(), test_point(), glbAllocParams()
    err = define_params(0.1 * asin(sqrt(0.307)), 0, 0.1 * asin(sqrt(0.55)), 0, 0.05 * 7.4e-5,
                        0.05 * 2.5e-3, 0.05)
    glbSetOscillationParameters(p0)
    glbSetRates()
    glbSetCentralValues(p0)
    glbSetInputErrors(err)
    glbSetStartingValues(p0)

    @testset "rates" begin
        for e in 1:glbNumOfExps()
            @test glbGetBinCentersListPtr(e) ≈ R["$name/$e/bincenters"] rtol = RTOL_EXACT
            for r in 1:glbGetNumberOfRules(e)
                k = "$name/$e/rule$r"
                # C leaves the rule rate outside the energy window uninitialised
                lo, hi = glbGetEnergyWindowBins(e, r)
                @test glbGetRuleRatePtr(e, r)[lo:hi] ≈ R["$k/rate"][lo:hi] rtol = RTOL_EXACT
                @test glbGetSignalRatePtr(e, r) ≈ R["$k/signal"] rtol = RTOL_EXACT
                @test glbGetBGRatePtr(e, r) ≈ R["$k/bg"] rtol = RTOL_EXACT
                @test glbTotalRuleRate(e, r, GLB_ALL, GLB_W_EFF, GLB_W_BG, GLB_W_COEFF, GLB_SIG) ≈
                      R["$k/total"] rtol = RTOL_EXACT
            end
            for c in 1:glbGetNumberOfChannels(e)
                @test glbGetChannelRatePtr(e, c, GLB_POST) ≈ R["$name/$e/channel$c/post"] rtol = RTOL_EXACT
            end
        end
    end

    @testset "chi2 without fit" begin
        @test glbChiSys(pt, GLB_ALL, GLB_ALL) ≈ R["$name/chisys"] rtol = RTOL_SYS atol = ATOL_SYS
        for e in 1:glbNumOfExps(), r in 1:glbGetNumberOfRules(e)
            @test glbChiSys(pt, e, r) ≈ R["$name/$e/rule$r/chisys"] rtol = RTOL_SYS atol = ATOL_SYS
        end
        glbSwitchSystematics(GLB_ALL, GLB_ALL, GLB_OFF)
        @test glbChiSys(pt, GLB_ALL, GLB_ALL) ≈ R["$name/chisys_off"] rtol = RTOL_SYS atol = ATOL_SYS
        glbSwitchSystematics(GLB_ALL, GLB_ALL, GLB_ON)
    end

    @testset "chi2 with fit" begin
        @test glbChiTheta13(pt, fit, GLB_ALL) ≈ R["$name/chitheta13"] rtol = RTOL_FIT atol = ATOL_FIT
        @test glbChiDelta(pt, fit, GLB_ALL) ≈ R["$name/chidelta"] rtol = RTOL_FIT atol = ATOL_FIT
        @test params_vector(fit) ≈ R["$name/chidelta_fit"] rtol = RTOL_FIT_PARAMS
        @test glbChiTheta13Delta(pt, nothing, GLB_ALL) ≈ R["$name/chitheta13delta"] rtol = RTOL_FIT atol = ATOL_FIT

        proj = glbAllocProjection()
        glbDefineProjection(proj, GLB_FIXED, GLB_FIXED, GLB_FREE, GLB_FIXED, GLB_FIXED, GLB_FREE)
        glbSetDensityProjectionFlag(proj, GLB_FREE, GLB_ALL)
        glbSetProjection(proj)
        @test glbChiNP(pt, nothing, GLB_ALL) ≈ R["$name/chinp"] rtol = RTOL_FIT atol = ATOL_FIT
    end

    if name != "DChooz"
        @testset "built-in chi2 functions" begin
            for (chi, errs) in ALT_CHI
                glbSetChiFunction(GLB_ALL, GLB_ALL, GLB_ON, chi, errs)
                rtol = chi == "chiSpectrumCalib" ? RTOL_CALIB : RTOL_SYS
                @test glbChiSys(pt, GLB_ALL, GLB_ALL) ≈ R["$name/alt/$chi"] rtol = rtol atol = ATOL_SYS
            end
            # At fixed nuisance parameters, with the fit rates of the last glbChiSys
            x = [0.01, -0.005, 0.02, 0.003]
            for e in 1:glbNumOfExps(), r in 1:glbGetNumberOfRules(e), (chi, errs) in ALT_CHI
                f = glbGetChiFunctionPtr(chi)
                @test f(e, r, length(errs), x, errs, nothing) ≈ R["$name/$e/rule$r/fixed/$chi"] rtol = RTOL_EXACT atol = ATOL_SYS
            end
        end
    end

    @testset "exposure" begin
        glbSetTargetMass(1, 2.0 * glbGetTargetMass(1))
        glbSetRunningTime(1, 1, 0.5 * glbGetRunningTime(1, 1))
        glbSetRates()
        for r in 1:glbGetNumberOfRules(1)
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
            probs(f) = [f(a, b, cp, E) for cp in (+1, -1) for E in energies for a in 1:3 for b in 1:3]
            @test probs((a, b, cp, E) -> glbVacuumProbability(a, b, cp, E, 1300.0)) ≈
                  R["prob/vacuum"] rtol = RTOL_EXACT atol = 1e-14
            @test probs((a, b, cp, E) -> glbConstantDensityProbability(a, b, cp, E, 1300.0, 2.8)) ≈
                  R["prob/constant"] rtol = RTOL_EXACT atol = 1e-14
            @test probs((a, b, cp, E) -> glbProfileProbability(1, a, b, cp, E)) ≈
                  R["prob/profile"] rtol = RTOL_EXACT atol = 1e-14

            # Unitarity
            for E in energies, cp in (+1, -1), a in 1:3
                @test sum(glbConstantDensityProbability(a, b, cp, E, 1300.0, 2.8) for b in 1:3) ≈ 1 atol = 1e-12
            end
        end

        @testset "energy scale shift" begin
            spectrum = [100.0 * exp(-0.5 * ((i - 8.0) / 3.0)^2) + i for i in 0:19]
            for (k, b) in enumerate((-0.3, -0.04, 0.01, 0.07, 0.3))
                @test glbShiftEnergyScale(b, spectrum, 20, 0.5, 10.5) ≈ R["shift/$k"] rtol = RTOL_EXACT
            end
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
            prob = (P, cp, E, psteps, L, rho, fs, ud) -> (calls[] += 1;
                                                          GLoBES.default_prob_func(P, cp, E, psteps, L, rho, fs, ud))
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
    end
end
