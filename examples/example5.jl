# GLoBES example 5 (port of example5.c)
#
# Sensitivity to sin²2θ13 of a reactor experiment with near and far detectors
# (Double Chooz), as a function of running time, with user-defined χ²
# functions for different sets of systematical errors.
#
# Run with `julia --project=<path to GLoBES.jl> example5.jl`; the output goes
# to test5a.dat ... test5e.dat in the current directory.
#
# The C version uses the GSL Brent root finder; `brent_root_*` below is a
# port of it, so that both versions iterate identically.

using GLoBES
include("myio.jl")

# If filenames given, write to file; for empty filenames, write to screen
const MYFILE1 = "test5a.dat"
const MYFILE2 = "test5b.dat"
const MYFILE3 = "test5c.dat"
const MYFILE4 = "test5d.dat"
const MYFILE5 = "test5e.dat"

# "True" oscillation parameters (cf. hep-ph/0405172v5)
const theta12 = asin(sqrt(0.3))
const theta13 = 0.0
const theta23 = π / 4
const deltacp = π / 2
const sdm = 7.9e-5
const ldm = 2.6e-3

# Experiment and analysis parameters
const min_runtime = 1e-2           # Minimum running time to consider [years]
const max_runtime = 1e4            # Maximum running time to consider [years]
const chi2_goal = 2.7              # Desired χ² value in root finder (2.7 = 90% C.L.)
const tSteps = 30                  # Number of data points for each curve
const logs22th13_precision = 0.005 # Desired precision of log(sin[th13]^2) in root finder

const MAX_SYS = 200
const sys_errors = zeros(MAX_SYS)    # Uncertainties of systematics params
const sys_startval = zeros(MAX_SYS)  # Starting values for systematics minimizer
const sigma_binbin = Ref(0.0)        # Bin-to-bin error

const EXP_FAR = 1
const EXP_NEAR = 2

# GLoBES parameter structures, allocated in main() after the experiments are loaded
const true_values = Ref{GlbParams}()
const test_values = Ref{GlbParams}()
const input_errors = Ref{GlbParams}()

#  H E L P E R   F U N C T I O N S

# Gauss likelihood (this is sufficient here due to the large event numbers
# in a reactor experiment; for other setups, one should use Poisson statistics)
likelihood(true_rate, fit_rate, sqr_sigma) =
    sqr_sigma > 0 ? (true_rate - fit_rate)^2 / sqr_sigma : 0.0

#  C H I ^ 2   F U N C T I O N S

# χ² for Double Chooz, including the following systematical errors:
#   x[1]: Flux normalization of reactor
#   x[2]: Fiducial mass error - far detector
#   x[3]: Fiducial mass error - near detector
#   x[4]: Energy calibration error - far detector
#   x[5]: Energy calibration error - near detector
function chiDCNorm(exp, rule, n_params, x, errors, user_data)
    n_bins = glbGetNumberOfBins(EXP_FAR)
    true_rates_N = glbGetRuleRatePtr(EXP_NEAR, 1)
    true_rates_F = glbGetRuleRatePtr(EXP_FAR, 1)

    # Request simulated energy interval and analysis energy window
    emin, emax = glbGetEminEmax(exp)
    ew_low, ew_high = glbGetEnergyWindowBins(exp, rule)

    # Apply energy calibration error
    signal_fit_rates_F =
        glbShiftEnergyScale(x[4], glbGetSignalFitRatePtr(EXP_FAR, 1), n_bins, emin, emax)
    signal_fit_rates_N =
        glbShiftEnergyScale(x[5], glbGetSignalFitRatePtr(EXP_NEAR, 1), n_bins, emin, emax)

    # Loop over all bins in energy window
    signal_norm_F = 1.0 + x[1] + x[2]
    signal_norm_N = 1.0 + x[1] + x[3]
    chi2 = 0.0
    for i = ew_low:ew_high
        # Statistical part of χ² for far detector
        # Normalization is affected by flux error x[1] and fiducial mass error x[2]
        fit_rate = signal_norm_F * signal_fit_rates_F[i]
        chi2 += likelihood(true_rates_F[i], fit_rate, true_rates_F[i])

        # Statistical part of χ² for near detector
        # Normalization is affected by flux error x[1] and fiducial mass error x[3]
        fit_rate = signal_norm_N * signal_fit_rates_N[i]
        chi2 += likelihood(true_rates_N[i], fit_rate, true_rates_N[i])
    end

    # Systematical part of χ² (= priors)
    for i = 1:n_params
        chi2 += (x[i] / errors[i])^2
    end

    # Save the systematics parameters as starting values for the next step
    for i = 1:n_params
        sys_startval[i] = x[i]
    end
    return chi2
end

# χ² for Double Chooz, including the following systematical errors:
#   x[1]: Flux normalization of reactor
#   x[2]: Fiducial mass error - far detector
#   x[3]: Fiducial mass error - near detector
#   x[4]: Energy calibration error - far detector
#   x[5]: Energy calibration error - near detector
#   x[6]...x[6+nBins-1]: Spectral error
# and the bin-to-bin error sigma_binbin. The calculation is based on
# Eq. (9) of hep-ph/0303232.
function chiDCSpectral(exp, rule, n_params, x, errors, user_data)
    n_bins = glbGetNumberOfBins(EXP_FAR)
    true_rates_N = glbGetRuleRatePtr(EXP_NEAR, 1)
    true_rates_F = glbGetRuleRatePtr(EXP_FAR, 1)

    # Request simulated energy interval and analysis energy window
    emin, emax = glbGetEminEmax(exp)
    ew_low, ew_high = glbGetEnergyWindowBins(exp, rule)

    # Apply energy calibration error
    signal_fit_rates_F =
        glbShiftEnergyScale(x[4], glbGetSignalFitRatePtr(EXP_FAR, 1), n_bins, emin, emax)
    signal_fit_rates_N =
        glbShiftEnergyScale(x[5], glbGetSignalFitRatePtr(EXP_NEAR, 1), n_bins, emin, emax)

    # Loop over all bins in energy window
    signal_norm_F = 1.0 + x[1] + x[2]
    signal_norm_N = 1.0 + x[1] + x[3]
    s2 = sigma_binbin[]^2
    chi2 = 0.0
    for i = ew_low:ew_high
        # Statistical part of χ² for far detector
        # Normalization is affected by flux error x[1], fiducial mass error x[2],
        # spectral perturbations x[5+i], and bin-to-bin error sigma_binbin
        fit_rate = (signal_norm_F + x[5 + i]) * signal_fit_rates_F[i]
        chi2 +=
            likelihood(true_rates_F[i], fit_rate, true_rates_F[i] * (1.0 + true_rates_F[i] * s2))

        # Statistical part of χ² for near detector
        # Normalization is affected by flux error x[1], fiducial mass error x[3],
        # spectral perturbations x[5+i], and bin-to-bin error sigma_binbin
        fit_rate = (signal_norm_N + x[5 + i]) * signal_fit_rates_N[i]
        chi2 +=
            likelihood(true_rates_N[i], fit_rate, true_rates_N[i] * (1.0 + true_rates_N[i] * s2))
    end

    # Systematical part of χ² (= priors)
    for i = 1:n_params
        chi2 += (x[i] / errors[i])^2
    end

    # Save the systematics parameters as starting values for the next step
    for i = 1:n_params
        sys_startval[i] = x[i]
    end
    return chi2
end

#  B R E N T   R O O T   F I N D E R   (as gsl_root_fsolver_brent)

mutable struct BrentSolver
    f::Function
    a::Float64
    b::Float64
    c::Float64
    d::Float64
    e::Float64
    fa::Float64
    fb::Float64
    fc::Float64
    root::Float64
    x_lower::Float64
    x_upper::Float64
    BrentSolver(f) = new(f, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
end

# Returns false if the endpoints do not straddle zero (GSL_EINVAL)
function brent_root_set!(s::BrentSolver, x_lower, x_upper)
    s.root = 0.5 * (x_lower + x_upper)
    s.x_lower = x_lower
    s.x_upper = x_upper
    f_lower = s.f(x_lower)
    f_upper = s.f(x_upper)
    s.a = x_lower
    s.fa = f_lower
    s.b = x_upper
    s.fb = f_upper
    s.c = x_upper
    s.fc = f_upper
    s.d = x_upper - x_lower
    s.e = x_upper - x_lower
    return !((f_lower < 0.0 && f_upper < 0.0) || (f_lower > 0.0 && f_upper > 0.0))
end

function brent_root_iterate!(s::BrentSolver)
    a, b, c, d, e = s.a, s.b, s.c, s.d, s.e
    fa, fb, fc = s.fa, s.fb, s.fc
    ac_equal = false

    if (fb < 0 && fc < 0) || (fb > 0 && fc > 0)
        ac_equal = true
        c = a
        fc = fa
        d = b - a
        e = b - a
    end
    if abs(fc) < abs(fb)
        ac_equal = true
        a = b
        b = c
        c = a
        fa = fb
        fb = fc
        fc = fa
    end

    tol = 0.5 * eps(Float64) * abs(b)
    m = 0.5 * (c - b)

    if fb == 0
        s.root = b
        s.x_lower = b
        s.x_upper = b
        return
    end
    if abs(m) <= tol
        s.root = b
        s.x_lower, s.x_upper = b < c ? (b, c) : (c, b)
        return
    end

    if abs(e) < tol || abs(fa) <= abs(fb)
        d = m
        e = m                      # use bisection
    else
        s_ = fb / fa                      # use inverse cubic interpolation
        if ac_equal
            p = 2 * m * s_
            q = 1 - s_
        else
            q = fa / fc
            r = fb / fc
            p = s_ * (2 * m * q * (q - r) - (b - a) * (r - 1))
            q = (q - 1) * (r - 1) * (s_ - 1)
        end
        if p > 0
            q = -q
        else
            p = -p
        end
        if 2 * p < min(3 * m * q - abs(tol * q), abs(e * q))
            e = d
            d = p / q
        else
            d = m
            e = m                  # interpolation failed, fall back to bisection
        end
    end

    a = b
    fa = fb
    b += abs(d) > tol ? d : (m > 0 ? +tol : -tol)
    fb = s.f(b)

    s.a, s.b, s.c, s.d, s.e = a, b, c, d, e
    s.fa, s.fb, s.fc = fa, fb, fc

    # Update the best estimate of the root and bounds on each iteration
    s.root = b
    if (fb < 0 && fc < 0) || (fb > 0 && fc > 0)
        c = a
    end
    s.x_lower, s.x_upper = b < c ? (b, c) : (c, b)
    return
end

# As gsl_root_test_interval
function root_test_interval(x_lower, x_upper, epsabs, epsrel)
    abs_lower, abs_upper = abs(x_lower), abs(x_upper)
    min_abs =
        (x_lower > 0 && x_upper > 0) || (x_lower < 0 && x_upper < 0) ? min(abs_lower, abs_upper) :
        0.0
    return abs(x_upper - x_lower) < epsabs + epsrel * min_abs
end

#  W R A P P E R   F O R   R O O T   F I N D E R

function DoChiSquare(x)
    # Set vector of test values
    thetheta13 = asin(sqrt(10^x)) / 2
    glbSetOscParams(test_values[], thetheta13, GLB_THETA_13)
    glbSetOscParams(test_values[], 200.0 / 2 * (x + 4) * π / 180, GLB_DELTA_CP)

    # Set starting values for systematics minimizer to the coordinates of
    # minimum in the last iteration. This accelerates the minimization and
    # prevents convergence problems.
    glbSetSysStartingValuesList(EXP_FAR, 1, GLB_ON, sys_startval)

    # Compute χ² for all loaded experiments and all rules
    # Correlations are unimportant in reactor experiments, so glbChiSys is sufficient
    chi2 = glbChiSys(test_values[], GLB_ALL, GLB_ALL)
    return chi2 - chi2_goal
end

# Compute sensitivity to sin²(2θ13) for different running times
function ComputeSensitivityCurve(s::BrentSolver)
    t_factor = (max_runtime / min_runtime)^(1.0 / tSteps)
    max_iter = 100                   # Maximum number of iterations allowed in root finder
    x_sens = ones(tSteps + 1)        # As x is always <= 0, 1 signals "no value"

    # Loop over all data points in this dataset, using log stepping
    for j = 0:tSteps
        # Set running time
        t = min_runtime * t_factor^j
        glbSetRunningTime(EXP_FAR, 1, t)
        glbSetRunningTime(EXP_NEAR, 1, t)

        # Calculate "true" event rates
        glbDefineParams(true_values[], theta12, theta13, theta23, deltacp, sdm, ldm)
        glbSetDensityParams(true_values[], 1.0, GLB_ALL)
        glbDefineParams(test_values[], theta12, theta13, theta23, deltacp, sdm, ldm)
        glbSetDensityParams(test_values[], 1.0, GLB_ALL)
        glbDefineParams(input_errors[], 0.1 * theta12, 0, 0.15 * theta23, 0, 0.05 * sdm, 0.05 * ldm)
        glbSetDensityParams(input_errors[], 0.05, GLB_ALL)
        glbSetOscillationParameters(true_values[])
        glbSetInputErrors(input_errors[])
        glbSetRates()

        # Determine sensitivity to sin²(2θ13) by using the Brent-Dekker
        # algorithm to find the θ13 for which χ² = 2.7 (90%)
        x_lo = -2.0
        x_hi = -0.3
        iter = 0

        # Start root finder. The initial search interval is guessed based on the
        # results from a neighboring grid point
        deviation = 0.1
        while true
            if j > 0
                x_lo = x_sens[j] - deviation
                x_hi = min(x_sens[j] + deviation, -0.001)
            end
            ok = brent_root_set!(s, x_lo, x_hi)
            deviation *= 2
            ok && break
        end

        # Iterate root finder (like the C version, iter is never incremented)
        x = 0.0
        while true
            brent_root_iterate!(s)
            x = s.root
            x_lo = s.x_lower
            x_hi = s.x_upper
            root_test_interval(x_lo, x_hi, logs22th13_precision, 0) && break
            iter < max_iter || break
        end

        # Save results
        x_sens[j + 1] = x
        AddToOutput2(t, x)
    end
end

#  M A I N   P R O G R A M

function main()
    n_bins = 62                       # Number of bins

    # Define χ² functions
    glbDefineChiFunction(chiDCNorm, 5, "chiDCNorm")
    glbDefineChiFunction(chiDCSpectral, 5 + n_bins, "chiDCSpectral")

    s = BrentSolver(DoChiSquare)

    # Load 2 experiments: DC far (#1) and near (#2) detectors
    glbClearExperimentList()
    glbInitExperiment(joinpath(@__DIR__, "D-Chooz_far.glb"))
    glbInitExperiment(joinpath(@__DIR__, "D-Chooz_near.glb"))
    if glbGetNumberOfBins(EXP_FAR) != n_bins || glbGetNumberOfBins(EXP_NEAR) != n_bins
        println("ERROR: Number of bins changed in AEDL file, but not in the code (or vice-versa).")
        return -1
    end

    # Initialize parameter vectors
    true_values[] = glbAllocParams()
    test_values[] = glbAllocParams()
    input_errors[] = glbAllocParams()
    glbDefineParams(true_values[], theta12, theta13, theta23, deltacp, sdm, ldm)
    glbSetDensityParams(true_values[], 1.0, GLB_ALL)
    glbDefineParams(test_values[], theta12, theta13, theta23, deltacp, sdm, ldm)
    glbSetDensityParams(test_values[], 1.0, GLB_ALL)
    glbDefineParams(input_errors[], 0.1 * theta12, 0, 0.15 * theta23, 0, 0.05 * sdm, 0.05 * ldm)
    glbSetDensityParams(input_errors[], 0.05, GLB_ALL)
    glbSetOscillationParameters(true_values[])
    glbSetInputErrors(input_errors[])

    # Calculate sensitivity curve without systematics
    InitOutput(MYFILE1, "Format: Running time   Log(10,s22th13) sens. \n")
    glbSwitchSystematics(GLB_ALL, GLB_ALL, GLB_OFF)
    ComputeSensitivityCurve(s)
    glbSwitchSystematics(GLB_ALL, GLB_ALL, GLB_ON)

    # Calculate sensitivity curve with normalization and energy calibration errors,
    # as defined in the AEDL files
    InitOutput(MYFILE2, "Format: Running time   Log(10,s22th13) sens. \n")
    ComputeSensitivityCurve(s)

    # Calculate sensitivity curve with the above + spectral error.
    # Since chiDCSpectral computes the complete χ² for the whole problem, it
    # must be called only for ONE of the two experiments
    InitOutput(MYFILE3, "Format: Running time   Log(10,s22th13) sens. \n")
    old_sys_errors = glbGetSysErrorsListPtr(EXP_FAR, 1, GLB_ON)   # Fill error array
    sys_dim = glbGetSysDimInExperiment(EXP_FAR, 1, GLB_ON)
    for i = 1:sys_dim                 # Normalization and energy calibration errors
        sys_errors[i] = old_sys_errors[i]
    end
    for i = (sys_dim + 1):(sys_dim + n_bins)
        sys_errors[i] = 0.02           # Spectral error
    end
    sigma_binbin[] = 0.0               # No bin-to-bin error for the moment
    glbSetChiFunction(EXP_FAR, 1, GLB_ON, "chiDCSpectral", sys_errors)
    glbSetChiFunction(EXP_NEAR, 1, GLB_ON, "chiZero", sys_errors)
    ComputeSensitivityCurve(s)

    # Calculate sensitivity curve with the above + bin-to-bin error
    # (uncorrelated between near and far detectors)
    InitOutput(MYFILE4, "Format: Running time   Log(10,s22th13) sens. \n")
    sigma_binbin[] = 0.005
    ComputeSensitivityCurve(s)

    # Calculate sensitivity curve with a different bin-to-bin error
    InitOutput(MYFILE5, "Format: Running time   Log(10,s22th13) sens. \n")
    sigma_binbin[] = 0.02
    ComputeSensitivityCurve(s)
    return 0
end

main()
