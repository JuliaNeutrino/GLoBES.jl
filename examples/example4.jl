# GLoBES example 4 (port of example4.c)
#
# Bar charts for the θ13 sensitivity: statistics only, with systematics, with
# correlations, and with correlations and the sgn(dm31²)-degeneracy, each for
# two values of dm31².
#
# Run with `julia --project=<path to GLoBES.jl> example4.jl`; the output goes
# to test4a.dat ... test4d.dat in the current directory.

using GLoBES
include("myio.jl")

# If filenames given, write to file; for empty filenames, write to screen
const MYFILE1 = "test4a.dat"
const MYFILE2 = "test4b.dat"
const MYFILE3 = "test4c.dat"
const MYFILE4 = "test4d.dat"

# Global oscillation parameters
const theta12 = asin(sqrt(0.8)) / 2
const theta13 = 0.0
const theta23 = π / 4
const deltacp = π / 2
const sdm = 7e-5
const ldm = 2e-3

# Simple sign-function
Sign(n) = n >= 0 ? +1.0 : -1.0

# Set simulated rates
function SetOscParams(theldm)
    true_values = glbAllocParams()
    glbDefineParams(true_values, theta12, theta13, theta23, deltacp, sdm, theldm)
    glbSetDensityParams(true_values, 1.0, GLB_ALL)
    glbSetOscillationParameters(true_values)
    glbSetRates()
end

# Calculate χ² with systematics
function CalcSystematics(theldm, thex)
    SetOscParams(theldm)
    thetheta13 = asin(sqrt(10^thex)) / 2
    test_values = glbAllocParams()
    glbDefineParams(test_values, theta12, thetheta13, theta23, 0.0, sdm, theldm)
    glbSetDensityParams(test_values, 1.0, GLB_ALL)
    return glbChiSys(test_values, GLB_ALL, GLB_ALL)
end

# Calculate χ² with statistics only
function CalcNoSystematics(theldm, thex)
    glbSwitchSystematics(GLB_ALL, GLB_ALL, GLB_OFF)
    res = CalcSystematics(theldm, thex)
    glbSwitchSystematics(GLB_ALL, GLB_ALL, GLB_ON)
    return res
end

# Calculate χ² with correlations
function CalcProjection(theldm, thex, start_vector)
    SetOscParams(theldm)

    thetheta13 = asin(sqrt(10^thex)) / 2
    thesign = Sign(glbGetOscParams(start_vector, GLB_DM_ATM))

    input_errors = glbAllocParams()
    central_values = glbAllocParams()
    minimum = glbAllocParams()

    glbDefineParams(input_errors, theta12 * 0.1, 0, 0, 0, sdm * 0.1, theldm / 3)
    glbSetDensityParams(input_errors, 0.05, GLB_ALL)
    # Set central values to respective ±dm31² to avoid falling into unwanted solution;
    # note that the error on dm31² should not be too small in order to avoid large
    # contributions from the priors
    glbDefineParams(central_values, theta12, theta13, theta23, deltacp, sdm, theldm * thesign)
    glbSetDensityParams(central_values, 1.0, GLB_ALL)
    glbSetCentralValues(central_values)
    glbSetInputErrors(input_errors)

    glbSetOscParams(start_vector, thetheta13, GLB_THETA_13)
    res = glbChiTheta13(start_vector, minimum, GLB_ALL)

    # Trick for next run: find position of minimum better by using deltacp of prior
    # calculation; note that this modifies the vectors in the main routine!
    glbSetOscParams(start_vector, glbGetOscParams(minimum, GLB_DELTA_CP), GLB_DELTA_CP)
    return res
end

# Find degeneracy
function FindDeg(deg_pos, theldm)
    SetOscParams(theldm)

    input_errors = glbAllocParams()
    central_values = glbAllocParams()

    glbDefineParams(central_values, theta12, theta13, theta23, deltacp, sdm, -theldm)
    glbSetDensityParams(central_values, 1.0, GLB_ALL)
    glbDefineParams(input_errors, theta12 * 0.1, 0, 0, 0, sdm * 0.1, theldm / 3)
    glbSetDensityParams(input_errors, 0.05, GLB_ALL)
    glbSetCentralValues(central_values)
    glbSetInputErrors(input_errors)
    return glbChiAll(central_values, deg_pos, GLB_ALL)
end

function main()
    # Initialize experiment NFstandard.glb
    glbInitExperiment(joinpath(@__DIR__, "NFstandard.glb"))

    original1 = glbAllocParams()
    original2 = glbAllocParams()
    degeneracy1 = glbAllocParams()
    degeneracy2 = glbAllocParams()

    # Compute 1st edges of bars: systematics off, fit value of deltacp=0
    InitOutput(MYFILE1, "Format: Log(10,s22th13)  chi^2:ldm=2e-3  chi^2:ldm=3.0e-3 \n")
    x = -7.0
    while x < -2.0 + 0.01
        res1 = CalcNoSystematics(2.0e-3, x)
        res2 = CalcNoSystematics(3.0e-3, x)
        AddToOutput(x, res1, res2)
        x += 5.0 / 50
    end

    # Compute 2nd edges of bars: systematics on, fit value of deltacp=0
    InitOutput(MYFILE2, "Format: Log(10,s22th13)  chi^2:ldm=2e-3  chi^2:ldm=3.0e-3 \n")
    x = -7.0
    while x < -2.0 + 0.01
        res1 = CalcSystematics(2.0e-3, x)
        res2 = CalcSystematics(3.0e-3, x)
        AddToOutput(x, res1, res2)
        x += 5.0 / 50
    end

    # Compute 3rd edges of bars: systematics on, projection onto s22th13 axis
    InitOutput(MYFILE3, "Format: Log(10,s22th13)  chi^2:ldm=2e-3  chi^2:ldm=3.0e-3 \n")
    glbDefineParams(original1, theta12, theta13, theta23, deltacp, sdm, 2.0e-3)
    glbSetDensityParams(original1, 1.0, GLB_ALL)
    glbDefineParams(original2, theta12, theta13, theta23, deltacp, sdm, 3.0e-3)
    glbSetDensityParams(original2, 1.0, GLB_ALL)
    x = -7.0
    while x < -2.0 + 0.01
        res1 = CalcProjection(2.0e-3, x, original1)
        res2 = CalcProjection(3.0e-3, x, original2)
        AddToOutput(x, res1, res2)
        x += 5.0 / 50
    end

    # Find sgn-degeneracies
    FindDeg(degeneracy1, 2.0e-3)
    FindDeg(degeneracy2, 3.0e-3)

    # Compute 4th edges of bars: systematics on, degeneracy, projection onto s22th13 axis
    InitOutput(MYFILE4, "Format: Log(10,s22th13)  chi^2:ldm=2e-3  chi^2:ldm=3.0e-3 \n")
    x = -7.0
    while x < -2.0 + 0.01
        res1 = CalcProjection(2.0e-3, x, degeneracy1)
        res2 = CalcProjection(3.0e-3, x, degeneracy2)
        AddToOutput(x, res1, res2)
        x += 5.0 / 30
    end
end

main()
