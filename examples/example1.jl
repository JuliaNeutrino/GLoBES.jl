# GLoBES example 1 (port of example1.c)
#
# Correlation between s22th13 and deltacp: χ² with systematics on a grid
# in log10(sin²2θ13) and δCP.
#
# Run from any directory with
#     julia --project=<path to GLoBES.jl> example1.jl
# The output goes to test1.dat in the current directory.

using GLoBES
include("myio.jl")

# If filename given, write to file; for empty filename write to screen
const MYFILE = "test1.dat"

function main()
    # Initialize experiment NFstandard.glb
    glbInitExperiment(joinpath(@__DIR__, "NFstandard.glb"))

    # Initialize output
    InitOutput(MYFILE, "Format: Log(10,s22th13)   deltacp   chi^2 \n")

    # Define standard oscillation parameters
    theta12 = asin(sqrt(0.8)) / 2
    theta13 = asin(sqrt(0.001)) / 2
    theta23 = π / 4
    deltacp = π / 2
    sdm = 7e-5
    ldm = 2e-3

    # Initialize parameter vector(s)
    true_values = glbAllocParams()
    test_values = glbAllocParams()

    glbDefineParams(true_values, theta12, theta13, theta23, deltacp, sdm, ldm)
    glbSetDensityParams(true_values, 1.0, GLB_ALL)
    glbDefineParams(test_values, theta12, theta13, theta23, deltacp, sdm, ldm)
    glbSetDensityParams(test_values, 1.0, GLB_ALL)

    # The simulated data are computed
    glbSetOscillationParameters(true_values)
    glbSetRates()

    # Iteration over all values to be computed
    x = -4.0
    while x < -2.0 + 0.01
        y = 0.0
        while y < 200.0 + 0.01
            # Set vector of test values
            thetheta13 = asin(sqrt(10^x)) / 2
            glbSetOscParams(test_values, thetheta13, GLB_THETA_13)
            glbSetOscParams(test_values, y * π / 180.0, GLB_DELTA_CP)

            # Compute χ² for all loaded experiments and all rules
            res = glbChiSys(test_values, GLB_ALL, GLB_ALL)
            AddToOutput(x, y, res)
            y += 200.0 / 50
        end
        x += 2.0 / 50
    end
end

main()
