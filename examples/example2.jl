# GLoBES example 2 (port of example2.c)
#
# Projection of the two- and n-dimensional χ² manifold onto the s22th13 axis:
# once minimizing over δCP only, once over all parameters but θ13.
#
# Run with `julia --project=<path to GLoBES.jl> example2.jl`; the output goes
# to test2.dat in the current directory.

using GLoBES
include("myio.jl")

# If filename given, write to file; for empty filename write to screen
const MYFILE = "test2.dat"

function main()
    # Initialize experiment NFstandard.glb
    glbInitExperiment(joinpath(@__DIR__, "NFstandard.glb"))

    # Initialize output
    InitOutput(MYFILE, "Format: Log(10,s22th13)   chi^2 one param   chi^2 all params \n")

    # Define standard oscillation parameters
    theta12 = asin(sqrt(0.8)) / 2
    theta13 = asin(sqrt(0.001)) / 2
    theta23 = π / 4
    deltacp = π / 2
    sdm = 7e-5
    ldm = 2e-3

    # Initialize parameter and projection vector(s)
    true_values = glbAllocParams()
    test_values = glbAllocParams()
    input_errors = glbAllocParams()
    theta13_projection = glbAllocProjection()

    glbDefineParams(true_values, theta12, theta13, theta23, deltacp, sdm, ldm)
    glbSetDensityParams(true_values, 1.0, GLB_ALL)
    glbDefineParams(test_values, theta12, theta13, theta23, deltacp, sdm, ldm)
    glbSetDensityParams(test_values, 1.0, GLB_ALL)

    # The simulated data are computed
    glbSetOscillationParameters(true_values)
    glbSetRates()

    # Set starting values and input errors for all projections
    glbDefineParams(input_errors, theta12 * 0.1, 0, 0, 0, sdm * 0.1, 0)
    glbSetDensityParams(input_errors, 0.05, GLB_ALL)
    glbSetCentralValues(true_values)
    glbSetInputErrors(input_errors)

    # Set two-parameter projection onto s22th13-axis: only deltacp free!
    glbDefineProjection(theta13_projection, GLB_FIXED, GLB_FIXED, GLB_FIXED,
                        GLB_FREE, GLB_FIXED, GLB_FIXED)
    glbSetDensityProjectionFlag(theta13_projection, GLB_FIXED, GLB_ALL)
    glbSetProjection(theta13_projection)

    # Iteration over all values to be computed
    x = -4.0
    while x < -2.0 + 0.001
        # Set vector of test=fit values
        thetheta13 = asin(sqrt(10^x)) / 2
        glbSetOscParams(test_values, thetheta13, GLB_THETA_13)

        # Guess fit value for deltacp in order to safely find minimum
        glbSetOscParams(test_values, 200.0 / 2 * (x + 4) * π / 180, GLB_DELTA_CP)

        # Compute χ² for two-parameter correlation: minimize over deltacp only
        res1 = glbChiNP(test_values, nothing, GLB_ALL)

        # Compute χ² for full correlation: minimize over all but theta13
        res2 = glbChiTheta13(test_values, nothing, GLB_ALL)

        AddToOutput(x, res1, res2)
        x += 2.0 / 50
    end
end

main()
