# GLoBES example 3 (port of example3.c)
#
# Find the sgn(dm31²)-degeneracy and compute its shape with glbChiSys.
#
# Run with `julia --project=<path to GLoBES.jl> example3.jl`; the output goes
# to test3.dat in the current directory.

using GLoBES
using Printf
include("myio.jl")

# If filename given, write to file; for empty filename write to screen
const MYFILE = "test3.dat"

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
    central_values = glbAllocParams()
    input_errors = glbAllocParams()
    deg_pos = glbAllocParams()

    # The simulated data are computed
    glbDefineParams(true_values, theta12, theta13, theta23, deltacp, sdm, ldm)
    glbSetDensityParams(true_values, 1.0, GLB_ALL)
    glbSetOscillationParameters(true_values)
    glbSetRates()

    # Find sgn-degeneracy
    glbDefineParams(central_values, theta12, theta13, theta23, deltacp, sdm, -ldm)
    glbSetDensityParams(central_values, 1.0, GLB_ALL)
    glbDefineParams(input_errors, theta12 * 0.1, 0, 0, 0, sdm * 0.1, ldm / 3)
    glbSetDensityParams(input_errors, 0.05, GLB_ALL)
    glbSetCentralValues(central_values)
    glbSetInputErrors(input_errors)
    CL = glbChiAll(central_values, deg_pos, GLB_ALL)

    @printf("Position of degeneracy: s22th13=%g, deltacp=%g; Confidence level: %g \n",
            glbGetOscParams(deg_pos, GLB_THETA_13), glbGetOscParams(deg_pos, GLB_DELTA_CP), CL)

    # If degeneracy at low enough confidence level: compute section
    if CL < 9.0
        x = -4.0
        while x < -2.0 + 0.01
            y = 0.0
            while y < 200.0 + 0.01
                # Set vector of test values
                thetheta13 = asin(sqrt(10^x)) / 2
                glbSetOscParams(deg_pos, thetheta13, GLB_THETA_13)
                glbSetOscParams(deg_pos, y * π / 180.0, GLB_DELTA_CP)

                # Compute χ² for all loaded experiments and all rules
                res = glbChiSys(deg_pos, GLB_ALL, GLB_ALL)
                AddToOutput(x, y, res)
                y += 200.0 / 50
            end
            x += 2.0 / 50
        end
    end
end

main()
