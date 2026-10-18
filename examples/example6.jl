# GLoBES example 6 (port of example6.c)
#
# Non-standard interactions and user-defined probability engines: wave packet
# decoherence in a reactor experiment, parameterized by the energy spread σ_E
# as a 7th oscillation parameter. Similar to Chapter 4 of hep-ph/0502147.
#
# Run with `julia --project=<path to GLoBES.jl> example6.jl`; the output goes
# to test6a.dat and test6b.dat in the current directory.

using GLoBES
include("myio.jl")

const GLB_SIGMA_E = 7        # Index of non-standard parameter sigma_E

# If filename given, write to file; for empty filename write to screen
const MYFILE1 = "test6a.dat"
const MYFILE2 = "test6b.dat"

#  U S E R - D E F I N E D   P R O B A B I L I T Y   E N G I N E

# Parameters of the engine, in GeV² for the mass splittings
mutable struct EngineParams
    th12::Float64
    th13::Float64
    th23::Float64
    deltacp::Float64
    sdm::Float64
    ldm::Float64
    sigma_E::Float64
end
const engine = EngineParams(0, 0, 0, 0, 0, 0, 0)

# Store oscillation parameters in internal data structures.
# For more sophisticated probability engines, this would be the right
# place to pre-compute the mixing matrix and parts of the Hamiltonian in
# order to speed up the calls to the actual probability matrix function.
function my_set_oscillation_parameters(p, user_data)
    engine.th12 = glbGetOscParams(p, GLB_THETA_12)
    engine.th13 = glbGetOscParams(p, GLB_THETA_13)
    engine.th23 = glbGetOscParams(p, GLB_THETA_23)
    engine.deltacp = glbGetOscParams(p, GLB_DELTA_CP)
    engine.sdm = glbGetOscParams(p, GLB_DM_21) * 1.0e-18   # Convert to GeV²
    engine.ldm = glbGetOscParams(p, GLB_DM_31) * 1.0e-18   # Convert to GeV²
    engine.sigma_E = glbGetOscParams(p, GLB_SIGMA_E)
    return 0
end

# Write oscillation parameters from internal data structures into p.
function my_get_oscillation_parameters(p, user_data)
    glbSetOscParams(p, engine.th12, GLB_THETA_12)
    glbSetOscParams(p, engine.th13, GLB_THETA_13)
    glbSetOscParams(p, engine.th23, GLB_THETA_23)
    glbSetOscParams(p, engine.deltacp, GLB_DELTA_CP)
    glbSetOscParams(p, engine.sdm * 1.0e18, GLB_DM_21)  # Convert to eV²
    glbSetOscParams(p, engine.ldm * 1.0e18, GLB_DM_31)  # Convert to eV²
    glbSetOscParams(p, engine.sigma_E, GLB_SIGMA_E)
    return 0
end

# Calculate oscillation probabilities.
# Since for our setup, only P_ee is required, all other entries of P are
# set to zero for simplicity. Furthermore, we neglect matter effects and
# the filter feature (parameter filter_sigma).
# The formula for P_ee is Eq. (36) from hep-ph/0502147.
#
# Parameters:
#   P:            The 3×3 matrix where the probabilities are to be stored
#   cp_sign:      +1 if probalities for neutrinos are requested, -1 for
#                 anti-neutrinos.
#   E:            The neutrino energy in GeV
#   psteps:       Number of constant density layers in the matter profile
#   length:       The lengths of these layers in km
#   density:      The individual densities of these layers in g/cm³
#   filter_sigma: Width of low-pass filter as given in the AEDL file
function my_probability_matrix(P, cp_sign, E, psteps, length, density, filter_sigma, user_data)
    # Set all probabilities to zero initially
    fill!(P, 0.0)

    # Calculate total baseline
    L = 0.0
    for i = 1:psteps
        L += length[i]
    end
    L = GLB_KM_TO_EV(L) * 1.0e9      # Convert to GeV⁻¹

    # Compute P_ee
    s12, c12 = sin(engine.th12), cos(engine.th12)
    s13, c13 = sin(engine.th13), cos(engine.th13)
    t = L / (4.0 * E)
    Delta21 = engine.sdm * t
    Delta31 = engine.ldm * t
    Delta32 = Delta31 - Delta21
    t = sqrt(2) * engine.sigma_E / E
    D21 = exp(-(Delta21 * t)^2)
    D31 = exp(-(Delta31 * t)^2)
    D32 = exp(-(Delta32 * t)^2)
    P[1, 1] =
        c13^4 * (1 - 2.0 * (s12 * c12)^2 * (1 - D21 * cos(2.0 * Delta21))) +
        2.0 *
        (s13 * c13)^2 *
        (D31 * c12^2 * cos(2.0 * Delta31) + D32 * s12^2 * cos(2.0 * Delta32)) +
        s13^4
    return 0
end

#  M A I N   P R O G R A M

function main()
    # Define standard oscillation parameters (cf. hep-ph/0405172v5)
    true_theta12 = asin(sqrt(0.3))
    true_theta13 = 0.0
    true_theta23 = π / 4
    true_deltacp = 0.0
    true_sdm = 7.9e-5
    true_ldm = 2.6e-3

    # Define one non-standard parameter sigma_E (wave packet energy spread
    # responsible for wave packet decoherence)
    true_sigma_E = 0.0

    # Register non-standard probability engine. This has to be done
    # before any calls to glbAllocParams() or glbAllocProjection()
    glbRegisterProbabilityEngine(
        7,      # Number of parameters
        my_probability_matrix,
        my_set_oscillation_parameters,
        my_get_oscillation_parameters,
        nothing,
    )

    # Initialize reactor experiment
    glbInitExperiment(joinpath(@__DIR__, "Reactor1.glb"))

    # Initialize parameter and projection vector(s)
    true_values = glbAllocParams()
    test_values = glbAllocParams()
    input_errors = glbAllocParams()
    th13_sigma_projection = glbAllocProjection()

    glbDefineParams(
        true_values,
        true_theta12,
        true_theta13,
        true_theta23,
        true_deltacp,
        true_sdm,
        true_ldm,
    )
    glbSetOscParams(true_values, true_sigma_E, GLB_SIGMA_E)   # Non-standard parameter
    glbSetDensityParams(true_values, 1.0, GLB_ALL)
    glbDefineParams(
        test_values,
        true_theta12,
        true_theta13,
        true_theta23,
        true_deltacp,
        true_sdm,
        true_ldm,
    )
    glbSetOscParams(test_values, true_sigma_E, GLB_SIGMA_E)   # Non-standard parameter
    glbSetDensityParams(test_values, 1.0, GLB_ALL)

    # The simulated data are computed
    glbSetOscillationParameters(true_values)
    glbSetRates()

    # Set starting values and input errors for all projections
    glbDefineParams(
        input_errors,
        0.1 * true_theta12,
        0.0,
        0.15 * true_theta23,
        0.0,
        0.05 * true_sdm,
        0.05 * true_ldm,
    )
    glbSetOscParams(input_errors, 0.0, GLB_SIGMA_E)      # Non-standard parameter
    glbSetDensityParams(input_errors, 0.05, GLB_ALL)
    glbSetCentralValues(true_values)
    glbSetInputErrors(input_errors)

    # Compute χ² without correlations
    InitOutput(MYFILE1, "Format: Log(10,s22th13)   sigma_E       chi^2 \n")
    x = 0.0
    while x < 0.05 + 0.001              # th13 loop
        y = 0.0
        while y < 0.010 + 0.001         # sigma_E loop
            # Set vector of test=fit values
            thetheta13 = asin(sqrt(x)) / 2.0
            glbSetOscParams(test_values, thetheta13, GLB_THETA_13)
            glbSetOscParams(test_values, y, GLB_SIGMA_E)

            # Compute χ² with systematics only
            res = glbChiSys(test_values, GLB_ALL, GLB_ALL)
            AddToOutput(x, y, res)
            y += 0.001
        end
        x += 0.005
    end

    # Set two-parameter projection onto s22th13 and sigma_E, marginalizing
    # over all other parameters except delta_CP and th23, which do not enter P_ee
    glbDefineProjection(
        th13_sigma_projection,
        GLB_FREE,
        GLB_FIXED,
        GLB_FIXED,
        GLB_FIXED,
        GLB_FREE,
        GLB_FREE,
    )
    glbSetDensityProjectionFlag(th13_sigma_projection, GLB_FIXED, GLB_ALL)
    glbSetProjectionFlag(th13_sigma_projection, GLB_FIXED, GLB_SIGMA_E)  # Non-standard parameter
    glbSetProjection(th13_sigma_projection)

    InitOutput(MYFILE2, "Format: Log(10,s22th13)   sigma_E       chi^2 \n")
    x = 0.0
    while x < 0.05 + 0.001              # th13 loop
        y = 0.0
        while y < 0.010 + 0.001         # sigma_E loop
            # Set vector of test=fit values
            thetheta13 = asin(sqrt(x)) / 2.0
            glbSetOscParams(test_values, thetheta13, GLB_THETA_13)
            glbSetOscParams(test_values, y, GLB_SIGMA_E)

            # Compute χ² with correlations
            res = glbChiNP(test_values, nothing, GLB_ALL)
            AddToOutput(x, y, res)
            y += 0.001
        end
        x += 0.005
    end
end

main()
