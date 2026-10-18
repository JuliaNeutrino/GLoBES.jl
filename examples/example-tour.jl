# GLoBES tour (port of example-tour.c)
#
# A walk through the main features: probabilities, χ² with systematics,
# correlations, projections, degeneracies and a two-experiment setup.
#
# Run with `julia --project=<path to GLoBES.jl> example-tour.jl`; the output
# goes to gl-tour.dat in the current directory.

using GLoBES
using Printf

function main()
    MYFILE = "gl-tour.dat"   # if empty, write to screen, otherwise to file name given here
    stream = isempty(MYFILE) ? stdout : open(MYFILE, "w")

    # Define my standard oscillation parameters
    theta12 = asin(sqrt(0.8)) / 2
    theta13 = asin(sqrt(0.001)) / 2
    theta23 = π / 4
    deltacp = π / 2
    sdm = 7e-5
    ldm = 2e-3

    # Initialize one experiment NFstandard.glb
    glbInitExperiment(joinpath(@__DIR__, "NFstandard.glb"))

    # Initialize a number of parameter vector(s)
    true_values = glbAllocParams()
    fit_values = glbAllocParams()
    central_values = glbAllocParams()
    input_errors = glbAllocParams()
    minimum = glbAllocParams()

    # Assign standard oscillation parameters
    glbDefineParams(true_values, theta12, theta13, theta23, deltacp, sdm, ldm)
    glbSetDensityParams(true_values, 1.0, GLB_ALL)

    # The simulated data are computed with "true_values"
    glbSetOscillationParameters(true_values)
    glbSetRates()

    # Return some low level information
    @printf(stream, "\nOscillation probabilities in vacuum: ")
    for i = 1:3
        @printf(stream, "1->%i: %g ", i, glbVacuumProbability(1, i, +1, 50, 3000))
    end
    @printf(stream, "\nOscillation probabilities in matter: ")
    for i = 1:3
        @printf(stream, "1->%i: %g ", i, glbProfileProbability(1, 1, i, +1, 50))
    end
    @printf(stream, "\n\n")

    # Set test/fit values slightly off true values at s22th13=0.0015
    glbCopyParams(true_values, fit_values)
    glbSetOscParams(fit_values, asin(sqrt(0.0015)) / 2, GLB_THETA_13)

    # Compute χ² with systematics only for all experiments and rules
    chi2 = glbChiSys(fit_values, GLB_ALL, GLB_ALL)
    @printf(stream, "chi2 with systematics only: %g\n\n", chi2)

    # What would be obtain from the appearance channels only?
    chi2 = glbChiSys(fit_values, 1, 1)
    @printf(stream, "This we would have from the CP-even appearance channel only: %g\n\n", chi2)

    # What would be obtain from the appearance channels only?
    chi2 =
        glbChiSys(fit_values, GLB_ALL, 1) +
        glbChiSys(fit_values, GLB_ALL, 2) +
        glbChiSys(fit_values, GLB_ALL, 3) +
        glbChiSys(fit_values, GLB_ALL, 4)
    @printf(stream, "The sum over all rules gives again: %g\n\n", chi2)

    # Prepare minimizers: Set errors for external parameters:
    # 10% for each of the solar parameters, 5% for the matter density
    glbDefineParams(input_errors, theta12 * 0.1, 0, 0, 0, sdm * 0.1, 0)
    glbSetDensityParams(input_errors, 0.05, GLB_ALL)
    glbSetCentralValues(true_values)
    glbSetInputErrors(input_errors)

    # Let's do the some with correlations included:
    chi2 = glbChiTheta13(fit_values, minimum, GLB_ALL)
    @printf(stream, "chi2 with correlations: %g \n", chi2)
    @printf(stream, "Position of minimum: theta12,theta13,theta23,delta,sdm,ldm,rho\n")
    glbPrintParams(stream, minimum)
    @printf(
        stream,
        "Note that s22theta13 is unchanged/kept fixed: %g! \n\n",
        sin(2 * glbGetOscParams(minimum, GLB_THETA_13))^2
    )

    # Check what the difference is if we keep in addition deltacp fixed.
    # This corresponds to projection onto theta13-deltacp-plane
    chi2 = glbChiTheta13Delta(fit_values, minimum, GLB_ALL)
    @printf(stream, "chi2 with correlations other than with deltacp: %g \n\n", chi2)

    # Similarly: what if we only take into account correlation with deltacp?
    # For this, we need an arbitrary user-defined projection:
    myprojection = glbAllocProjection()
    glbDefineProjection(
        myprojection,
        GLB_FIXED,
        GLB_FIXED,
        GLB_FIXED,
        GLB_FREE,
        GLB_FIXED,
        GLB_FIXED,
    )
    glbSetDensityProjectionFlag(myprojection, GLB_FREE, GLB_ALL)
    glbSetProjection(myprojection)
    chi2 = glbChiNP(fit_values, minimum, GLB_ALL)
    @printf(stream, "chi2 with correlation only with deltacp: %g \n\n", chi2)

    # Now: compare with statistics χ² only
    glbSwitchSystematics(GLB_ALL, GLB_ALL, GLB_OFF)
    chi2 = glbChiSys(fit_values, GLB_ALL, GLB_ALL)
    glbSwitchSystematics(GLB_ALL, GLB_ALL, GLB_ON)
    @printf(stream, "chi2 with statistics only: %g\n\n", chi2)

    # Find position of sgn-degeneracy
    glbDefineParams(input_errors, theta12 * 0.1, 0, 0, 0, sdm * 0.1, ldm / 3)
    glbSetDensityParams(input_errors, 0.05, GLB_ALL)
    glbDefineParams(central_values, theta12, theta13, theta23, deltacp, sdm, -ldm)
    glbSetDensityParams(central_values, 1.0, GLB_ALL)
    glbSetCentralValues(central_values)
    glbSetInputErrors(input_errors)
    chi2 = glbChiAll(central_values, minimum, GLB_ALL)
    @printf(stream, "chi2 at minimum: %g \n", chi2)
    @printf(stream, "Position of minimum: theta12,theta13,theta23,delta,sdm,ldm,rho\n")
    glbPrintParams(stream, minimum)

    # Now: multi-experiment setup
    @printf(stream, "\nNOW: TWO-EXPERIMENT SETUP NuFact@3000km+NuFact@7500km\n\n")

    # Initialize two experiments NFstandard.glb
    glbClearExperimentList()
    glbInitExperiment(joinpath(@__DIR__, "NFstandard.glb"))
    glbInitExperiment(joinpath(@__DIR__, "NFstandard.glb"))

    # Change baseline of second experiment to magic baseline
    lengths, densities = glbAverageDensityProfile(7500)
    @printf(stream, "Magic baseline length: %g, Density: %g\n\n", lengths[1], densities[1])
    glbSetProfileDataInExperiment(2, lengths, densities)

    # Initialize a number of parameter vector(s) again, because they depend
    # on the number of experiments
    true_values = glbAllocParams()
    fit_values = glbAllocParams()
    central_values = glbAllocParams()
    input_errors = glbAllocParams()
    minimum = glbAllocParams()

    # Assign standard oscillation parameters
    glbDefineParams(true_values, theta12, theta13, theta23, deltacp, sdm, ldm)
    glbSetDensityParams(true_values, 1.0, GLB_ALL)

    # The simulated data are computed with "true_values"
    glbSetOscillationParameters(true_values)
    glbSetRates()

    # Set test/fit values slightly off true values at s22th13=0.0015
    glbCopyParams(true_values, fit_values)
    glbSetOscParams(fit_values, asin(sqrt(0.0015)) / 2, GLB_THETA_13)

    # Compute χ² with systematics only for all experiments and rules
    chi2 = glbChiSys(fit_values, GLB_ALL, GLB_ALL)
    @printf(stream, "chi2 with systematics for all exps: %g\n", chi2)

    # Compute χ² with systematics only for each experiment
    chi2 = glbChiSys(fit_values, 1, GLB_ALL)
    @printf(stream, "chi2 with systematics for 3000km: %g\n", chi2)
    chi2b = glbChiSys(fit_values, 2, GLB_ALL)
    @printf(stream, "chi2 with systematics for 7500km: %g\n", chi2b)
    @printf(stream, "The two add again to: %g\n\n", chi2 + chi2b)

    # Compute χ² with correlations for each exp. and combination.
    # The sum of the two χ² is not equal to the χ² of the combination anymore,
    # since the minimum may have a different position and the priors are only added once.
    # Note that there are now two densities in the vectors!
    # Note that a minimum at a negative value of theta13 is unphysical. However, if there
    # can be none at a positive value found at small enough χ², there is none (Magic Baseline)
    glbDefineParams(input_errors, theta12 * 0.1, 0, 0, 0, sdm * 0.1, 0)
    glbSetDensityParams(input_errors, 0.05, GLB_ALL)
    glbSetCentralValues(true_values)
    glbSetInputErrors(input_errors)
    chi2 = glbChiTheta13(fit_values, minimum, 1)
    @printf(stream, "chi2 with correlations for 3000km: %g \n", chi2)
    glbPrintParams(stream, minimum)
    chi2b = glbChiTheta13(fit_values, minimum, 2)
    @printf(stream, "\nchi2 with correlations for 7500km: %g \n", chi2b)
    glbPrintParams(stream, minimum)
    chi2sum = glbChiTheta13(fit_values, minimum, GLB_ALL)
    @printf(stream, "\nchi2 with correlations for combination: %g \n", chi2sum)
    glbPrintParams(stream, minimum)
    @printf(
        stream,
        "\nThe sum of the two chi2s is %g, whereas the total chi2 is %g !\n\n",
        chi2 + chi2b,
        chi2sum
    )

    # Find sgn-degeneracy for 3000km experiments and test if it is still there for comb.
    glbDefineParams(input_errors, theta12 * 0.1, theta13, theta23, deltacp, sdm * 0.1, ldm / 3)
    glbSetDensityParams(input_errors, 0.05, GLB_ALL)
    glbDefineParams(central_values, theta12, theta13, theta23, deltacp, sdm, -ldm)
    glbSetDensityParams(central_values, 1.0, GLB_ALL)
    glbSetInputErrors(input_errors)
    glbSetCentralValues(central_values)
    chi2 = glbChiAll(central_values, minimum, 1)
    @printf(stream, "chi2 at minimum, L=3000km: %g \n", chi2)
    glbPrintParams(stream, minimum)
    chi2 = glbChiAll(minimum, minimum, GLB_ALL)
    @printf(stream, "\nchi2 for combination at minimum of L=3000km: %g \n", chi2)
    glbPrintParams(stream, minimum)

    stream === stdout || close(stream)
end

main()
