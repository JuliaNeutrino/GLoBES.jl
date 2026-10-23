"""
    GLoBES

A Julia implementation of GLoBES, the General Long Baseline Experiment
Simulator (P. Huber, J. Kopp, M. Lindner, M. Rolinec, W. Winter).

It reads the same AEDL experiment description files and provides the GLoBES
C API under the same names (`glbInitExperiment`, `glbSetOscillationParameters`,
`glbSetRates`, `glbChiSys`, `glbChiTheta13`, ...). Differences to the C API:

* Experiment, rule and channel numbers are 1-based. The parameter index
  constants (`GLB_THETA_12`, ..., `GLB_DM_31`) are 1-based too, so code that
  uses the constants reads the same as in C.
* Functions that fill output pointers in C return values instead
  (e.g. `glbGetEminEmax(exp)` returns `(emin, emax)`).
* Rate "pointers" are returned as the live Julia vectors.

If you use this code for a publication, cite the GLoBES papers
(hep-ph/0407333 and hep-ph/0701187), as the original authors ask.
"""
module GLoBES

using LinearAlgebra
using Printf
using SpecialFunctions: erf
import Downloads
import p7zip_jll

const GLB_VERSION = "3.0.11"

# ---------------------------------------------------------------------------
#  Constants
# ---------------------------------------------------------------------------

const GLB_THETA_12 = 1
const GLB_THETA_13 = 2
const GLB_THETA_23 = 3
const GLB_DELTA_CP = 4
const GLB_DM_21 = 5
const GLB_DM_31 = 6
const GLB_DM_SOL = GLB_DM_21
const GLB_DM_ATM = GLB_DM_31

const GLB_OFF = 0
const GLB_ON = 1

const GLB_W_EFF = 1
const GLB_WO_EFF = 2
const GLB_W_BG = 3
const GLB_WO_BG = 4
const GLB_W_COEFF = 5
const GLB_WO_COEFF = 6

const GLB_FIXED = 0
const GLB_FREE = 1
const GLB_UNDEFINED = -2

const GLB_PRE = 1
const GLB_POST = 2

const GLB_MIN_NESTED_POWELL = 0
const GLB_MIN_POWELL = 1
const GLB_MIN_DEFAULT = GLB_MIN_NESTED_POWELL

const GLB_EFF = 1
const GLB_BG = 2
const GLB_SIG = 3

const GLB_ALL = -1

const GLB_EARTH_RADIUS = 6371.0
const GLB_EV_TO_KM_FACTOR = 1.973269804e-10   # ħc in eV km (C GLoBES: 1.9747235e-10)
GLB_KM_TO_EV(x) = x / GLB_EV_TO_KM_FACTOR
GLB_EV_TO_KM(x) = x * GLB_EV_TO_KM_FACTOR

const GLB_NAN = NaN

# ---------------------------------------------------------------------------
#  Error handling
# ---------------------------------------------------------------------------

"""
    GLoBESError

An error of a GLoBES.jl function: an impossible request or a missing experiment.

```julia
simulated_params(; setup = Setup())   # throws a GLoBESError: nothing is simulated in a new setup
```
"""
struct GLoBESError <: Exception
    msg::String
end
Base.showerror(io::IO, e::GLoBESError) = print(io, "GLoBES error: ", e.msg)

const verbosity = Ref(1)

"""
    glbSetVerbosityLevel(level)

0 = silent, 1 = errors (default), 2 = errors and warnings.
"""
glbSetVerbosityLevel(level::Integer) = (verbosity[] = level; 0)

# The glb* functions report errors like the C library: print and carry on with
# a default or an error code. The high-level API runs them in `strict` mode,
# where errors throw a GLoBESError instead.
const STRICT = Ref(0)

function glb_error(msg)
    STRICT[] > 0 && throw(GLoBESError(msg))
    verbosity[] >= 1 && println(stderr, "GLoBES: ERROR: ", msg)
    return nothing
end

function strict(f)
    STRICT[] += 1
    try
        return f()
    finally
        STRICT[] -= 1
    end
end

glb_warning(msg) = (verbosity[] >= 2 && println(stderr, "GLoBES: Warning: ", msg); nothing)
glb_fatal(msg) = throw(GLoBESError(msg))

# ---------------------------------------------------------------------------
#  Data files: current directory, then the directories in GLB_PATH
# ---------------------------------------------------------------------------

const glb_path_vector = String[]

function glb_setup_path!()
    empty!(glb_path_vector)
    push!(glb_path_vector, pwd())
    p = get(ENV, "GLB_PATH", "")
    for d in split(p, ':')
        isempty(d) || push!(glb_path_vector, d)
    end
end

"""
    glb_find_file(name, extra_dir=nothing) -> path or nothing

Locate a data file like GLoBES does: absolute path, then the current
directory and the `GLB_PATH` directories. `extra_dir` (the directory of the
including AEDL file) is tried after the current directory.
"""
function glb_find_file(name::AbstractString, extra_dir = nothing)
    isabspath(name) && return isfile(name) ? String(name) : nothing
    dirs = String[pwd()]
    extra_dir !== nothing && push!(dirs, extra_dir)
    append!(dirs, glb_path_vector)
    for d in dirs
        f = joinpath(d, name)
        isfile(f) && return f
    end
    return nothing
end

include("types.jl")
include("probability.jl")
include("profile.jl")
include("fluxes.jl")
include("smear.jl")

# ---------------------------------------------------------------------------
#  Global state (the C library keeps all of this in static variables)
# ---------------------------------------------------------------------------

builtin_prior(p, ud) = builtin_prior_impl(p)
noop_hook(p, ud) = 0

mutable struct GlobalState
    experiments::Vector{ExperimentData}
    oscp::Int
    engine::StandardEngine
    prob_func::Function
    set_func::Function
    get_func::Function
    prob_user_data::Any
    P::Matrix{Float64}
    Pbar::Matrix{Float64}
    sys_list::Vector{ChiFunction}
    minimizer::Int
    prior_func::Function
    starting_values_hook::Function
    input_errors_hook::Function
    prior_user_data::Any
    input_errors::GlbParams
    central_values::GlbParams
    para_tab::Vector{Int}
    index_tab::Vector{Int}
    n_free::Int
    n_fix::Int
    fix_params::Vector{Float64}
    count::Int
    single_experiment::Int
    single_rule::Int
    last_priors::Float64
    initialized::Bool
    truth::Union{Nothing, GlbParams}   # parameters of the simulated data (high-level API)
    extra_names::Vector{Symbol}       # names of non-standard engine parameters (high-level API)
end

new_state() = GlobalState(
    ExperimentData[],
    6,
    StandardEngine(),
    (a...) -> 0,
    (a...) -> 0,
    (a...) -> nothing,
    nothing,
    zeros(3, 3),
    zeros(3, 3),
    ChiFunction[],
    GLB_MIN_DEFAULT,
    builtin_prior,
    noop_hook,
    noop_hook,
    nothing,
    GlbParams(fill(NaN, 6), Float64[], 0),
    GlbParams(fill(NaN, 6), Float64[], 0),
    Int[],
    Int[],
    0,
    0,
    Float64[],
    0,
    GLB_ALL,
    GLB_ALL,
    0.0,
    false,
    nothing,
    Symbol[],
)

"""
    Setup()

An independent GLoBES simulation: its own experiments, simulated data,
probability engine, χ² functions, priors and minimizer settings.

Most scripts never need more than the default setup, which every function
uses unless told otherwise. Use [`with_setup`](@ref) to work with another.
"""
mutable struct Setup
    state::GlobalState
end

const CURRENT_SETUP = Ref{Setup}(Setup(new_state()))

# The engine reads and changes the state of the current setup
@inline state() = CURRENT_SETUP[].state

include("experiment.jl")
include("aedl.jl")
include("systematics.jl")
include("minimize.jl")
include("api.jl")
include("highlevel.jl")
include("mathematica.jl")
include("sarah.jl")
include("spectrum.jl")
include("mechanisms.jl")
include("newphysics.jl")
include("nufit.jl")
include("oscmodel.jl")
include("latex.jl")
include("analysis.jl")
include("runfile.jl")

function __init__()
    glbInit()
end

# ---------------------------------------------------------------------------
#  Exports
# ---------------------------------------------------------------------------

export GLB_THETA_12,
    GLB_THETA_13,
    GLB_THETA_23,
    GLB_DELTA_CP,
    GLB_DM_21,
    GLB_DM_31,
    GLB_DM_SOL,
    GLB_DM_ATM,
    GLB_OFF,
    GLB_ON,
    GLB_W_EFF,
    GLB_WO_EFF,
    GLB_W_BG,
    GLB_WO_BG,
    GLB_W_COEFF,
    GLB_WO_COEFF,
    GLB_FIXED,
    GLB_FREE,
    GLB_PRE,
    GLB_POST,
    GLB_MIN_NESTED_POWELL,
    GLB_MIN_POWELL,
    GLB_MIN_DEFAULT,
    GLB_EFF,
    GLB_BG,
    GLB_SIG,
    GLB_ALL,
    GLB_EARTH_RADIUS,
    GLB_EV_TO_KM_FACTOR,
    GLB_KM_TO_EV,
    GLB_EV_TO_KM,
    GLB_NAN

export GlbParams, GlbProjection, ExperimentData, GLoBESError, AEDLError

# High-level API
export Setup,
    current_setup,
    with_setup,
    OscParams,
    with_param,
    Experiment,
    load_experiments,
    experiments,
    clear_experiments!,
    name,
    baseline,
    baseline!,
    target_mass,
    target_mass!,
    running_time,
    running_time!,
    source_power,
    source_power!,
    energy_range,
    bin_centers,
    bin_widths,
    density_profile,
    density_profile!,
    prem_profile,
    average_profile,
    rule_names,
    channel_names,
    energy_window,
    energy_window!,
    energy_window_bins,
    simulate!,
    simulated_params,
    event_rates,
    total_events,
    channel_rates,
    fit_rates,
    true_rates,
    probability,
    chi2,
    marginalize,
    FitResult,
    minimizer!,
    chi2_scan,
    systematics!,
    chi2_function,
    chi2_function!,
    systematic_errors,
    systematic_errors!,
    define_chi2,
    poisson_chi2,
    energy_scale_shift,
    probability_engine!,
    standard_engine!,
    flux,
    cross_section

# Models, NuFIT, reports
export BSMModel,
    load_sarah,
    check_model,
    add_gauge!,
    add_global!,
    add_fermion!,
    add_scalar!,
    add_lagrangian!,
    add_vevs!,
    add_parameter!,
    spectrum,
    read_leshouches,
    leshouches_values,
    write_leshouches,
    ScotogenicLoop,
    TypeISeesaw,
    neutrino_mass_matrix,
    neutrino_spectrum,
    casas_ibarra,
    orthogonal_matrix,
    pmns_matrix,
    light_masses,
    br_lfv,
    scotogenic_masses,
    NewPhysics,
    ScalarNSI,
    MatterNSI,
    Decoherence,
    use_model!,
    nufit,
    nufit_params,
    nufit_priors,
    nufit_ranges,
    nufit_chi2,
    nufit_prior,
    nufit_download,
    OscillationModel,
    load_model,
    model_from_string,
    model_parameters,
    hamiltonian,
    effective_mixing,
    effective_masses,
    effective_parameters,
    model_probability,
    run_model,
    ModelRun,
    latex,
    latex_expression,
    latex_table,
    TexReport,
    add!,
    section!,
    figure!,
    table!,
    write_tex,
    compile_tex,
    probability_band,
    mass_ordering_sensitivity,
    cp_violation_sensitivity,
    chi2_profile,
    plot_probability,
    plot_events,
    plot_curves,
    plot_contour,
    plot_nufit,
    plot_grid,
    paper_theme!

export glbInit,
    glbInitExperiment,
    glbClearExperimentList,
    glbNumOfExps,
    glbGetExperiment,
    glbSetVerbosityLevel,
    glbVersionOfExperiment,
    glbAllocParams,
    glbFreeParams,
    glbDefineParams,
    glbCopyParams,
    glbPrintParams,
    glbSetOscParams,
    glbGetOscParams,
    glbSetDensityParams,
    glbGetDensityParams,
    glbSetIteration,
    glbGetIteration,
    glbAllocProjection,
    glbFreeProjection,
    glbDefineProjection,
    glbCopyProjection,
    glbPrintProjection,
    glbSetProjectionFlag,
    glbGetProjectionFlag,
    glbSetDensityProjectionFlag,
    glbGetDensityProjectionFlag,
    glbSetProjection,
    glbGetProjection,
    glbSetOscillationParameters,
    glbGetOscillationParameters,
    glbSetRates,
    glbSetCentralValues,
    glbGetCentralValues,
    glbSetStartingValues,
    glbGetStartingValues,
    glbSetInputErrors,
    glbGetInputErrors,
    glbChiSys,
    glbChiNP,
    glbChiTheta12,
    glbChiTheta13,
    glbChiTheta23,
    glbChiDelta,
    glbChiDm21,
    glbChiDm31,
    glbChiTheta13Delta,
    glbChiAll,
    glbSelectMinimizer,
    glbRegisterPriorFunction,
    glbRegisterProbabilityEngine,
    glbGetNumOfOscParams,
    glbDefineChiFunction,
    glbSetChiFunction,
    glbGetChiFunction,
    glbGetChiFunctionPtr,
    glbGetChiFunctionPtrInExperiment,
    glbSetErrorDim,
    glbGetErrorDim,
    glbGetSysDim,
    glbGetSysDimInExperiment,
    glbSwitchSystematics,
    glbGetSysOnOffState,
    glbSetSignalErrors,
    glbGetSignalErrors,
    glbSetBGErrors,
    glbGetBGErrors,
    glbSetBGCenters,
    glbGetBGCenters,
    glbSetSignalStartingValues,
    glbGetSignalStartingValues,
    glbSetBGStartingValues,
    glbGetBGStartingValues,
    glbSetSysErrorsList,
    glbGetSysErrorsListPtr,
    glbSetSysStartingValuesList,
    glbGetSysStartingValuesListPtr,
    glbShiftEnergyScale,
    glbChiSpectrumTilt,
    glbChiNoSysSpectrum,
    glbChiSpectrumOnly,
    glbChiTotalRatesTilt,
    glbChiNoSysTotalRates,
    glbChiSpectrumCalib,
    glbChiZero,
    glb_likelihood,
    glb_prior,
    glbGetRuleRatePtr,
    glbGetSignalRatePtr,
    glbGetBGRatePtr,
    glbGetSignalFitRatePtr,
    glbGetBGFitRatePtr,
    glbGetChannelRatePtr,
    glbGetChannelFitRatePtr,
    glbGetEfficiencyPtr,
    glbGetBackgroundPtr,
    glbGetChannelRates,
    glbGetUserData,
    glbTotalRuleRate,
    glbShowRuleRates,
    glbShowChannelRates,
    glbShowChannelProbs,
    glbGetNumberOfBins,
    glbGetNumberOfSamplingPoints,
    glbGetNumberOfChannels,
    glbGetNumberOfRules,
    glbGetNumberOfFluxes,
    glbGetBinSizeListPtr,
    glbGetBinCentersListPtr,
    glbGetSamplingStepsizeListPtr,
    glbGetSamplingPointsListPtr,
    glbGetEminEmax,
    glbGetEnergyWindow,
    glbSetEnergyWindow,
    glbGetEnergyWindowBins,
    glbGetLengthOfRule,
    glbGetChannelInRule,
    glbGetCoefficientInRule,
    glbGetNormalizationInRule,
    glbNameToValue,
    glbValueToName,
    glbSetTargetMass,
    glbGetTargetMass,
    glbSetSourcePower,
    glbGetSourcePower,
    glbSetRunningTime,
    glbGetRunningTime,
    glbSetFilterState,
    glbGetFilterState,
    glbSetFilterStateInExperiment,
    glbGetFilterStateInExperiment,
    glbSetFilter,
    glbGetFilter,
    glbSetFilterInExperiment,
    glbGetFilterInExperiment,
    glbFlux,
    glbXSection,
    glbLoadProfileData,
    glbStaceyProfile,
    glbAverageDensityProfile,
    glbGetProfileData,
    glbGetProfileDataInExperiment,
    glbSetProfileDataInExperiment,
    glbSetBaselineInExperiment,
    glbGetBaselineInExperiment,
    glbGetProfileTypeInExperiment,
    glbVacuumProbability,
    glbConstantDensityProbability,
    glbProfileProbability,
    glbFilteredConstantDensityProbability,
    glbDefineAEDLVariable,
    glbDefineAEDLList,
    glbClearAEDLVariables

end # module
