# Core data structures (port of glb_types.h)
#
# Integer "unset" markers follow the C code (-1) so that the defaulting
# logic of glbDefaultExp can be reproduced faithfully.

const GLB_OLD_NORM = -99

mutable struct Flux
    builtin::Int
    file_name::Union{Nothing,String}
    time::Float64
    parent_energy::Float64
    stored_muons::Float64
    target_power::Float64
    norm::Float64
    gamma::Float64
    end_point::Float64
    storage::Matrix{Float64}     # 501 x 7: E, nu_e, nu_mu, nu_tau, nubar_e, nubar_mu, nubar_tau
end
Flux() = Flux(-1, nothing, -1, -1, -1, -1, -1, -1, -1, zeros(0, 0))

mutable struct XSec
    builtin::Int
    file_name::Union{Nothing,String}
    storage::Matrix{Float64}     # 1001 x 7: log10(E), sigma/E for e,mu,tau, ebar,mubar,taubar
end
XSec() = XSec(-1, nothing, zeros(0, 0))

# Energy resolution function sigma(E, params)
const SigmaFunction = Function

mutable struct SmearData
    type::Int
    numofbins::Int
    simbins::Int
    e_min::Float64
    e_max::Float64
    e_sim_min::Float64
    e_sim_max::Float64
    sigma::Union{Nothing,Vector{Float64}}
    num_of_params::Int
    sig_f::Union{Nothing,Function}
    binsize::Union{Nothing,Vector{Float64}}
    simbinsize::Union{Nothing,Vector{Float64}}
    bincenter::Vector{Float64}
    simbincenter::Vector{Float64}
    corr_fac::Float64            # options->corr_fac (filter value for type 2)
end
SmearData() = SmearData(-1, -1, -1, -1, -1, -1, -1, nothing, -1, nothing,
                        nothing, nothing, Float64[], Float64[], -1)

"""
    ChiFunction(func, dim, name, user_data)

A registered χ² function. `func(exp, rule, n_params, x, errors, user_data)`
returns the χ² for the given experiment/rule and nuisance parameters `x`.
"""
mutable struct ChiFunction
    func::Function
    dim::Int
    name::String
    user_data::Any
end

mutable struct Channel
    flux::Int       # 1-based flux index
    cp::Int         # +1 neutrino, -1 antineutrino
    initial::Int    # 1,2,3 or 11,12,13 for NOSC_
    final::Int
    xsec::Int       # 1-based cross section index
    energy::Int     # 1-based smearing index
end

is_nosc(c::Channel) = c.initial > 9 || c.final > 9

mutable struct Rule
    signal_coeff::Vector{Float64}
    signal_channels::Vector{Int}       # 1-based
    bg_coeff::Vector{Float64}
    bg_channels::Vector{Int}
    sys_on_off::Int
    sys_on::Union{Nothing,ChiFunction}
    sys_off::Union{Nothing,ChiFunction}
    sys_on_string::Union{Nothing,String}
    sys_off_string::Union{Nothing,String}
    sys_on_errors::Union{Nothing,Vector{Float64}}
    sys_on_startvals::Union{Nothing,Vector{Float64}}
    sys_off_errors::Union{Nothing,Vector{Float64}}
    sys_off_startvals::Union{Nothing,Vector{Float64}}
    signal_errors::Vector{Float64}     # [norm, tilt] for GLoBES 2.0.11-style chi^2 functions
    signal_startvals::Vector{Float64}
    bg_errors::Vector{Float64}
    bg_startvals::Vector{Float64}
    bg_centers::Vector{Float64}        # [norm, tilt]
    energy_window::Vector{Float64}     # [low, high] in GeV
    energy_window_bins::Vector{Int}    # [low, high], 1-based inclusive
    errordim_sys_on::Int
    errordim_sys_off::Int
end
Rule() = Rule(Float64[], Int[], Float64[], Int[], GLB_ON, nothing, nothing, nothing, nothing,
              nothing, nothing, nothing, nothing, [-1.0, -1.0], [0.0, 0.0], [-1.0, -1.0],
              [0.0, 0.0], [1.0, 0.0], [-1.0, -1.0], [-1, -1], -1, -1)

"""
    Experiment

Everything GLoBES knows about one experiment: what the AEDL file defined,
the precomputed smearing matrices, and the rate buffers. Indices are 1-based.
"""
mutable struct Experiment
    version::Union{Nothing,String}
    filename::String
    names::Vector{Tuple{String,String,Int}}   # (name, context, 1-based value)

    fluxes::Vector{Flux}
    xsecs::Vector{XSec}

    binsize::Union{Nothing,Vector{Float64}}
    simbinsize::Union{Nothing,Vector{Float64}}

    density_profile_type::Int
    baseline::Float64
    emin::Float64
    emax::Float64
    numofbins::Int
    targetmass::Float64

    channels::Vector{Channel}
    rules::Vector{Rule}

    smear_data::Vector{SmearData}
    smear::Vector{Vector{Vector{Float64}}}  # [sm][bin] -> nonzero entries
    lowrange::Vector{Vector{Int}}           # [sm][bin] -> 1-based first sampling point
    uprange::Vector{Vector{Int}}
    user_smear_given::Vector{Bool}

    simtresh::Float64
    simbeam::Float64
    simbins::Int

    filter_state::Int
    filter_value::Float64

    psteps::Int
    lengthtab::Vector{Float64}
    densitytab::Vector{Float64}
    densitybuffer::Vector{Float64}

    user_pre_smearing_channel::Vector{Union{Nothing,Vector{Float64}}}
    user_post_smearing_channel::Vector{Union{Nothing,Vector{Float64}}}
    user_pre_smearing_background::Vector{Union{Nothing,Vector{Float64}}}
    user_post_smearing_background::Vector{Union{Nothing,Vector{Float64}}}

    energy_tab::Vector{Float64}

    # Rate buffers
    chrb_0::Vector{Vector{Float64}}   # true pre-smearing channel rates
    chrb_1::Vector{Vector{Float64}}   # fitted pre-smearing channel rates
    chra_0::Vector{Vector{Float64}}   # true post-smearing channel rates
    chra_1::Vector{Vector{Float64}}   # fitted post-smearing channel rates
    SignalRates::Vector{Vector{Float64}}
    BackgroundRates::Vector{Vector{Float64}}
    rates0::Vector{Vector{Float64}}
    rates1::Vector{Vector{Float64}}
    rates1BG::Vector{Vector{Float64}}
end

function Experiment()
    Experiment(nothing, "", Tuple{String,String,Int}[], Flux[], XSec[], nothing, nothing,
               -1, -1, -1, -1, -1, -1, Channel[], Rule[],
               SmearData[], Vector{Vector{Float64}}[], Vector{Int}[], Vector{Int}[], Bool[],
               -1, -1, -1, -1, -1, -1, Float64[], Float64[], Float64[],
               Union{Nothing,Vector{Float64}}[], Union{Nothing,Vector{Float64}}[],
               Union{Nothing,Vector{Float64}}[], Union{Nothing,Vector{Float64}}[],
               Float64[], Vector{Float64}[], Vector{Float64}[], Vector{Float64}[], Vector{Float64}[],
               Vector{Float64}[], Vector{Float64}[], Vector{Float64}[], Vector{Float64}[], Vector{Float64}[])
end

Base.show(io::IO, e::Experiment) =
    print(io, "Experiment(\"", basename(e.filename), "\", L=", e.baseline, " km, ",
          length(e.channels), " channels, ", length(e.rules), " rules, ", e.numofbins, " bins)")

# Oscillation / density parameter vector (glb_params)
"""
    GlbParams

Oscillation parameters (`osc`, length = number of oscillation parameters,
6 for the standard engine), matter density scaling factors (`density`, one per
experiment), and the iteration count of the last minimization.
"""
mutable struct GlbParams
    osc::Vector{Float64}
    density::Vector{Float64}
    iterations::Int
end

"""
    GlbProjection

Marks every oscillation parameter and density scaling as `GLB_FREE` or
`GLB_FIXED` for the χ² projections.
"""
mutable struct GlbProjection
    osc::Vector{Int}
    density::Vector{Int}
end
