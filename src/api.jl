# User-facing API (glb_wrapper.c, globes.c, glb_glb.c, glb_prior.c)

# ---------------------------------------------------------------------------
#  Initialization and experiments
# ---------------------------------------------------------------------------

"""
    glbInit()

(Re)initialize GLoBES: remove all experiments, install the standard
probability engine, the default prior and the built-in χ² functions.
Called automatically when the package is loaded.
"""
function glbInit(name::AbstractString = "GLoBES.jl")
    st = state()
    glb_setup_path!()
    empty!(st.experiments)
    glbRegisterProbabilityEngine(6, nothing, nothing, nothing, nothing)
    glbSelectMinimizer(GLB_MIN_DEFAULT)
    glbRegisterPriorFunction(nothing, nothing, nothing, nothing)
    register_builtin_chi_functions!()
    init_minimizer!()
    st.initialized = true
    return nothing
end

function init_minimizer!()
    st = state()
    st.input_errors = glbAllocParams()
    st.central_values = glbAllocParams()
    reset_projection!()
end

"""
    glbInitExperiment(filename) -> status

Read an AEDL file and append its experiment(s) to the experiment list.
Returns 0 on success, -1 if defaults had to be used, -2 on failure.
The file is searched in the current directory and in `GLB_PATH`.

Like the C library, this resets projection, central values and input errors.
"""
function glbInitExperiment(filename::AbstractString)
    path = glb_find_file(filename)
    if path === nothing
        glb_error("Could not open file $filename")
        return -2
    end
    exps = try
        parse_aedl_file(path)
    catch err
        (err isa AEDLError && STRICT[] == 0) || rethrow()
        glb_error(sprint(showerror, err))
        return -2
    end
    length(state().experiments) + length(exps) > 32 && glb_fatal("Too many experiments!")
    k = 0
    for e in exps
        k += default_exp!(e)
        push!(state().experiments, e)
    end
    init_minimizer!()
    return k != 0 ? -1 : 0
end

"""
    glbClearExperimentList()

Remove all experiments.
"""
function glbClearExperimentList()
    empty!(state().experiments)
    init_minimizer!()
    return nothing
end

"""
    glbNumOfExps() -> Int

Number of loaded experiments (`glb_num_of_exps` in C).
"""
glbNumOfExps() = length(state().experiments)

"""
    glbGetExperiment(exp) -> ExperimentData

Direct access to the internal experiment structure.
"""
glbGetExperiment(exp::Integer) = state().experiments[exp]

function getexp(exp)
    (1 <= exp <= length(state().experiments)) ||
        (glb_error("Invalid value for experiment number"); return nothing)
    return state().experiments[exp]
end

function getrule(exp, rule)
    e = getexp(exp)
    e === nothing && return nothing
    (1 <= rule <= length(e.rules)) || (glb_error("Invalid rule number"); return nothing)
    return e.rules[rule]
end

glbVersionOfExperiment(exp) = (e = getexp(exp); e === nothing ? nothing : e.version)

# ---------------------------------------------------------------------------
#  Parameters and projections
# ---------------------------------------------------------------------------

"""
    glbAllocParams() -> GlbParams

New parameter vector with all entries `NaN` (sized for the current
probability engine and number of experiments).
"""
glbAllocParams() = GlbParams(fill(NaN, state().oscp), fill(NaN, length(state().experiments)), 0)
glbFreeParams(p) = nothing

"""
    glbDefineParams(p, θ12, θ13, θ23, δ, Δm²21, Δm²31)
"""
function glbDefineParams(p::GlbParams, th12, th13, th23, delta, dms, dma)
    p.osc[1:6] .= (th12, th13, th23, delta, dms, dma)
    return p
end

function glbSetOscParams(p::GlbParams, v, which)
    (1 <= which <= state().oscp) ||
        (glb_error("glbSetOscParams: Oscillation list length mismatch"); return nothing)
    p.osc[which] = v
    return p
end

function glbGetOscParams(p::GlbParams, which)
    (1 <= which <= length(p.osc)) ||
        (glb_error("glbGetOscParams: Oscillation list length mismatch"); return 0.0)
    return p.osc[which]
end

function glbSetDensityParams(p::GlbParams, v, which)
    if which == GLB_ALL
        fill!(p.density, v)
    elseif 1 <= which <= length(p.density)
        p.density[which] = v
    else
        glb_error("glbSetDensityParams: Density list length mismatch")
        return nothing
    end
    return p
end

function glbGetDensityParams(p::GlbParams, which)
    (1 <= which <= length(p.density)) ||
        (glb_error("glbGetDensityParams: Density list length mismatch"); return 0.0)
    return p.density[which]
end

glbSetIteration(p::GlbParams, it) = (p.iterations = it; p)
glbGetIteration(p::GlbParams) = p.iterations

function glbCopyParams(src::GlbParams, dest::GlbParams)
    src === dest && return dest
    dest.iterations = src.iterations
    dest.osc = copy(src.osc)
    dest.density = copy(src.density)
    return dest
end
glbCopyParams(src::GlbParams) = glbCopyParams(src, glbAllocParams())

function glbPrintParams(io::IO, p::GlbParams)
    println(io, join((@sprintf("%g", x) for x in p.osc), " "), " ")
    println(io, join((@sprintf("%g", x) for x in p.density), " "), " ")
    println(io, "Iterations: ", p.iterations)
end
glbPrintParams(p::GlbParams) = glbPrintParams(stdout, p)

Base.show(io::IO, p::GlbParams) =
    print(io, "GlbParams(osc=", p.osc, ", density=", p.density, ", iterations=", p.iterations, ")")

glbAllocProjection() = GlbProjection(
    fill(GLB_UNDEFINED, state().oscp),
    fill(GLB_UNDEFINED, length(state().experiments)),
)
glbFreeProjection(p) = nothing

"""
    glbDefineProjection(proj, θ12, θ13, θ23, δ, Δm²21, Δm²31)

Each argument is `GLB_FREE` or `GLB_FIXED`.
"""
function glbDefineProjection(p::GlbProjection, a, b, c, d, e, f)
    p.osc[1:6] .= (a, b, c, d, e, f)
    return p
end

function glbSetProjectionFlag(p::GlbProjection, flag, which)
    (1 <= which <= length(p.osc)) ||
        (glb_error("glbSetProjectionFlag: list length mismatch"); return nothing)
    p.osc[which] = flag
    return p
end
glbGetProjectionFlag(p::GlbProjection, which) = p.osc[which]

function glbSetDensityProjectionFlag(p::GlbProjection, flag, which)
    if which == GLB_ALL
        fill!(p.density, flag)
    elseif 1 <= which <= length(p.density)
        p.density[which] = flag
    else
        glb_error("glbSetDensityProjectionFlag: Density list length mismatch")
        return nothing
    end
    return p
end
glbGetDensityProjectionFlag(p::GlbProjection, which) = p.density[which]

function glbCopyProjection(src::GlbProjection, dest::GlbProjection)
    dest.osc = copy(src.osc)
    dest.density = copy(src.density)
    return dest
end

function glbPrintProjection(io::IO, p::GlbProjection)
    fl(x) = x == GLB_FREE ? "free" : x == GLB_FIXED ? "fixed" : "undefined"
    names = ("Theta12", "Theta13", "Theta23", "deltacp", "dm21", "dm31")
    for (i, x) in enumerate(p.osc)
        println(io, i <= 6 ? names[i] : "Param $i", ": ", fl(x))
    end
    for (i, x) in enumerate(p.density)
        println(io, "Density $i: ", fl(x))
    end
end
glbPrintProjection(p::GlbProjection) = glbPrintProjection(stdout, p)

Base.show(io::IO, p::GlbProjection) =
    print(io, "GlbProjection(osc=", p.osc, ", density=", p.density, ")")

# ---------------------------------------------------------------------------
#  Oscillation engine
# ---------------------------------------------------------------------------

"""
    glbSetOscillationParameters(p)

Set the oscillation parameters used by the probability engine.
"""
function glbSetOscillationParameters(p::GlbParams)
    any(isnan, p.osc) &&
        (glb_error("glbSetOscillationParameters: parameters incompletely defined"); return -1)
    return state().set_func(p, state().prob_user_data)
end

"""
    glbGetOscillationParameters(p)

Copy the current oscillation parameters into `p`.
"""
function glbGetOscillationParameters(p::GlbParams)
    state().get_func(p, state().prob_user_data)
    return 0
end

glbGetNumOfOscParams() = state().oscp

default_prob_func(P, cp, E, ps, l, d, fs, ud) =
    probability_matrix!(P, state().engine, cp, E, ps, l, d, fs)
default_set_func(p, ud) = set_oscillation_parameters!(state().engine, p.osc)
function default_get_func(p, ud)
    v = get_oscillation_parameters(state().engine)
    p.osc[1:6] .= v
    return 0
end

"""
    glbRegisterProbabilityEngine(n_params, prob_func, set_func, get_func, user_data)

Replace the oscillation engine. `prob_func(P, cp_sign, E, psteps, lengths,
densities, filter_sigma, user_data)` must fill the 3×3 matrix `P[α,β]` and
return 0; `set_func(p::GlbParams, user_data)` and `get_func(p, user_data)`
set/get the engine's parameters. `nothing` restores the default for that hook.
"""
function glbRegisterProbabilityEngine(n_params, prob_func, set_func, get_func, user_data = nothing)
    st = state()
    st.oscp = (n_params === nothing || n_params <= 0) ? 6 : n_params
    st.prob_func = prob_func === nothing ? default_prob_func : prob_func
    st.set_func = set_func === nothing ? default_set_func : set_func
    st.get_func = get_func === nothing ? default_get_func : get_func
    st.prob_user_data = user_data
    init_minimizer!()
    return 0
end

# ---------------------------------------------------------------------------
#  Priors
# ---------------------------------------------------------------------------

sprior(x, center, sigma) = abs(sigma) < 1e-12 ? 0.0 : (x - center)^2 / sigma^2

# The default prior: Gaussian terms for all free parameters, using the
# central values and input errors set with glbSetCentralValues/glbSetInputErrors
function builtin_prior_impl(p::GlbParams)
    st = state()
    res = 0.0
    N = st.oscp
    length(st.para_tab) < N + length(st.experiments) && return 0.0
    for i = 1:N
        st.para_tab[i] == GLB_FREE &&
            (res += sprior(p.osc[i], st.central_values.osc[i], st.input_errors.osc[i]))
    end
    for i in eachindex(st.experiments)
        st.para_tab[N + i] == GLB_FREE &&
            (res += sprior(p.density[i], st.central_values.density[i], st.input_errors.density[i]))
    end
    return res
end

"""
    glbRegisterPriorFunction(prior, starting, error, user_data)

Install a user-defined prior `prior(p::GlbParams, user_data) -> Float64`.
`starting` and `error` are called by `glbSetCentralValues`/`glbSetInputErrors`
with the new values. `nothing` for all three restores the built-in prior.
"""
function glbRegisterPriorFunction(prior, starting, error, user_data = nothing)
    st = state()
    st.prior_func = prior === nothing ? builtin_prior : prior
    st.starting_values_hook = starting === nothing ? noop_hook : starting
    st.input_errors_hook = error === nothing ? noop_hook : error
    st.prior_user_data = user_data
    return 0
end

"""
    glbSetCentralValues(p)

Central values of the Gaussian priors for free parameters.
"""
function glbSetCentralValues(p::GlbParams)
    glbCopyParams(p, state().central_values)
    return state().starting_values_hook(state().central_values, state().prior_user_data)
end
glbGetCentralValues(p::GlbParams) = (glbCopyParams(state().central_values, p); 0)
const glbSetStartingValues = glbSetCentralValues
const glbGetStartingValues = glbGetCentralValues

"""
    glbSetInputErrors(p)

Widths of the Gaussian priors (0 = no prior) for free parameters.
"""
function glbSetInputErrors(p::GlbParams)
    glbCopyParams(p, state().input_errors)
    return state().input_errors_hook(state().input_errors, state().prior_user_data)
end
glbGetInputErrors(p::GlbParams) = (glbCopyParams(state().input_errors, p); 0)

# ---------------------------------------------------------------------------
#  Event rates
# ---------------------------------------------------------------------------

"""
    glbSetRates()

Compute the "true" event rates of all experiments for the current
oscillation parameters.
"""
function glbSetRates()
    for e in state().experiments
        set_profile_scaling!(e, 1.0)
        set_rates!(e)
    end
    return nothing
end

"""
    glbSetNewRates()

Compute the fitted event rates of all experiments for the current
oscillation parameters.
"""
function glbSetNewRates()
    for e in state().experiments
        set_profile_scaling!(e, 1.0)
        set_new_rates!(e)
    end
    return nothing
end

glbGetRuleRatePtr(exp, rule) =
    (getrule(exp, rule) === nothing ? nothing : state().experiments[exp].rates0[rule])
glbGetSignalRatePtr(exp, rule) =
    (getrule(exp, rule) === nothing ? nothing : state().experiments[exp].SignalRates[rule])
glbGetBGRatePtr(exp, rule) =
    (getrule(exp, rule) === nothing ? nothing : state().experiments[exp].BackgroundRates[rule])
glbGetSignalFitRatePtr(exp, rule) =
    (getrule(exp, rule) === nothing ? nothing : state().experiments[exp].rates1[rule])
glbGetBGFitRatePtr(exp, rule) =
    (getrule(exp, rule) === nothing ? nothing : state().experiments[exp].rates1BG[rule])

function getchannel_vec(exp, channel, a, b, pre_post)
    e = getexp(exp)
    e === nothing && return nothing
    (1 <= channel <= length(e.channels)) || (glb_error("Invalid channel number"); return nothing)
    pre_post == GLB_PRE && return getfield(e, a)[channel]
    pre_post == GLB_POST && return getfield(e, b)[channel]
    glb_error("pre_post must be GLB_PRE or GLB_POST")
    return nothing
end

glbGetChannelRatePtr(exp, channel, pre_post) =
    getchannel_vec(exp, channel, :chrb_0, :chra_0, pre_post)
glbGetChannelFitRatePtr(exp, channel, pre_post) =
    getchannel_vec(exp, channel, :chrb_1, :chra_1, pre_post)
glbGetEfficiencyPtr(exp, channel, pre_post) =
    getchannel_vec(exp, channel, :user_pre_smearing_channel, :user_post_smearing_channel, pre_post)
glbGetBackgroundPtr(exp, channel, pre_post) = getchannel_vec(
    exp,
    channel,
    :user_pre_smearing_background,
    :user_post_smearing_background,
    pre_post,
)

"""
    glbGetChannelRates(exp, channel, smearing) -> Vector

Copy of the true rates of a channel, `GLB_PRE` or `GLB_POST` smearing.
"""
glbGetChannelRates(exp, channel, smearing) =
    (v = glbGetChannelRatePtr(exp, channel, smearing); v === nothing ? nothing : copy(v))

"""
    glbGetUserData(exp, channel, smearing, bgeff) -> Vector

Copy of the efficiencies (`GLB_EFF`) or backgrounds (`GLB_BG`) of a channel.
"""
function glbGetUserData(exp, channel, smearing, bgeff)
    v =
        bgeff == GLB_EFF ? glbGetEfficiencyPtr(exp, channel, smearing) :
        bgeff == GLB_BG ? glbGetBackgroundPtr(exp, channel, smearing) :
        (glb_error("bgeff must be GLB_EFF or GLB_BG"); nothing)
    return v === nothing ? nothing : copy(v)
end

glbGetLengthOfRule(exp, rule, signal) = (
    r = getrule(exp, rule);
    r === nothing ? -1 : signal == GLB_SIG ? length(r.signal_channels) : length(r.bg_channels)
)
glbGetChannelInRule(exp, rule, pos, signal) = (
    r = getrule(exp, rule);
    r === nothing ? -1 : signal == GLB_SIG ? r.signal_channels[pos] : r.bg_channels[pos]
)
glbGetCoefficientInRule(exp, rule, pos, signal) = (
    r = getrule(exp, rule);
    r === nothing ? -1.0 : signal == GLB_SIG ? r.signal_coeff[pos] : r.bg_coeff[pos]
)
glbGetNormalizationInRule(exp, rule, signal) = (
    r = getrule(exp, rule);
    r === nothing ? -1.0 : signal == GLB_SIG ? 1.0 : signal == GLB_BG ? r.bg_centers[1] : -1.0
)

function rule_component_rates(exp, rule, pos, effi, bgi, coeffi, signal)
    n = glbGetLengthOfRule(exp, rule, signal)
    positions = pos == GLB_ALL ? (1:n) : (pos:pos)
    res = Vector{Vector{Float64}}()
    for i in positions
        ch = glbGetChannelInRule(exp, rule, i, signal)
        coeff =
            glbGetCoefficientInRule(exp, rule, i, signal) *
            glbGetNormalizationInRule(exp, rule, signal)
        rates = copy(glbGetChannelRatePtr(exp, ch, GLB_POST))
        bgi == GLB_WO_BG && (rates .-= glbGetBackgroundPtr(exp, ch, GLB_POST))
        effi == GLB_WO_EFF && (rates ./= glbGetEfficiencyPtr(exp, ch, GLB_POST))
        coeffi == GLB_W_COEFF && (rates .*= coeff)
        push!(res, rates)
    end
    return res
end

"""
    glbTotalRuleRate(exp, rule, pos, effi, bgi, coeffi, signal) -> Float64

Total (all bins) signal (`GLB_SIG`) or background (`GLB_BG`) events of a rule,
optionally restricted to the channel at position `pos` of the rule.
`effi = GLB_W_EFF/GLB_WO_EFF`, `bgi = GLB_W_BG/GLB_WO_BG`,
`coeffi = GLB_W_COEFF/GLB_WO_COEFF`.
"""
function glbTotalRuleRate(exp, rule, pos, effi, bgi, coeffi, signal)
    getrule(exp, rule) === nothing && return -1.0
    return sum(sum, rule_component_rates(exp, rule, pos, effi, bgi, coeffi, signal); init = 0.0)
end

function print_channel_table(io, energy, res)
    println(io)
    for i in eachindex(energy)
        print(io, @sprintf("%6.4g", energy[i]), "\t")
        for r in res
            print(io, @sprintf("%12.6g", r[i]), "\t")
        end
        println(io)
    end
    print(io, "----------------------", "----------------"^max(length(res) - 1, 0))
    print(io, "\nTotal:\t")
    for r in res
        print(io, @sprintf("%12.6g", sum(r)), "\t")
    end
    println(io)
end

"""
    glbShowRuleRates(io, exp, rule, pos, effi, bgi, coeffi, signal)

Print the binned rates of a rule (all components for `pos = GLB_ALL`).
"""
function glbShowRuleRates(io::IO, exp, rule, pos, effi, bgi, coeffi, signal)
    getrule(exp, rule) === nothing && return -1
    res = rule_component_rates(exp, rule, pos, effi, bgi, coeffi, signal)
    print_channel_table(io, state().experiments[exp].smear_data[1].bincenter, res)
    return 0
end
glbShowRuleRates(exp, rule, pos, effi, bgi, coeffi, signal) =
    glbShowRuleRates(stdout, exp, rule, pos, effi, bgi, coeffi, signal)

"""
    glbShowChannelRates(io, exp, channel, smearing, effi, bgi)

Print the binned rates of one channel or of all channels (`GLB_ALL`).
"""
function glbShowChannelRates(io::IO, exp, channel, smearing, effi, bgi)
    e = getexp(exp)
    e === nothing && return -1
    chans = channel == GLB_ALL ? (1:length(e.channels)) : (channel:channel)
    res = Vector{Vector{Float64}}()
    for c in chans
        r = glbGetChannelRates(exp, c, smearing)
        bgi == GLB_WO_BG && (r .-= glbGetBackgroundPtr(exp, c, smearing))
        effi == GLB_WO_EFF && (r ./= glbGetEfficiencyPtr(exp, c, smearing))
        push!(res, r)
    end
    sd = e.smear_data[e.channels[first(chans)].energy]
    energy = smearing == GLB_PRE ? sd.simbincenter : sd.bincenter
    print_channel_table(io, energy, res)
    return 0
end
glbShowChannelRates(exp, channel, smearing, effi, bgi) =
    glbShowChannelRates(stdout, exp, channel, smearing, effi, bgi)

"""
    glbShowChannelProbs(io, exp, channel, smearing, effi, bgi)

Print the oscillation probability of the channel(s) at the bin energies.
"""
function glbShowChannelProbs(io::IO, exp, channel, smearing, effi = GLB_W_EFF, bgi = GLB_W_BG)
    e = getexp(exp)
    e === nothing && return -1
    chans = channel == GLB_ALL ? (1:length(e.channels)) : (channel:channel)
    sd = e.smear_data[e.channels[first(chans)].energy]
    energy = smearing == GLB_PRE ? sd.simbincenter : sd.bincenter
    res = Vector{Vector{Float64}}()
    for c in chans
        ch = e.channels[c]
        push!(
            res,
            [
                is_nosc(ch) ? 1.0 : glbProfileProbability(exp, ch.initial, ch.final, ch.cp, E) for
                E in energy
            ],
        )
    end
    print_channel_table(io, energy, res)
    return 0
end
glbShowChannelProbs(exp, channel, smearing, effi = GLB_W_EFF, bgi = GLB_W_BG) =
    glbShowChannelProbs(stdout, exp, channel, smearing, effi, bgi)

# ---------------------------------------------------------------------------
#  Experiment properties
# ---------------------------------------------------------------------------

glbGetNumberOfBins(exp) = (e = getexp(exp); e === nothing ? -1 : e.numofbins)
glbGetNumberOfSamplingPoints(exp) = (e = getexp(exp); e === nothing ? -1 : e.simbins)
glbGetNumberOfChannels(exp) = (e = getexp(exp); e === nothing ? -1 : length(e.channels))
glbGetNumberOfRules(exp) = (e = getexp(exp); e === nothing ? -1 : length(e.rules))
glbGetNumberOfFluxes(exp) = (e = getexp(exp); e === nothing ? -1 : length(e.fluxes))
glbGetBinSizeListPtr(exp) = (e = getexp(exp); e === nothing ? nothing : e.smear_data[1].binsize)
glbGetBinCentersListPtr(exp) =
    (e = getexp(exp); e === nothing ? nothing : e.smear_data[1].bincenter)
glbGetSamplingStepsizeListPtr(exp) =
    (e = getexp(exp); e === nothing ? nothing : e.smear_data[1].simbinsize)
glbGetSamplingPointsListPtr(exp) =
    (e = getexp(exp); e === nothing ? nothing : e.smear_data[1].simbincenter)

"""
    glbGetEminEmax(exp) -> (emin, emax)
"""
glbGetEminEmax(exp) = (e = getexp(exp); e === nothing ? nothing : (e.emin, e.emax))

"""
    glbGetEnergyWindow(exp, rule) -> (low, high)
"""
glbGetEnergyWindow(exp, rule) =
    (r = getrule(exp, rule); r === nothing ? nothing : (r.energy_window[1], r.energy_window[2]))

"""
    glbGetEnergyWindowBins(exp, rule) -> (first_bin, last_bin)   (1-based, inclusive)
"""
glbGetEnergyWindowBins(exp, rule) = (
    r = getrule(exp, rule);
    r === nothing ? nothing : (r.energy_window_bins[1], r.energy_window_bins[2])
)

"""
    glbSetEnergyWindow(exp, rule, low, high)

Change the analysis energy window of a rule (`exp`, `rule` may be `GLB_ALL`).
"""
function glbSetEnergyWindow(exp, rule, low, high)
    foreach_exp_rule(exp, rule) do i, j
        e = state().experiments[i]
        r = e.rules[j]
        lo, hi = max(min(low, high), e.emin), min(max(low, high), e.emax)
        r.energy_window .= (lo, hi)
        bc = e.smear_data[1].bincenter
        k = 1
        while k <= e.numofbins && bc[k] <= lo
            k += 1
        end
        r.energy_window_bins[1] = k
        while k <= e.numofbins && bc[k] < hi
            k += 1
        end
        r.energy_window_bins[2] = k - 1
        0
    end
end

"""
    glbNameToValue(exp, context, name) -> Int

1-based index of the AEDL object `name` (e.g. `"#CC"`) in `context`
(`"rule"`, `"channel"`, `"energy"`, `"cross"`, `"flux"`/`"nuflux"`), or -1.
"""
function glbNameToValue(exp, context, name)
    e = getexp(exp)
    e === nothing && return -1
    for (n, c, v) in e.names
        n == name && c == context && return v
    end
    return -1
end

function glbValueToName(exp, context, value)
    e = getexp(exp)
    e === nothing && return nothing
    for (n, c, v) in e.names
        v == value && c == context && return n
    end
    return nothing
end

function for_exps(f, exp)
    if exp == GLB_ALL
        foreach(f, state().experiments)
    else
        e = getexp(exp)
        e === nothing && return -1
        f(e)
    end
    return 0
end

function glbSetTargetMass(exp, mass)
    mass > 0 || (glb_error("Target mass has to be positive"); return -1)
    for_exps(e -> (e.targetmass = mass), exp)
end
glbGetTargetMass(exp) = (e = getexp(exp); e === nothing ? -1.0 : e.targetmass)

function set_flux_field(exp, flux_ident, field, value)
    for_exps(exp) do e
        fl = flux_ident == GLB_ALL ? e.fluxes : [e.fluxes[flux_ident]]
        foreach(f -> setfield!(f, field, Float64(value)), fl)
    end
end

function glbSetSourcePower(exp, flux_ident, power)
    power <= 0 && (glb_error("Power must be positive"); return -1)
    set_flux_field(exp, flux_ident, :target_power, power)
end
glbGetSourcePower(exp, flux_ident) =
    (e = getexp(exp); e === nothing ? -1.0 : e.fluxes[flux_ident].target_power)

function glbSetRunningTime(exp, flux_ident, time)
    time <= 0 && (glb_error("Time must be positive"); return -1)
    set_flux_field(exp, flux_ident, :time, time)
end
glbGetRunningTime(exp, flux_ident) =
    (e = getexp(exp); e === nothing ? -1.0 : e.fluxes[flux_ident].time)

glbSetFilterStateInExperiment(exp, on_off) = for_exps(e -> (e.filter_state = on_off), exp)
glbGetFilterStateInExperiment(exp) = (e = getexp(exp); e === nothing ? -1 : e.filter_state)
glbSetFilterState(on_off) = glbSetFilterStateInExperiment(GLB_ALL, on_off)
glbGetFilterState() = glbGetFilterStateInExperiment(1)
glbSetFilterInExperiment(exp, filter) = for_exps(e -> (e.filter_value = filter), exp)
glbGetFilterInExperiment(exp) = (e = getexp(exp); e === nothing ? -1.0 : e.filter_value)
glbSetFilter(filter) = glbSetFilterInExperiment(GLB_ALL, filter)
glbGetFilter() = glbGetFilterInExperiment(1)

"""
    glbFlux(exp, flux_ident, E, L, flavour, anti)

Flux of an experiment's source (per GeV, including running time and power).
"""
function glbFlux(exp, flux_ident, energy, distance, flavour, anti)
    e = getexp(exp)
    e === nothing && return -1.0
    return glb_flux_calc(energy, distance, flavour, anti, e.fluxes[flux_ident])
end

"""
    glbXSection(exp, xsec_ident, E, flavour, anti)

Cross section in 1e-38 cm².
"""
function glbXSection(exp, xsec_ident, energy, flavour, anti)
    e = getexp(exp)
    e === nothing && return -1.0
    return glb_xsec_calc(energy, flavour, anti, e.xsecs[xsec_ident])
end

# ---------------------------------------------------------------------------
#  Matter profiles
# ---------------------------------------------------------------------------

glbGetProfileData() = glbGetProfileDataInExperiment(1)

"""
    glbGetProfileDataInExperiment(exp) -> (lengths, densities)
"""
function glbGetProfileDataInExperiment(exp)
    e = getexp(exp)
    e === nothing && return nothing
    e.psteps > 0 || (glb_error("No profile data to get"); return nothing)
    return copy(e.lengthtab), copy(e.densitytab)
end

glbGetProfileTypeInExperiment(exp) = (e = getexp(exp); e === nothing ? -1 : e.density_profile_type)

"""
    glbSetProfileDataInExperiment(exp, lengths, densities)

Use an arbitrary layered matter profile (sets the profile type to 3).
"""
function glbSetProfileDataInExperiment(exp, lengths, densities)
    e = getexp(exp)
    e === nothing && return -1
    isempty(lengths) && (glb_error("No profile data to set"); return -1)
    e.psteps = length(lengths)
    e.density_profile_type = 3
    e.lengthtab = Float64.(collect(lengths))
    e.densitytab = Float64.(collect(densities))
    e.densitybuffer = copy(e.densitytab)
    e.baseline = sum(e.lengthtab)
    return 0
end

"""
    glbSetBaselineInExperiment(exp, L)

Change the baseline, recomputing the matter profile according to its type.
"""
function glbSetBaselineInExperiment(exp, baseline)
    e = getexp(exp)
    e === nothing && return -1
    t = e.density_profile_type
    if t == 1
        baseline > 0 || (glb_error("Baseline must be a positive number"); return -1)
        e.psteps = 1
        e.lengthtab, e.densitytab = glbAverageDensityProfile(baseline)
    elseif t == 2
        baseline > 0 || (glb_error("Baseline must be a positive number"); return -1)
        e.lengthtab, e.densitytab = glbStaceyProfile(baseline, e.psteps)
    elseif t == 3
        ttl = sum(e.lengthtab)
        e.lengthtab .*= baseline / ttl
    else
        glb_error("glbSetBaselineInExperiment failed")
        return -1
    end
    e.baseline = baseline
    e.densitybuffer = copy(e.densitytab)
    return 0
end
glbGetBaselineInExperiment(exp) = (e = getexp(exp); e === nothing ? -1.0 : e.baseline)

# ---------------------------------------------------------------------------
#  Oscillation probabilities
# ---------------------------------------------------------------------------

function prob_helper(initial, final, cp_sign, E, psteps, l, d, fs)
    P = zeros(3, 3)
    status = state().prob_func(P, cp_sign, E, psteps, l, d, fs, state().prob_user_data)
    status != 0 && (glb_error("Calculation of oscillation probabilities failed."); return -1.0)
    return P[initial, final]
end

"""
    glbVacuumProbability(initial, final, cp_sign, E, L)

Vacuum oscillation probability (flavours 1=e, 2=μ, 3=τ; `cp_sign=±1`,
E in GeV, L in km) for the current oscillation parameters.
"""
glbVacuumProbability(initial, final, cp_sign, E, L) =
    prob_helper(initial, final, cp_sign, E, 1, [L], [0.0], -1.0)

"""
    glbConstantDensityProbability(initial, final, cp_sign, E, L, rho)

Oscillation probability in matter of constant density `rho` (g/cm³).
"""
glbConstantDensityProbability(initial, final, cp_sign, E, L, rho) =
    prob_helper(initial, final, cp_sign, E, 1, [L], [rho], -1.0)

"""
    glbProfileProbability(exp, initial, final, cp_sign, E)

Probability with the matter profile of experiment `exp` (filter ignored).
"""
function glbProfileProbability(exp, initial, final, cp_sign, E)
    e = getexp(exp)
    e === nothing && return -1.0
    set_profile_scaling!(e, 1.0)
    return prob_helper(initial, final, cp_sign, E, e.psteps, e.lengthtab, e.densitybuffer, -1.0)
end

"""
    glbFilteredConstantDensityProbability(exp, initial, final, cp_sign, E)

Like `glbProfileProbability`, including the experiment's low-pass filter.
"""
function glbFilteredConstantDensityProbability(exp, initial, final, cp_sign, E)
    e = getexp(exp)
    e === nothing && return -1.0
    set_profile_scaling!(e, 1.0)
    fs = e.filter_state == GLB_ON ? e.filter_value : -1.0
    return prob_helper(initial, final, cp_sign, E, e.psteps, e.lengthtab, e.densitybuffer, fs)
end
