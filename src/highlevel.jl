# The high-level Julia API.
#
# It is built on the glb* functions, which mirror the C library, but uses
# names instead of index constants, keyword arguments, return values instead
# of status codes, and exceptions instead of printed errors. Every function
# works on the current setup (see `Setup`), or on the one passed as `setup`.

# ---------------------------------------------------------------------------
#  Setups
# ---------------------------------------------------------------------------

"""
    Setup()

Create a new, empty setup: no experiments, the standard probability engine,
the built-in χ² functions and priors.
"""
function Setup()
    s = Setup(new_state())
    with_setup(glbInit, s)
    return s
end

"""
    current_setup() -> Setup

The setup that functions use when no `setup` keyword is given.
"""
current_setup() = CURRENT_SETUP[]

"""
    with_setup(f, setup)

Call `f()` with `setup` as the current setup, and restore the previous one
afterwards. Setups are not thread-safe: use one setup per thread at a time.

```julia
far = Setup()
with_setup(far) do
    Experiment("NOvA.glb"; BASELINE=1500)
    simulate!(truth)
end
```
"""
function with_setup(f, s::Setup)
    old = CURRENT_SETUP[]
    CURRENT_SETUP[] = s
    try
        return f()
    finally
        CURRENT_SETUP[] = old
    end
end

# Run f in setup s with errors raised as exceptions
in_setup(f, s::Setup) = with_setup(() -> strict(f), s)

function Base.show(io::IO, s::Setup)
    n = length(s.state.experiments)
    print(
        io,
        "Setup(",
        n,
        n == 1 ? " experiment" : " experiments",
        s === CURRENT_SETUP[] ? ", current)" : ")",
    )
end

# ---------------------------------------------------------------------------
#  Oscillation parameters
# ---------------------------------------------------------------------------

# Canonical names of the standard parameters and accepted ASCII spellings
const PARAM_NAMES = (:θ12, :θ13, :θ23, :δ, :Δm21, :Δm31)
const PARAM_ALIASES = Dict{Symbol, Symbol}(
    :theta12 => :θ12,
    :theta13 => :θ13,
    :theta23 => :θ23,
    :delta => :δ,
    :deltacp => :δ,
    :δCP => :δ,
    :dcp => :δ,
    :dm21 => :Δm21,
    :dm31 => :Δm31,
    :Δm²21 => :Δm21,
    :Δm²31 => :Δm31,
    :rho => :ρ,
)

canonical(name::Symbol) = get(PARAM_ALIASES, name, name)

"""
    OscParams(; θ12, θ13, θ23, δ=0, Δm21, Δm31, ρ=1, extra=())
    OscParams(p::OscParams; changes...)

Neutrino oscillation parameters: mixing angles and the CP phase in radians,
mass squared differences Δm²21 and Δm²31 in eV², and the matter density
scaling ρ (1 = the density profile of the experiment). ASCII names work too:
`theta12`, `theta13`, `theta23`, `delta`, `dm21`, `dm31`, `rho`.

`ρ` is one number for all experiments or a vector with one entry each.
Additional parameters of a user-defined probability engine go in `extra`, as
a named tuple. The second form copies `p` with some parameters changed.

```julia
truth = OscParams(θ12=asin(sqrt(0.307)), θ13=asin(sqrt(0.022)), θ23=asin(sqrt(0.55)),
                  δ=1.2π, Δm21=7.4e-5, Δm31=2.5e-3)
test = OscParams(truth; δ=0, Δm31=-2.4e-3)
truth.θ23, truth.sin²θ23, truth.sin²2θ13
```
"""
struct OscParams
    values::Vector{Float64}        # θ12, θ13, θ23, δ, Δm21, Δm31, extra parameters...
    extra_names::Vector{Symbol}
    densities::Vector{Float64}     # ρ: one entry, or one per experiment
end

function OscParams(; kwargs...)
    kw = Dict(canonical(k) => v for (k, v) in kwargs)
    missing_names = [n for n in (:θ12, :θ13, :θ23, :Δm21, :Δm31) if !haskey(kw, n)]
    isempty(missing_names) || throw(
        ArgumentError("OscParams needs $(join(missing_names, ", ")); δ defaults to 0 and ρ to 1"),
    )
    extra = get(kw, :extra, NamedTuple())
    for k in keys(kw)
        k in PARAM_NAMES ||
            k in (:ρ, :extra) ||
            throw(
                ArgumentError(
                    "unknown oscillation parameter $k; use `extra` for non-standard parameters",
                ),
            )
    end
    vals = Float64[kw[:θ12], kw[:θ13], kw[:θ23], get(kw, :δ, 0.0), kw[:Δm21], kw[:Δm31]]
    append!(vals, Float64.(values(extra)))
    return OscParams(vals, collect(Symbol, keys(extra)), Float64.(vcat(get(kw, :ρ, 1.0))))
end

function OscParams(p::OscParams; kwargs...)
    vals = copy(p.values)
    ρ = copy(p.densities)
    extra_names = copy(p.extra_names)
    for (k, v) in kwargs
        k = canonical(k)
        if k == :ρ
            ρ = Float64.(vcat(v))
        elseif k == :extra
            for (n, x) in pairs(v)
                n in extra_names ? (vals[6 + findfirst(==(n), extra_names)] = x) :
                (push!(extra_names, n); push!(vals, x))
            end
        else
            vals[param_index(OscParams(vals, extra_names, ρ), k)] = v
        end
    end
    return OscParams(vals, extra_names, ρ)
end

function param_index(p::OscParams, name::Symbol)
    name = canonical(name)
    i = findfirst(==(name), PARAM_NAMES)
    i !== nothing && return i
    j = findfirst(==(name), p.extra_names)
    j !== nothing && return 6 + j
    throw(ArgumentError("unknown oscillation parameter $name"))
end

# Mixing angles are often quoted as sin²θ or sin²2θ
const DERIVED = Dict{Symbol, Tuple{Symbol, Function, Function}}(
    # name => (angle, angle -> value, value -> angle)
    Symbol("sin²θ12") => (:θ12, θ -> sin(θ)^2, x -> asin(sqrt(x))),
    Symbol("sin²θ13") => (:θ13, θ -> sin(θ)^2, x -> asin(sqrt(x))),
    Symbol("sin²θ23") => (:θ23, θ -> sin(θ)^2, x -> asin(sqrt(x))),
    Symbol("sin²2θ12") => (:θ12, θ -> sin(2θ)^2, x -> asin(sqrt(x)) / 2),
    Symbol("sin²2θ13") => (:θ13, θ -> sin(2θ)^2, x -> asin(sqrt(x)) / 2),
    Symbol("sin²2θ23") => (:θ23, θ -> sin(2θ)^2, x -> asin(sqrt(x)) / 2),
)
for (ascii, name) in (
    :sin2theta12 => "sin²θ12",
    :sin2theta13 => "sin²θ13",
    :sin2theta23 => "sin²θ23",
    :sin22theta12 => "sin²2θ12",
    :sin22theta13 => "sin²2θ13",
    :sin22theta23 => "sin²2θ23",
)
    PARAM_ALIASES[ascii] = Symbol(name)
end

function Base.getproperty(p::OscParams, name::Symbol)
    name in fieldnames(OscParams) && return getfield(p, name)
    name = canonical(name)
    name == :ρ && return length(p.densities) == 1 ? p.densities[1] : copy(p.densities)
    if haskey(DERIVED, name)
        angle, f, _ = DERIVED[name]
        return f(getfield(p, :values)[param_index(p, angle)])
    end
    return getfield(p, :values)[param_index(p, name)]
end

Base.propertynames(p::OscParams) = (PARAM_NAMES..., p.extra_names..., :ρ, keys(DERIVED)...)

"""
    with_param(p::OscParams, name, value) -> OscParams

Copy of `p` with one parameter changed. `name` can also be one of
`:sin²θ12`, `:sin²θ13`, `:sin²θ23`, `:sin²2θ12`, `:sin²2θ13`, `:sin²2θ23`
(ASCII: `:sin2theta13`, `:sin22theta13`, ...), which sets the angle.
"""
function with_param(p::OscParams, name::Symbol, value)
    name = canonical(name)
    if haskey(DERIVED, name)
        angle, _, inv = DERIVED[name]
        return OscParams(p; (angle => inv(value),)...)
    end
    return OscParams(p; (name => value,)...)
end

Base.:(==)(a::OscParams, b::OscParams) =
    a.values == b.values && a.extra_names == b.extra_names && a.densities == b.densities
Base.isapprox(a::OscParams, b::OscParams; kw...) =
    a.extra_names == b.extra_names &&
    isapprox(a.values, b.values; kw...) &&
    isapprox(a.densities, b.densities; kw...)

function Base.show(io::IO, p::OscParams)
    v = p.values
    @printf(io, "OscParams(θ12=%.5g, θ13=%.5g, θ23=%.5g, δ=%.5g, Δm21=%.5g, Δm31=%.5g", v[1:6]...)
    for (i, n) in enumerate(p.extra_names)
        @printf(io, ", %s=%.5g", n, v[6 + i])
    end
    ρ =
        length(p.densities) == 1 ? @sprintf("%.5g", p.densities[1]) :
        "[" * join((@sprintf("%.5g", x) for x in p.densities), ", ") * "]"
    print(io, ", ρ=", ρ, ")")
end

function Base.show(io::IO, ::MIME"text/plain", p::OscParams)
    v = p.values
    println(io, "OscParams")
    @printf(io, "  θ12 = %.6g  (sin²θ12 = %.6g)\n", v[1], sin(v[1])^2)
    @printf(io, "  θ13 = %.6g  (sin²2θ13 = %.6g)\n", v[2], sin(2v[2])^2)
    @printf(io, "  θ23 = %.6g  (sin²θ23 = %.6g)\n", v[3], sin(v[3])^2)
    @printf(io, "  δ   = %.6g  (%.4g π)\n", v[4], v[4] / π)
    @printf(io, "  Δm21 = %.6g eV²\n", v[5])
    @printf(io, "  Δm31 = %.6g eV²\n", v[6])
    for (i, n) in enumerate(p.extra_names)
        @printf(io, "  %s = %.6g\n", n, v[6 + i])
    end
    print(io, "  ρ   = ", length(p.densities) == 1 ? p.densities[1] : p.densities)
end

# Conversion to and from the parameter vectors of the engine
function to_glb(p::OscParams)
    st = state()
    length(p.values) == st.oscp || throw(
        GLoBESError(
            "the probability engine has $(st.oscp) parameters, these OscParams have $(length(p.values))",
        ),
    )
    g = glbAllocParams()
    g.osc .= p.values
    n = length(st.experiments)
    if length(p.densities) == 1
        fill!(g.density, p.densities[1])
    elseif length(p.densities) == n
        g.density .= p.densities
    else
        throw(GLoBESError("ρ has $(length(p.densities)) entries for $n experiments"))
    end
    return g
end

function from_glb(g::GlbParams, like::OscParams)
    ρ = isempty(g.density) ? copy(like.densities) : copy(g.density)
    length(ρ) > 1 && all(==(ρ[1]), ρ) && length(like.densities) == 1 && (ρ = ρ[1:1])
    return OscParams(copy(g.osc), copy(like.extra_names), ρ)
end

# ---------------------------------------------------------------------------
#  Experiments
# ---------------------------------------------------------------------------

"""
    Experiment(file; setup=current_setup(), variables...)

Load an experiment from an AEDL file and add it to the setup. Keyword
arguments define AEDL variables used in the file, numbers or lists:

```julia
nova = Experiment("NOvA.glb")
nf = Experiment("NFvar.glb"; BASELINE=3000, emax=50)
```

The file is looked up in the current directory and in the directories listed
in the environment variable `GLB_PATH`. Files it includes and its flux and
cross section files are also looked up next to it. A file that defines
several experiments (separated by `#NEXT#`) needs [`load_experiments`](@ref).
"""
struct Experiment
    setup::Setup
    data::ExperimentData
end

function Experiment(file::AbstractString; setup::Setup = current_setup(), variables...)
    exps = load_experiments(file; setup, variables...)
    length(exps) == 1 ||
        throw(GLoBESError("$file defines $(length(exps)) experiments; use load_experiments"))
    return exps[1]
end

"""
    load_experiments(file; setup=current_setup(), variables...) -> Vector{Experiment}

Load all experiments defined in an AEDL file; see [`Experiment`](@ref).
"""
function load_experiments(file::AbstractString; setup::Setup = current_setup(), variables...)
    in_setup(setup) do
        vars, lists = copy(aedl_presym_vars), copy(aedl_presym_lists)
        try
            for (k, v) in variables
                v isa AbstractVector ?
                glbDefineAEDLList("%" * lstrip(String(k), '%'), Float64.(v)) :
                glbDefineAEDLVariable(String(k), Float64(v))
            end
            n = length(state().experiments)
            status = glbInitExperiment(file)
            status == -2 && throw(GLoBESError("could not load $file"))
            return [Experiment(setup, e) for e in state().experiments[(n + 1):end]]
        finally
            copy!(aedl_presym_vars, vars)
            copy!(aedl_presym_lists, lists)
        end
    end
end

"""
    experiments(; setup=current_setup()) -> Vector{Experiment}

All experiments of a setup, in the order they were loaded.
"""
experiments(; setup::Setup = current_setup()) =
    [Experiment(setup, e) for e in setup.state.experiments]

"""
    clear_experiments!(; setup=current_setup())

Remove all experiments from a setup.
"""
clear_experiments!(; setup::Setup = current_setup()) =
    (in_setup(glbClearExperimentList, setup); nothing)

# Position of an experiment in its setup, as used by the glb* functions
function index(e::Experiment)
    i = findfirst(x -> x === e.data, e.setup.state.experiments)
    i === nothing && throw(GLoBESError("experiment $(name(e)) is no longer part of its setup"))
    return i
end

Base.:(==)(a::Experiment, b::Experiment) = a.data === b.data
Base.hash(e::Experiment, h::UInt) = hash(objectid(e.data), h)

# Run f(i) in the experiment's setup, where i is its index
on(f, e::Experiment) = in_setup(() -> f(index(e)), e.setup)

"""
    name(e::Experiment) -> String

The name of the AEDL file the experiment was loaded from.
"""
name(e::Experiment) = splitext(basename(e.data.filename))[1]

function Base.show(io::IO, e::Experiment)
    d = e.data
    print(
        io,
        "Experiment(\"",
        name(e),
        "\", L = ",
        d.baseline,
        " km, ",
        d.numofbins,
        " bins from ",
        d.emin,
        " to ",
        d.emax,
        " GeV, ",
        length(d.rules),
        length(d.rules) == 1 ? " rule)" : " rules)",
    )
end

function Base.show(io::IO, ::MIME"text/plain", e::Experiment)
    d = e.data
    println(io, "Experiment \"", name(e), "\" (", d.filename, ")")
    println(io, "  baseline     ", d.baseline, " km")
    println(io, "  target mass  ", d.targetmass, " kt")
    @printf(io, "  energy bins  %d from %g to %g GeV\n", d.numofbins, d.emin, d.emax)
    println(io, "  channels     ", length(d.channels))
    print(io, "  rules        ", join(rule_names(e), ", "))
end

"""
    baseline(e) -> Float64
    baseline!(e, L)

Baseline in km. Changing it recomputes the matter density profile in the way
the AEDL file defines it (average density, PREM layers, or rescaled layers).
"""
baseline(e::Experiment) = e.data.baseline
baseline!(e::Experiment, L::Real) = (on(i -> glbSetBaselineInExperiment(i, Float64(L)), e); e)

"""
    target_mass(e) -> Float64
    target_mass!(e, mass)

Fiducial detector mass in kt (the unit used in AEDL files).
"""
target_mass(e::Experiment) = e.data.targetmass
target_mass!(e::Experiment, m::Real) = (on(i -> glbSetTargetMass(i, Float64(m)), e); e)

"""
    running_time(e; flux=1) -> Float64
    running_time!(e, years; flux=:all)

Running time in years of one flux, or of all fluxes.
"""
running_time(e::Experiment; flux::Integer = 1) = e.data.fluxes[flux].time
function running_time!(e::Experiment, t::Real; flux = :all)
    on(e) do i
        for f in (flux === :all ? eachindex(e.data.fluxes) : (flux,))
            glbSetRunningTime(i, f, Float64(t))
        end
    end
    return e
end

"""
    source_power(e; flux=1) -> Float64
    source_power!(e, power; flux=:all)

Source power (stored muons or ions per year, or beam power in MW; see the
`@power` and `@stored_muons` entries of the AEDL file).
"""
source_power(e::Experiment; flux::Integer = 1) = e.data.fluxes[flux].target_power
function source_power!(e::Experiment, p::Real; flux = :all)
    on(e) do i
        for f in (flux === :all ? eachindex(e.data.fluxes) : (flux,))
            glbSetSourcePower(i, f, Float64(p))
        end
    end
    return e
end

"""
    energy_range(e) -> (emin, emax)
    bin_centers(e) -> Vector{Float64}
    bin_widths(e) -> Vector{Float64}

The reconstructed energy range and bins of an experiment, in GeV.
"""
energy_range(e::Experiment) = (e.data.emin, e.data.emax)
bin_centers(e::Experiment) = copy(e.data.smear_data[1].bincenter)
bin_widths(e::Experiment) = copy(e.data.smear_data[1].binsize)

"""
    density_profile(e) -> (lengths, densities)
    density_profile!(e, lengths, densities)

The matter density profile: layers with lengths in km and densities in g/cm³.
See also [`prem_profile`](@ref) and [`average_profile`](@ref).
"""
density_profile(e::Experiment) = (copy(e.data.lengthtab), copy(e.data.densitytab))
function density_profile!(e::Experiment, lengths::AbstractVector, densities::AbstractVector)
    length(lengths) == length(densities) ||
        throw(ArgumentError("lengths and densities differ in length"))
    on(i -> glbSetProfileDataInExperiment(i, Float64.(lengths), Float64.(densities)), e)
    return e
end

"""
    prem_profile(L, layers) -> (lengths, densities)
    average_profile(L) -> (lengths, densities)

Matter density along a baseline of `L` km through the Earth (PREM): in
`layers` steps of equal length, or a single layer with the average density.
"""
prem_profile(L::Real, layers::Integer) = glbStaceyProfile(Float64(L), layers)
average_profile(L::Real) = glbAverageDensityProfile(Float64(L))
@doc (@doc prem_profile) average_profile

"""
    rule_names(e) -> Vector{String}
    channel_names(e) -> Vector{String}

Names of the rules and channels defined in the AEDL file, without the leading
`#` (their numbers if they have no name).
"""
rule_names(e::Experiment) = names_in(e, "rule", length(e.data.rules))
channel_names(e::Experiment) = names_in(e, "channel", length(e.data.channels))

function names_in(e::Experiment, context, n)
    out = string.(1:n)
    for (nm, c, v) in e.data.names
        c == context && 1 <= v <= n && (out[v] = lstrip(nm, '#'))
    end
    return out
end

# A rule given by number or name
rule_index(e::Experiment, r::Integer) =
    1 <= r <= length(e.data.rules) ? Int(r) : throw(GLoBESError("$(name(e)) has no rule $r"))
function rule_index(e::Experiment, r::Union{AbstractString, Symbol})
    i = findfirst(==(lstrip(String(r), '#')), rule_names(e))
    i === nothing &&
        throw(GLoBESError("$(name(e)) has no rule \"$r\"; it has $(join(rule_names(e), ", "))"))
    return i
end
channel_index(e::Experiment, c::Integer) =
    1 <= c <= length(e.data.channels) ? Int(c) : throw(GLoBESError("$(name(e)) has no channel $c"))
function channel_index(e::Experiment, c::Union{AbstractString, Symbol})
    i = findfirst(==(lstrip(String(c), '#')), channel_names(e))
    i === nothing && throw(GLoBESError("$(name(e)) has no channel \"$c\""))
    return i
end

# ---------------------------------------------------------------------------
#  Simulated data
# ---------------------------------------------------------------------------

"""
    simulate!(truth::OscParams; setup=current_setup())

Compute the simulated ("true") event rates of all experiments for the
oscillation parameters `truth`. Every χ² compares against these rates.
"""
function simulate!(truth::OscParams; setup::Setup = current_setup())
    in_setup(setup) do
        isempty(state().experiments) && throw(GLoBESError("no experiments loaded"))
        g = to_glb(truth)
        glbSetOscillationParameters(g)
        glbSetRates()
        state().truth = g
    end
    return nothing
end

"""
    simulated_params(; setup=current_setup()) -> OscParams

The parameters passed to the last [`simulate!`](@ref).
"""
function simulated_params(; setup::Setup = current_setup())
    g = setup.state.truth
    g === nothing && throw(GLoBESError("no data simulated yet; call simulate!"))
    return from_glb(g, OscParams(copy(g.osc), engine_param_names(setup), [1.0]))
end

engine_param_names(setup::Setup) = setup.state.extra_names

"""
    event_rates(e; rule=1) -> NamedTuple
    event_rates(e, params::OscParams; rule=1) -> NamedTuple

Event rates of one rule per energy bin: those of the simulated data, or those
predicted for `params`. Returns `(energy, width, signal, background, total)`;
`total` is signal plus background (with the background normalisation of the
rule) inside the rule's energy window and zero outside.
"""
function event_rates(e::Experiment; rule = 1)
    on(e) do i
        r = rule_index(e, rule)
        d = e.data
        state().truth === nothing && throw(GLoBESError("no data simulated yet; call simulate!"))
        (
            energy = bin_centers(e),
            width = bin_widths(e),
            signal = copy(d.SignalRates[r]),
            background = d.rules[r].bg_centers[1] .* d.BackgroundRates[r],
            total = copy(d.rates0[r]),
        )
    end
end

function event_rates(e::Experiment, params::OscParams; rule = 1)
    on(e) do i
        r = rule_index(e, rule)
        d = e.data
        g = to_glb(params)
        saved = glbAllocParams()
        glbGetOscillationParameters(saved)
        try
            glbSetOscillationParameters(g)
            set_profile_scaling!(d, g.density[i])
            set_new_rates!(d)
            sig = copy(d.rates1[r])
            bg = d.rules[r].bg_centers[1] .* d.rates1BG[r]
            lo, hi = d.rules[r].energy_window_bins
            total = zeros(length(sig))
            total[lo:hi] .= sig[lo:hi] .+ bg[lo:hi]
            (
                energy = bin_centers(e),
                width = bin_widths(e),
                signal = sig,
                background = bg,
                total = total,
            )
        finally
            glbSetOscillationParameters(saved)
        end
    end
end

"""
    total_events(e; rule=1) -> Float64
    total_events(e, params; rule=1) -> Float64

Number of events of one rule in its energy window, signal plus background.
"""
total_events(e::Experiment; rule = 1) = sum(event_rates(e; rule).total)
total_events(e::Experiment, params::OscParams; rule = 1) = sum(event_rates(e, params; rule).total)

"""
    channel_rates(e, channel; smeared=true) -> Vector{Float64}

Simulated event rates of one channel, after (`smeared=true`, per energy bin)
or before energy smearing (per sampling point).
"""
channel_rates(e::Experiment, channel; smeared::Bool = true) = on(
    i -> copy(glbGetChannelRatePtr(i, channel_index(e, channel), smeared ? GLB_POST : GLB_PRE)),
    e,
)

"""
    fit_rates(e; rule=1) -> (signal, background)

The predicted signal and background rates of the χ² evaluation in progress,
without systematics. For use inside user-defined χ² functions (see
[`define_chi2`](@ref)); `background` is not yet multiplied by the background
normalisation. The vectors are the engine's buffers: do not modify them.
"""
function fit_rates(e::Experiment; rule = 1)
    r = rule_index(e, rule)
    return (signal = e.data.rates1[r], background = e.data.rates1BG[r])
end

"""
    true_rates(e; rule=1) -> Vector{Float64}

The simulated rates of a rule inside its energy window, signal plus
background, as used by the χ² functions. The engine's buffer: do not modify.
"""
true_rates(e::Experiment; rule = 1) = e.data.rates0[rule_index(e, rule)]

"""
    energy_window(e; rule=1) -> (emin, emax)
    energy_window!(e, emin, emax; rule=:all)
    energy_window_bins(e; rule=1) -> UnitRange

The energy window of the analysis of a rule, in GeV, and the bins in it.
"""
energy_window(e::Experiment; rule = 1) = Tuple(e.data.rules[rule_index(e, rule)].energy_window)
function energy_window!(e::Experiment, lo::Real, hi::Real; rule = :all)
    on(
        i -> glbSetEnergyWindow(
            i,
            rule === :all ? GLB_ALL : rule_index(e, rule),
            Float64(lo),
            Float64(hi),
        ),
        e,
    )
    return e
end
energy_window_bins(e::Experiment; rule = 1) =
    (b = e.data.rules[rule_index(e, rule)].energy_window_bins; b[1]:b[2])

# ---------------------------------------------------------------------------
#  Probabilities
# ---------------------------------------------------------------------------

const FLAVOURS = Dict{Symbol, Int}(:e => 1, :μ => 2, :mu => 2, :τ => 3, :tau => 3)
flavour(f::Integer) = 1 <= f <= 3 ? Int(f) : throw(ArgumentError("flavour must be 1, 2 or 3"))
flavour(f::Symbol) =
    haskey(FLAVOURS, f) ? FLAVOURS[f] : throw(ArgumentError("unknown flavour $f; use :e, :μ or :τ"))

"""
    probability(params, from, to; E, L, ρ=0, antineutrino=false)
    probability(params, from, to, e::Experiment; E, antineutrino=false)

Oscillation probability P(ν_from → ν_to) at neutrino energy `E` in GeV
(a number or a vector): over `L` km of constant matter density `ρ` in g/cm³
(vacuum for ρ = 0), or along the baseline and matter profile of an
experiment. Flavours are `:e`, `:μ` (`:mu`), `:τ` (`:tau`) or 1, 2, 3.

```julia
probability(truth, :μ, :e; E=2.0, L=1300, ρ=2.8)
probability(truth, :μ, :e, dune; E=0.5:0.1:5, antineutrino=true)
```
"""
function probability(
    params::OscParams,
    from,
    to;
    E,
    L::Real,
    ρ::Real = 0.0,
    antineutrino::Bool = false,
    setup::Setup = current_setup(),
)
    a, b, cp = flavour(from), flavour(to), antineutrino ? -1 : +1
    with_engine_params(params, setup) do
        P(En) =
            ρ == 0 ? glbVacuumProbability(a, b, cp, En, Float64(L)) :
            glbConstantDensityProbability(a, b, cp, En, Float64(L), Float64(ρ))
        E isa Number ? P(Float64(E)) : [P(Float64(x)) for x in E]
    end
end

function probability(params::OscParams, from, to, e::Experiment; E, antineutrino::Bool = false)
    a, b, cp = flavour(from), flavour(to), antineutrino ? -1 : +1
    with_engine_params(params, e.setup) do
        i = index(e)
        P(En) = glbProfileProbability(i, a, b, cp, En)
        E isa Number ? P(Float64(E)) : [P(Float64(x)) for x in E]
    end
end

# Evaluate f with the engine's oscillation parameters set to params
function with_engine_params(f, params::OscParams, setup::Setup)
    in_setup(setup) do
        saved = glbAllocParams()
        glbGetOscillationParameters(saved)
        g = glbAllocParams()
        g.osc .= params.values
        length(params.values) == state().oscp || throw(
            GLoBESError(
                "the probability engine has $(state().oscp) parameters, these OscParams have $(length(params.values))",
            ),
        )
        glbSetOscillationParameters(g)
        try
            return f()
        finally
            all(isnan, saved.osc) || glbSetOscillationParameters(saved)
        end
    end
end

# ---------------------------------------------------------------------------
#  χ² with systematics only
# ---------------------------------------------------------------------------

"""
    chi2(test::OscParams; experiment=:all, rule=:all, systematics=true, setup=current_setup())

χ² of the simulated data against the prediction for `test`, minimised over
the systematics (nuisance) parameters of each rule but not over oscillation
parameters (`glbChiSys`). `experiment` is an [`Experiment`](@ref) or `:all`,
`rule` a rule number or name or `:all`. `systematics=false` uses the rules'
χ² functions without systematics (statistics only, by default).
"""
function chi2(
    test::OscParams;
    experiment = :all,
    rule = :all,
    systematics::Bool = true,
    setup::Setup = current_setup(),
)
    setup = experiment === :all ? setup : experiment.setup
    in_setup(setup) do
        check_simulated()
        g = to_glb(test)
        exp_i = experiment === :all ? GLB_ALL : index(experiment)
        rule_i =
            rule === :all ? GLB_ALL : experiment === :all ? Int(rule) : rule_index(experiment, rule)
        systematics && return glbChiSys(g, exp_i, rule_i)
        saved = [r.sys_on_off for d in state().experiments for r in d.rules]
        try
            glbSwitchSystematics(GLB_ALL, GLB_ALL, GLB_OFF)
            return glbChiSys(g, exp_i, rule_i)
        finally
            k = 0
            for d in state().experiments, r in d.rules
                r.sys_on_off = saved[k += 1]
            end
        end
    end
end

check_simulated() =
    state().truth === nothing && throw(GLoBESError("no data simulated yet; call simulate!"))

# ---------------------------------------------------------------------------
#  Fits: χ² minimised over oscillation parameters
# ---------------------------------------------------------------------------

"""
    FitResult

Result of [`marginalize`](@ref): `chi2` is the minimum, `params` the
parameters at the minimum, `free` the parameters that were fitted and
`iterations` the number of χ² evaluations.
"""
struct FitResult
    chi2::Float64
    params::OscParams
    free::Vector{Symbol}
    iterations::Int
end

function Base.show(io::IO, ::MIME"text/plain", r::FitResult)
    @printf(
        io,
        "FitResult: χ² = %.6g after %d evaluations, free: %s\n",
        r.chi2,
        r.iterations,
        join(r.free, ", ")
    )
    show(io, MIME"text/plain"(), r.params)
end

"""
    marginalize(test; free=nothing, fixed=nothing, priors=(;), centers=simulated_params(),
                prior=nothing, experiment=:all, setup=current_setup()) -> FitResult

Minimise χ² over some oscillation parameters, starting from `test`, with the
others fixed at their values in `test` ("projection" in the GLoBES manual).

* `free` or `fixed`: the parameters to fit, or those to keep; by default all
  are fitted. Names are `:θ12, :θ13, :θ23, :δ, :Δm21, :Δm31`, `:ρ` for the
  matter densities of all experiments, and the names of extra parameters.
* `priors`: Gaussian priors as standard deviations, `(θ12=0.03, ρ=0.05)`;
  they add `((x - center)/σ)²` to χ² for fitted parameters.
* `centers`: the centres of the priors, the simulated parameters by default.
* `prior`: a function `params::OscParams -> Float64` added to χ², for
  priors that are not Gaussian.
* `experiment`: fit one experiment only.

```julia
r = marginalize(OscParams(truth; δ=0); fixed=(:θ13, :δ), priors=(θ12=0.1truth.θ12, ρ=0.05))
r.chi2, r.params.θ23
```
"""
function marginalize(
    test::OscParams;
    free = nothing,
    fixed = nothing,
    priors = NamedTuple(),
    centers = nothing,
    prior = nothing,
    experiment = :all,
    setup::Setup = current_setup(),
)
    setup = experiment === :all ? setup : experiment.setup
    in_setup(setup) do
        check_simulated()
        st = state()
        n_exp = length(st.experiments)
        g = to_glb(test)
        free_flags = projection_flags(test, free, fixed, n_exp)
        c = centers === nothing ? glbCopyParams(st.truth) : to_glb(centers)
        err = glbAllocParams()
        fill!(err.osc, 0.0)
        fill!(err.density, 0.0)
        for (k, σ) in pairs(priors)
            k = canonical(k)
            if k == :ρ
                fill!(err.density, Float64(σ))
            else
                err.osc[param_index(test, k)] = σ
            end
        end

        # Save what the glb* layer has set, so that it is unchanged afterwards
        saved_proj =
            (copy(st.para_tab), copy(st.index_tab), st.n_free, st.n_fix, copy(st.fix_params))
        saved_centers = glbCopyParams(st.central_values)
        saved_errors = glbCopyParams(st.input_errors)
        saved_prior =
            (st.prior_func, st.starting_values_hook, st.input_errors_hook, st.prior_user_data)
        try
            proj = glbAllocProjection()
            proj.osc .= free_flags[1:st.oscp]
            proj.density .= free_flags[(st.oscp + 1):end]
            glbSetProjection(proj)
            glbSetCentralValues(c)
            glbSetInputErrors(err)
            if prior !== nothing
                user = p -> Float64(prior(from_glb(p, test)))
                glbRegisterPriorFunction(
                    (p, ud) -> builtin_prior_impl(p) + user(p),
                    nothing,
                    nothing,
                )
                glbSetCentralValues(c)
                glbSetInputErrors(err)
            end
            out = glbAllocParams()
            exp_i = experiment === :all ? GLB_ALL : index(experiment)
            χ² = glbChiNP(g, out, exp_i)
            if experiment !== :all             # only this experiment's density was fitted
                out.density .= g.density
                out.density[exp_i] = glbGetDensityParams(out, exp_i)
            end
            names = [n for (n, f) in zip(free_names(test, n_exp), free_flags) if f == GLB_FREE]
            return FitResult(χ², from_glb(out, test), unique(names), abs(out.iterations))
        finally
            st.prior_func, st.starting_values_hook, st.input_errors_hook, st.prior_user_data =
                saved_prior
            st.para_tab, st.index_tab, st.n_free, st.n_fix, st.fix_params = saved_proj
            glbCopyParams(saved_centers, st.central_values)
            glbCopyParams(saved_errors, st.input_errors)
        end
    end
end

free_names(test::OscParams, n_exp) = [PARAM_NAMES..., test.extra_names..., fill(:ρ, n_exp)...]

function projection_flags(test::OscParams, free, fixed, n_exp)
    free === nothing ||
        fixed === nothing ||
        throw(ArgumentError("give either free or fixed, not both"))
    n = length(test.values)
    names = Set(canonical.(collect(something(free, fixed, ()))))
    for k in names
        k == :ρ || param_index(test, k)          # throws for unknown names
    end
    flag(k) =
        free !== nothing ? (k in names ? GLB_FREE : GLB_FIXED) :
        fixed !== nothing ? (k in names ? GLB_FIXED : GLB_FREE) : GLB_FREE
    return [[flag(k) for k in free_names(test, 0)]; fill(flag(:ρ), n_exp)]
end

"""
    minimizer!(kind; setup=current_setup())

`:nested` (default): Powell's method over the oscillation parameters, with the
systematics minimised separately for every rule at each step. `:hybrid`: one
Powell minimisation over oscillation and systematics parameters together.
"""
function minimizer!(kind::Symbol; setup::Setup = current_setup())
    id =
        kind === :nested ? GLB_MIN_NESTED_POWELL :
        kind === :hybrid ? GLB_MIN_POWELL :
        throw(ArgumentError("minimizer must be :nested or :hybrid"))
    in_setup(() -> glbSelectMinimizer(id), setup)
    return nothing
end

# ---------------------------------------------------------------------------
#  Scans
# ---------------------------------------------------------------------------

"""
    chi2_scan(base, x => xs; kwargs...) -> Vector{Float64}
    chi2_scan(base, x => xs, y => ys; kwargs...) -> Matrix{Float64}

χ² on a grid in one or two parameters, with the other parameters as in
`base`. Parameter names are as in [`with_param`](@ref), including
`:sin²2θ13` and `:sin²θ23`. Without keyword arguments the χ² is minimised
over systematics only ([`chi2`](@ref)); with `free`, `fixed`, `priors`, ...
it is also minimised over oscillation parameters ([`marginalize`](@ref)),
with the scanned parameters always fixed. The matrix has `length(xs)` rows.

```julia
grid = chi2_scan(truth, :sin²2θ13 => range(0.07, 0.11, 41), :δ => range(0, 2π, 61))
proj = chi2_scan(truth, :δ => range(0, 2π, 61); fixed=(:θ12, :Δm21), priors=(ρ=0.05,))
```
"""
function chi2_scan(base::OscParams, x::Pair{Symbol}; kwargs...)
    return [scan_point(with_param(base, x.first, v), (x.first,); kwargs...) for v in x.second]
end

function chi2_scan(base::OscParams, x::Pair{Symbol}, y::Pair{Symbol}; kwargs...)
    return [
        scan_point(
            with_param(with_param(base, x.first, u), y.first, v),
            (x.first, y.first);
            kwargs...,
        ) for u in x.second, v in y.second
    ]
end

function scan_point(
    p::OscParams,
    scanned;
    free = nothing,
    fixed = nothing,
    experiment = :all,
    rule = :all,
    systematics = true,
    setup = current_setup(),
    fit_kwargs...,
)
    if free === nothing && fixed === nothing && isempty(fit_kwargs)
        return chi2(p; experiment, rule, systematics, setup)
    end
    angles = Tuple(
        haskey(DERIVED, canonical(s)) ? DERIVED[canonical(s)][1] : canonical(s) for s in scanned
    )
    if free !== nothing
        free = Tuple(setdiff(canonical.(collect(free)), angles))
    else
        fixed = Tuple(union(canonical.(collect(something(fixed, ()))), angles))
    end
    return marginalize(p; free, fixed, experiment, setup, fit_kwargs...).chi2
end

# ---------------------------------------------------------------------------
#  Systematics
# ---------------------------------------------------------------------------

"""
    systematics!(e, on::Bool; rule=:all)
    systematics!(on::Bool; setup=current_setup())

Switch the systematics of rules on or off. With systematics off, a rule uses
its `@sys_off_function` (no nuisance parameters, by default).
"""
function systematics!(e::Experiment, on_::Bool; rule = :all)
    on(
        i -> glbSwitchSystematics(
            i,
            rule === :all ? GLB_ALL : rule_index(e, rule),
            on_ ? GLB_ON : GLB_OFF,
        ),
        e,
    )
    return e
end
systematics!(on_::Bool; setup::Setup = current_setup()) =
    (in_setup(() -> glbSwitchSystematics(GLB_ALL, GLB_ALL, on_ ? GLB_ON : GLB_OFF), setup); nothing)

"""
    chi2_function(e; rule=1, systematics=true) -> String
    chi2_function!(e, name; rule=:all, errors, systematics=true)

The χ² function a rule uses with systematics on (or off), by name, and
changing it. `errors` are the widths of the priors of its nuisance
parameters. Built in: `"chiSpectrumTilt"`, `"chiSpectrumCalib"`,
`"chiSpectrumOnly"`, `"chiTotalRatesTilt"`, `"chiNoSysSpectrum"`,
`"chiNoSysTotalRates"`, `"chiZero"`; others with [`define_chi2`](@ref).
"""
chi2_function(e::Experiment; rule = 1, systematics::Bool = true) =
    on(i -> glbGetChiFunction(i, rule_index(e, rule), systematics ? GLB_ON : GLB_OFF), e)

function chi2_function!(
    e::Experiment,
    name::AbstractString;
    rule = :all,
    errors = Float64[],
    systematics::Bool = true,
)
    on(e) do i
        st = glbSetChiFunction(
            i,
            rule === :all ? GLB_ALL : rule_index(e, rule),
            systematics ? GLB_ON : GLB_OFF,
            String(name),
            Float64.(errors),
        )
        st == 0 || throw(GLoBESError("could not set χ² function $name"))
    end
    return e
end
@doc (@doc chi2_function) chi2_function!

"""
    systematic_errors(e; rule=1, systematics=true) -> Vector{Float64}
    systematic_errors!(e, errors; rule=:all, systematics=true)

Prior widths of the nuisance parameters of a rule's χ² function.
"""
systematic_errors(e::Experiment; rule = 1, systematics::Bool = true) =
    on(i -> copy(glbGetSysErrorsListPtr(i, rule_index(e, rule), systematics ? GLB_ON : GLB_OFF)), e)
function systematic_errors!(e::Experiment, errors; rule = :all, systematics::Bool = true)
    on(
        i -> glbSetSysErrorsList(
            i,
            rule === :all ? GLB_ALL : rule_index(e, rule),
            systematics ? GLB_ON : GLB_OFF,
            Float64.(errors),
        ),
        e,
    )
    return e
end

"""
    define_chi2(f, name, n_params; setup=current_setup())

Define a χ² function for use in AEDL files (`@sys_on_function = "name"`) or
with [`chi2_function!`](@ref). It is called as `f(e::Experiment, rule::Int, x,
errors)` with the nuisance parameters `x` (minimised over) and their prior
widths `errors`, and returns χ² including the priors of `x`. Inside, use
[`true_rates`](@ref) and [`fit_rates`](@ref), and [`poisson_chi2`](@ref).

```julia
define_chi2("normOnly", 1) do e, rule, x, errors
    data = true_rates(e; rule)
    s, b = fit_rates(e; rule)
    sum(poisson_chi2(data[i], (1 + x[1]) * (s[i] + b[i])) for i in energy_window_bins(e; rule)) +
        (x[1] / errors[1])^2
end
```
"""
function define_chi2(f, name::AbstractString, n_params::Integer; setup::Setup = current_setup())
    wrapper =
        (exp, rule, n, x, errors, ud) ->
            Float64(f(Experiment(setup, setup.state.experiments[exp]), rule, x, errors))
    in_setup(() -> glbDefineChiFunction(wrapper, n_params, String(name)), setup)
    return nothing
end

"""
    poisson_chi2(observed, predicted) -> Float64

The Poisson log-likelihood χ², `2(predicted - observed + observed log(observed/predicted))`,
of one bin, with the treatment of empty bins and vanishing predictions
described in the manual.
"""
poisson_chi2(observed::Real, predicted::Real) =
    glb_likelihood(Float64(observed), Float64(predicted))

"""
    energy_scale_shift(rates, b; centers, widths) -> Vector{Float64}

Rates of a binned spectrum after an energy calibration error `b` (the true
energy is `(1+b)` times the reconstructed one). See the manual for the
method.
"""
energy_scale_shift(
    rates::AbstractVector,
    b::Real;
    centers::AbstractVector,
    widths::AbstractVector,
) = shift_energy_scale(Float64(b), rates, centers, widths)
energy_scale_shift(e::Experiment, rates::AbstractVector, b::Real) = shift_energy_scale(
    Float64(b),
    rates,
    e.data.smear_data[1].bincenter,
    e.data.smear_data[1].binsize,
)

# ---------------------------------------------------------------------------
#  Non-standard physics
# ---------------------------------------------------------------------------

"""
    probability_engine!(P!, extra; setup=current_setup())
    standard_engine!(; setup=current_setup())

Replace the oscillation probabilities with your own, which may depend on
additional parameters `extra` (names, e.g. `(:σE,)`). `P!` is called as
`P!(P, params::OscParams, E, antineutrino, lengths, densities)` and fills
the 3×3 matrix `P[from, to]` for energy `E` in GeV along a matter profile
of layers with `lengths` in km and `densities` in g/cm³.

Experiments loaded before the change keep working. Create `OscParams` with
`extra=(σE=0.0,)` to set the additional parameters.
"""
function probability_engine!(P!, extra = (); setup::Setup = current_setup())
    names = collect(Symbol, extra)
    current = Ref{OscParams}()
    template = OscParams(zeros(6 + length(names)), names, [1.0])
    set = (p, ud) -> (current[] = OscParams(copy(p.osc), names, [1.0]); 0)
    get = (p, ud) -> (isassigned(current) && (p.osc .= current[].values); 0)
    prob =
        (P, cp, E, psteps, l, d, fs, ud) ->
            (P!(P, current[], E, cp < 0, view(l, 1:psteps), view(d, 1:psteps)); 0)
    in_setup(setup) do
        glbRegisterProbabilityEngine(length(template.values), prob, set, get)
        state().extra_names = names
        resize_params!()
    end
    return nothing
end

function standard_engine!(; setup::Setup = current_setup())
    in_setup(setup) do
        glbRegisterProbabilityEngine(6, nothing, nothing, nothing)
        state().extra_names = Symbol[]
        resize_params!()
    end
    return nothing
end

# Simulated data of the old engine no longer fit the parameter vectors
resize_params!() = (state().truth = nothing)

# ---------------------------------------------------------------------------
#  Fluxes and cross sections
# ---------------------------------------------------------------------------

"""
    flux(e, E; flux=1, flavour=:e, antineutrino=false, L=baseline(e)) -> Float64

Neutrino flux of an experiment at energy `E` in GeV and distance `L` in km,
in the units of the flux tables (see the manual).
"""
flux(
    e::Experiment,
    E::Real;
    flux::Integer = 1,
    flavour = :e,
    antineutrino::Bool = false,
    L::Real = baseline(e),
) = on(
    i -> glbFlux(i, flux, Float64(E), Float64(L), GLoBES.flavour(flavour), antineutrino ? -1 : +1),
    e,
)

"""
    cross_section(e, E; xsec=1, flavour=:e, antineutrino=false) -> Float64

Total cross section per target nucleon, in 10⁻³⁸ cm², of the cross section
table `xsec` of an experiment at energy `E` in GeV. (The tables hold the cross
section divided by energy, σ/E.)
"""
cross_section(e::Experiment, E::Real; xsec::Integer = 1, flavour = :e, antineutrino::Bool = false) =
    on(i -> glbXSection(i, xsec, Float64(E), GLoBES.flavour(flavour), antineutrino ? -1 : +1), e)
