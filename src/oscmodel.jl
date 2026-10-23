# Oscillation models defined in files written like SARAH model files: each
# ingredient of the effective Hamiltonian is one assignment, and every
# ingredient that is left out takes its built-in default.
#
#   MixingMatrix = R23 . Ud . R13 . Conj[Ud] . R12;
#   MassSquared = Diag[{0, dm21, dm31}];
#   MassHierarchy[[1]] = {NH, {m1 -> mlightest, m2 -> Sqrt[m1^2 + dm21], m3 -> Sqrt[m1^2 + dm31]}};
#   MassMatrix = U . Diag[{m1, m2, m3}] . Dagger[U] + Sqrt[Abs[dm31]] eta;
#   AdditionalTerms[[1]] = {NSI, A eps};
#   Decoherence = Gamma0 (E/GeV)^n;
#
# H = M M†/2E (or U MassSquared U†/2E without a mass matrix) + Σ AdditionalTerms.
# All quantities are in eV: E is the neutrino energy, A = √2 G_F N_e the
# matter potential, L the length of a layer (in eV⁻¹); `GeV`, `km`, ... are
# units. For antineutrinos H → H* with A → -A and CP → -1.

# The defaults, written in the model language itself
const DEFAULT_MODEL_SOURCE = """
R12 = Rotation[1, 2, theta12];
R13 = Rotation[1, 3, theta13];
R23 = Rotation[2, 3, theta23];
Ud = PhaseMatrix[3, deltaCP];
PMNS = R23 . Ud . R13 . Conj[Ud] . R12;

MixingMatrix = PMNS;
MassSquared = Diag[{0, dm21, dm31}];
MassHierarchy[[1]] = {NH, {m1 -> mlightest, m2 -> Sqrt[m1^2 + dm21], m3 -> Sqrt[m1^2 + dm31]}};
MassHierarchy[[2]] = {IH, {m3 -> mlightest, m1 -> Sqrt[m3^2 - dm31], m2 -> Sqrt[m1^2 + dm21]}};
AdditionalTerms[[1]] = {Matter, A Diag[{1, 0, 0}]};
"""

# Ingredients with one definition (`MixingMatrix = ...;`) and the symbol that
# holds their value in other definitions
const SINGLE_BLOCKS = (
    (:MixingMatrix, :U),
    (:MassSquared, :M2),
    (:MassMatrix, :M),
    (:Decoherence, :Gamma),
    (:Hamiltonian, :H),
)
# Ingredients with named entries (`AdditionalTerms[[i]] = {name, definition};`)
const LIST_BLOCKS = (:MassHierarchy, :AdditionalTerms)

function main_symbol(kind)
    i = findfirst(k -> k[1] === kind, SINGLE_BLOCKS)
    return i === nothing ? (kind === :AdditionalTerms ? :H : nothing) : SINGLE_BLOCKS[i][2]
end

"""
    ModelBlock

An ingredient of a model file: its `kind` (`:MixingMatrix`, `:MassSquared`,
`:MassMatrix`, `:MassHierarchy`, `:AdditionalTerms`, `:Decoherence` or
`:Hamiltonian`), `name` (the first entry of `MassHierarchy[[i]]` and
`AdditionalTerms[[i]]`, the kind otherwise), the expression that defines it
and, for a mass hierarchy, the masses it defines (`locals`, in `order`).
"""
struct ModelBlock
    kind::Symbol
    name::Symbol
    definition::Any
    locals::Dict{Symbol, Any}
    order::Vector{Symbol}
    off::Bool
    default::Bool                  # one of the built-in blocks
end

block_label(b::ModelBlock) = b.kind in LIST_BLOCKS ? "$(b.kind)[$(b.name)]" : String(b.kind)

"""
    ModelParameter

A parameter of an oscillation model: its `name` in the model files, the
name in `OscParams` (`key`), its default value (in the model's units, eV),
LaTeX name, display unit and scan range, and where the default comes from
(`:nufit`, `:file` or `:auto`).
"""
mutable struct ModelParameter
    name::Symbol
    key::Symbol
    value::Float64
    latex::String
    unit::String
    range::Union{Nothing, Tuple{Float64, Float64}}
    description::String
    source::Symbol
end

"""
    OscillationModel

A neutrino oscillation model read from model files with
[`load_model`](@ref): flavours, parameters, blocks and settings. Use it as
the probability engine with [`use_model!`](@ref), make parameters with
`OscParams(model; ...)` and compute [`hamiltonian`](@ref),
[`effective_mixing`](@ref) and [`effective_masses`](@ref).

```julia
model = load_model("examples/models/SNSI")   # model.flavours, model.parameters, model.blocks
```
"""
mutable struct OscillationModel
    name::String
    latex_name::String
    authors::String
    date::String
    flavours::Vector{Symbol}
    globals::Dict{Symbol, Any}
    global_order::Vector{Symbol}
    blocks::Vector{ModelBlock}
    parameters::Vector{ModelParameter}
    hierarchy::Symbol                       # a MassHierarchy name, or :Auto
    nufit::NuFit
    dir::String
    files::Dict{Symbol, String}             # :model, :parameters, :particles, :run
    warnings::Vector{String}
end

nflavours(m::OscillationModel) = length(m.flavours)

function Base.show(io::IO, m::OscillationModel)
    print(
        io,
        "OscillationModel(\"",
        m.name,
        "\", ",
        nflavours(m),
        " flavours, ",
        length(m.parameters),
        " parameters)",
    )
end

function Base.show(io::IO, ::MIME"text/plain", m::OscillationModel)
    println(io, "OscillationModel \"", m.name, "\"", isempty(m.authors) ? "" : " by $(m.authors)")
    println(io, "  flavours:   ", join(m.flavours, ", "))
    for b in m.blocks
        b.off && continue
        println(io, "  ", rpad(block_label(b), 28), b.default ? "(default)" : definition_string(b))
    end
    println(io, "  hierarchy:  ", m.hierarchy)
    print(io, "  parameters: ", join((p.name for p in m.parameters), ", "))
end

function definition_string(b::ModelBlock)
    b.definition === nothing && return "{" * join(b.order, ", ") * "}"
    s = sprint(show_mexpr, b.definition)
    length(s) > 50 ? s[1:47] * "..." : s
end

# ---------------------------------------------------------------------------
#  Reading model files
# ---------------------------------------------------------------------------

const STANDARD_PARAMETERS = (
    (:theta12, :θ12, "\\theta_{12}", "°"),
    (:theta13, :θ13, "\\theta_{13}", "°"),
    (:theta23, :θ23, "\\theta_{23}", "°"),
    (:deltaCP, :δ, "\\delta_{\\mathrm{CP}}", "°"),
    (:dm21, :Δm21, "\\Delta m^2_{21}", "eV²"),
    (:dm31, :Δm31, "\\Delta m^2_{31}", "eV²"),
    (:mlightest, :m_lightest, "m_{\\mathrm{lightest}}", "eV"),
)
const STANDARD_ALIASES = Dict(:delta => :deltaCP, :dCP => :deltaCP, :delta_CP => :deltaCP)

# Symbols that are not parameters
const MODEL_CONSTANTS = Dict(
    :eV => 1.0,
    :keV => 1e3,
    :MeV => 1e6,
    :GeV => 1e9,
    :TeV => 1e12,
    :km => 1 / GLB_EV_TO_KM_FACTOR,
    :m => 1e-3 / GLB_EV_TO_KM_FACTOR,
    :cm => 1e-5 / GLB_EV_TO_KM_FACTOR,
    :Pi => π,
    :Degree => π / 180,
    :GF => G_FERMI * 1e-18,
)
const DYNAMIC_SYMBOLS = (:E, :A, :L, :CP, :rho)
const MODEL_FUNCTIONS = Set([
    :Plus,
    :Times,
    :Dot,
    :Power,
    :List,
    :Sqrt,
    :Exp,
    :MatrixExp,
    :Log,
    :Sin,
    :Cos,
    :Tan,
    :ArcSin,
    :ArcCos,
    :ArcTan,
    :Abs,
    :Arg,
    :Re,
    :Im,
    :Conj,
    :Conjugate,
    :Dagger,
    :ConjugateTranspose,
    :HC,
    :Transpose,
    :Inverse,
    :Det,
    :Tr,
    :Diag,
    :DiagonalMatrix,
    :IdentityMatrix,
    :Rotation,
    :Sign,
    :Eigenvalues,
    :Eigenvectors,
    :Max,
    :Min,
    :Range,
    :Rule,
    :Set,
    :Off,
    :On,
    :Block,
    :ZeroMatrix,
    :MatrixPower,
    :Commutator,
    :PhaseMatrix,
])

"""
    load_model(dir; model, parameters, particles, run) -> OscillationModel

Read an oscillation model from a directory like a SARAH model:

    MyModel/
      MyModel.m        the model: mixing matrix, masses, hierarchies, additional terms
      parameters.m     parameter values (NuFIT by default), LaTeX names, units, ranges
      particles.m      optional: Flavours = {e, mu, tau, s}
      Run.m            what to calculate, plot, fit and write to LaTeX (see run_model)

Any ingredient left out takes the default: the PDG mixing matrix, Δm² from
NuFIT, the NH and IH hierarchies and the standard matter potential. An
empty directory, or `load_model()`, is the standard three-flavour model.

```julia
model = load_model("examples/models/SNSI")   # directory with SNSI.m, parameters.m, particles.m, Run.m
model = load_model("examples/models/SNSI/SNSI.m"; parameters = "examples/models/SNSI/parameters.m")   # the model file, the others by name
standard = load_model()   # the standard three-flavour model
```
"""
function load_model(
    dir::AbstractString = "";
    model = nothing,
    parameters = nothing,
    particles = nothing,
    run = nothing,
)
    dir = isempty(dir) ? "" : expanduser(dir)
    files = Dict{Symbol, String}()
    if !isempty(dir) && isfile(dir)
        files[:model] = dir
        dir = dirname(dir)
    elseif !isempty(dir)
        isdir(dir) || throw(ArgumentError("no model directory $dir"))
        main = joinpath(dir, basename(rstrip(dir, '/')) * ".m")
        others = filter(readdir(dir)) do f
            endswith(f, ".m") &&
                !(lowercase(f) in ("parameters.m", "particles.m", "run.m", "main.m"))
        end
        if isfile(main)
            files[:model] = main
        elseif length(others) == 1
            files[:model] = joinpath(dir, only(others))
        elseif length(others) > 1
            throw(
                ArgumentError("cannot tell which file in $dir is the model: $(join(others, ", "))"),
            )
        end
        for (k, names) in (
            (:parameters, ("parameters.m",)),
            (:particles, ("particles.m",)),
            (:run, ("Run.m", "run.m", "Main.m", "main.m")),
        )
            for n in names
                isfile(joinpath(dir, n)) && (files[k] = joinpath(dir, n); break)
            end
        end
    end
    for (k, v) in ((:model, model), (:parameters, parameters), (:particles, particles), (:run, run))
        v === nothing || (files[k] = expanduser(v))
    end
    m = new_model(dir, files)
    sources = [(k, read(files[k], String)) for k in (:particles, :model) if haskey(files, k)]
    for (_, src) in sources
        read_model_statements!(m, parse_mathematica(src))
    end
    add_default_blocks!(m)
    params_src =
        haskey(files, :parameters) ? parse_mathematica(read(files[:parameters], String)) : Any[]
    setup_parameters!(m, params_src)
    validate_model!(m)
    return m
end

new_model(dir = "", files = Dict{Symbol, String}()) = OscillationModel(
    "Standard",
    "",
    "",
    "",
    [:e, :mu, :tau],
    Dict{Symbol, Any}(),
    Symbol[],
    ModelBlock[],
    ModelParameter[],
    :Auto,
    nufit(),
    dir,
    files,
    String[],
)

"""
    model_from_string(src; parameters="") -> OscillationModel

A model from source text instead of files (`src` as in the model file,
`parameters` as in `parameters.m`).

```julia
liv = model_from_string("AdditionalTerms[[2]] = {LIV, c E Diag[{1, -1, 0}]};"; parameters = "ParameterDefinitions = {{c, {Value -> 0}}};")   # model file and parameters.m as text
```
"""
function model_from_string(src::AbstractString; parameters::AbstractString = "")
    m = new_model()
    m.name = "Model"
    read_model_statements!(m, parse_mathematica(src))
    add_default_blocks!(m)
    setup_parameters!(m, parse_mathematica(parameters))
    validate_model!(m)
    return m
end

function read_model_statements!(m::OscillationModel, statements; default = false)
    for st in statements
        if !(is_call(st, :Set) || is_call(st, :SetDelayed))
            push!(m.warnings, "ignored statement $(sprint(show_mexpr, st))")
            continue
        end
        lhs, rhs = st.args
        if lhs isa Symbol && any(k -> k[1] === lhs, SINGLE_BLOCKS)
            add_block!(
                m,
                ModelBlock(lhs, lhs, rhs, Dict{Symbol, Any}(), Symbol[], rhs === :Off, default),
            )
        elseif is_call(lhs, :Part) && lhs.args[1] in LIST_BLOCKS
            push!(m.blocks, list_block(lhs.args[1], rhs, default))
            replace_named!(m)
        elseif lhs isa Symbol
            name = String(lhs)
            if startswith(name, "Model`")
                key = name[7:end]
                key == "Name" && (m.name = rhs)
                key == "NameLaTeX" && (m.latex_name = rhs)
                key == "Authors" && (m.authors = rhs)
                key == "Date" && (m.date = rhs)
            elseif lhs === :Flavours || lhs === :Flavors
                m.flavours = Symbol[f isa Symbol ? f : Symbol(f) for f in rhs.args]
            elseif lhs === :Hierarchy || lhs === :MassOrdering
                m.hierarchy = rhs
            else
                default && haskey(m.globals, lhs) && continue
                m.globals[lhs] = rhs
                lhs in m.global_order || push!(m.global_order, lhs)
            end
        else
            throw(
                ArgumentError(
                    "cannot assign to $(sprint(show_mexpr, lhs)); the ingredients are " *
                    join((String(k[1]) * " = ..." for k in SINGLE_BLOCKS), ", ") *
                    ", MassHierarchy[[i]] = {name, {m1 -> ..., ...}} and AdditionalTerms[[i]] = {name, H}",
                ),
            )
        end
    end
end

# A model sets each single ingredient once: a later definition replaces an earlier one
function add_block!(m::OscillationModel, b::ModelBlock)
    filter!(x -> x.kind !== b.kind, m.blocks)
    push!(m.blocks, b)
end

# MassHierarchy[[i]] = {name, {m1 -> ..., ...}} and AdditionalTerms[[i]] = {name, H}
function list_block(kind::Symbol, rhs, default::Bool)
    is_call(rhs, :List) && length(rhs.args) == 2 && rhs.args[1] isa Symbol || throw(
        ArgumentError(
            kind === :MassHierarchy ?
            "write MassHierarchy[[i]] = {name, {m1 -> ..., m2 -> ..., m3 -> ...}}" :
            "write AdditionalTerms[[i]] = {name, H}",
        ),
    )
    name, def = rhs.args
    kind === :AdditionalTerms &&
        return ModelBlock(kind, name, def, Dict{Symbol, Any}(), Symbol[], def === :Off, default)
    is_call(def, :List) && all(r -> is_call(r, :Rule) && r.args[1] isa Symbol, def.args) ||
        throw(ArgumentError("the masses of MassHierarchy $name are rules, {m1 -> ..., m2 -> ...}"))
    locals = Dict{Symbol, Any}(r.args[1] => r.args[2] for r in def.args)
    return ModelBlock(
        kind,
        name,
        nothing,
        locals,
        Symbol[r.args[1] for r in def.args],
        false,
        default,
    )
end

# A list entry replaces an earlier one of the same kind and name
function replace_named!(m::OscillationModel)
    b = m.blocks[end]
    i = findfirst(x -> x.kind === b.kind && x.name === b.name, m.blocks)
    i < length(m.blocks) && (m.blocks[i] = pop!(m.blocks))
end

# Built-in blocks that the model does not replace
function add_default_blocks!(m::OscillationModel)
    defaults = new_model()
    read_model_statements!(defaults, parse_mathematica(DEFAULT_MODEL_SOURCE); default = true)
    for k in defaults.global_order
        haskey(m.globals, k) && continue
        m.globals[k] = defaults.globals[k]
        push!(m.global_order, k)
    end
    has(kind) = any(b -> b.kind === kind, m.blocks)
    for b in defaults.blocks
        if b.kind === :MassSquared
            (has(:MassSquared) || has(:MassMatrix) || has(:Hamiltonian)) && continue
        elseif b.kind in LIST_BLOCKS
            any(x -> x.kind === b.kind && x.name === b.name, m.blocks) && continue
            b.kind === :AdditionalTerms && has(:Hamiltonian) && continue
        else
            has(b.kind) && continue
        end
        push!(m.blocks, b)
    end
    return m
end

# ---------------------------------------------------------------------------
#  Parameters
# ---------------------------------------------------------------------------

function setup_parameters!(m::OscillationModel, statements)
    defs = Dict{Symbol, Dict{Symbol, Any}}()
    order = Symbol[]
    nf_opts = Dict{Symbol, Any}()
    for st in statements
        is_call(st, :Set) || continue
        lhs, rhs = st.args
        if lhs === :ParameterDefinitions
            for entry in rhs.args
                name = canonical_model_name(entry.args[1])
                defs[name] = options(entry.args[2])
                push!(order, name)
            end
        elseif lhs === :NuFit
            nf_opts = options(rhs)
        elseif lhs === :Hierarchy || lhs === :MassOrdering
            m.hierarchy = rhs
        elseif lhs isa Symbol
            # name = value
            name = canonical_model_name(lhs)
            defs[name] = Dict{Symbol, Any}(:Value => rhs)
            push!(order, name)
        end
    end
    version = String(get(nf_opts, :Version, NUFIT_LATEST))
    sk = get(nf_opts, :SK, :False) === :True
    ordering = get(nf_opts, :Ordering, nothing)
    ordering === nothing && (ordering = m.hierarchy === :IH ? :inverted : :normal)
    ordering =
        ordering in (:Inverted, :IH, :IO, :inverted) ? :inverted :
        ordering in (:Best, :best) ? :best : :normal
    m.nufit = nufit(version; ordering, sk)
    best = nufit_params(m.nufit)
    nufit_value = Dict(
        :theta12 => best.θ12,
        :theta13 => best.θ13,
        :theta23 => best.θ23,
        :deltaCP => best.δ,
        :dm21 => best.Δm21,
        :dm31 => best.Δm31,
        :mlightest => 0.0,
    )
    params = ModelParameter[]
    for (name, key, latex, unit) in STANDARD_PARAMETERS
        push!(
            params,
            ModelParameter(
                name,
                key,
                nufit_value[name],
                latex,
                unit,
                nothing,
                "",
                name === :mlightest ? :default : :nufit,
            ),
        )
    end
    used = free_symbols(m)
    for name in union(order, sort!(collect(setdiff(used, order))))
        any(p -> p.name === name, params) && continue
        push!(
            params,
            ModelParameter(name, canonical(name), 0.0, guess_latex(name), "", nothing, "", :auto),
        )
        haskey(defs, name) ||
            push!(m.warnings, "parameter $name is not defined in parameters.m; it is 0 by default")
    end
    for (name, opts) in defs
        i = findfirst(p -> p.name === name, params)
        i === nothing && (
            push!(m.warnings, "parameter $name in parameters.m is not used by the model");
            continue
        )
        p = params[i]
        if haskey(opts, :Value) && opts[:Value] !== :NuFit
            p.value =
                Float64(real(model_eval(opts[:Value], ModelContext(m, Dict{Symbol, Float64}()), 0)))
            p.source = :file
        end
        haskey(opts, :LaTeX) && (p.latex = opts[:LaTeX])
        haskey(opts, :Unit) && (p.unit = string(opts[:Unit]))
        haskey(opts, :Description) && (p.description = opts[:Description])
        if haskey(opts, :Range)
            r = model_eval(opts[:Range], ModelContext(m, Dict{Symbol, Float64}()), 0)
            p.range = (Float64(real(r[1])), Float64(real(r[2])))
        end
    end
    m.parameters = params
    return m
end

canonical_model_name(s::Symbol) = get(STANDARD_ALIASES, s, s)

# A LaTeX name from an ASCII parameter name: eta_etau → \eta_{e\tau}
function guess_latex(name::Symbol)
    greek = Dict(
        "alpha" => "\\alpha",
        "beta" => "\\beta",
        "gamma" => "\\gamma",
        "Gamma" => "\\Gamma",
        "delta" => "\\delta",
        "Delta" => "\\Delta",
        "epsilon" => "\\epsilon",
        "eps" => "\\epsilon",
        "eta" => "\\eta",
        "theta" => "\\theta",
        "lambda" => "\\lambda",
        "mu" => "\\mu",
        "nu" => "\\nu",
        "phi" => "\\phi",
        "sigma" => "\\sigma",
        "tau" => "\\tau",
        "omega" => "\\omega",
        "rho" => "\\rho",
    )
    tr(s) = replace(
        s,
        r"(Gamma|Delta|alpha|beta|gamma|delta|epsilon|eps|eta|theta|lambda|mu|nu|phi|sigma|tau|omega|rho)" =>
            g -> greek[g] * " ",
    )
    s = String(name)
    parts = split(s, '_'; limit = 2)
    base = parts[1]
    m = match(r"^([A-Za-z]+?)(\d+)$", base)
    m !== nothing &&
        length(parts) == 1 &&
        return strip(tr(m.captures[1])) * "_{" * m.captures[2] * "}"
    length(parts) == 1 && return strip(tr(base))
    return strip(tr(base)) * "_{" * subscript_latex(parts[2]) * "}"
end

# Flavour subscripts: etau → e\tau, mumu → \mu\mu; otherwise Greek names
function subscript_latex(sub)
    flavours = Dict("e" => "e", "mu" => "\\mu", "tau" => "\\tau", "s" => "s")
    if occursin(r"^(e|mu|tau|s|\d)+$", sub)
        return join(get(flavours, t.match, t.match) for t in eachmatch(r"mu|tau|e|s|\d", sub))
    end
    return replace(
        strip(
            replace(
                sub,
                r"(Gamma|Delta|alpha|beta|gamma|delta|epsilon|eps|eta|theta|lambda|mu|nu|phi|sigma|tau|omega|rho)" =>
                    g -> "\\" * (g == "eps" ? "epsilon" : g) * " ",
            ),
        ),
        " " => "",
    )
end

# Symbols that the model uses but does not define
function free_symbols(m::OscillationModel)
    defined = Set{Symbol}(keys(m.globals))
    for b in m.blocks
        union!(defined, keys(b.locals))
        b.kind === :AdditionalTerms && push!(defined, b.name)
    end
    union!(
        defined,
        (:U, :M2, :M, :H, :Gamma, :MixingMatrix, :MassSquared, :MassMatrix, :Hvac, :I, :N),
    )
    union!(
        defined,
        keys(MODEL_CONSTANTS),
        DYNAMIC_SYMBOLS,
        (first(p) for p in STANDARD_PARAMETERS),
        keys(STANDARD_ALIASES),
    )
    found = Set{Symbol}()
    walk(e) =
        if e isa Symbol
            e in defined || push!(found, e)
        elseif e isa MExpr
            e.head isa Symbol && e.head in MODEL_FUNCTIONS || walk(e.head)
            foreach(walk, e.args)
        end
    foreach(walk, values(m.globals))
    for b in m.blocks
        b.off && continue
        b.definition === nothing || walk(b.definition)
        foreach(walk, values(b.locals))
    end
    return Set(canonical_model_name(s) for s in found)
end

function validate_model!(m::OscillationModel)
    m.hierarchy === :Auto ||
        any(b -> b.kind === :MassHierarchy && b.name === m.hierarchy, m.blocks) ||
        throw(ArgumentError("no MassHierarchy named $(m.hierarchy)"))
    return m
end

"""
    model_parameters(model) -> Vector{ModelParameter}

The parameters of a model with their defaults.

```julia
model_parameters(model)   # name, key in OscParams, default value, LaTeX name, unit and range of each parameter
```
"""
model_parameters(m::OscillationModel) = m.parameters

# ---------------------------------------------------------------------------
#  Evaluation
# ---------------------------------------------------------------------------

struct DependsOnEnergy <: Exception end

# Values of model symbols for one parameter point and, optionally, energy and layer
mutable struct ModelContext
    model::OscillationModel
    values::Dict{Symbol, Float64}          # by model name
    dynamic::Dict{Symbol, Any}             # E, A, L, CP, rho (absent: static evaluation)
    cache::Dict{Tuple{Int, Symbol}, Any}
    parent::Union{Nothing, ModelContext}   # static context whose cache can be used
    active::Set{Tuple{Int, Symbol}}        # for detecting circular definitions
end

ModelContext(m::OscillationModel, values, dynamic = Dict{Symbol, Any}(); parent = nothing) =
    ModelContext(
        m,
        values,
        dynamic,
        Dict{Tuple{Int, Symbol}, Any}(),
        parent,
        Set{Tuple{Int, Symbol}}(),
    )

function selected_hierarchy(ctx::ModelContext)
    m = ctx.model
    h = m.hierarchy
    if h === :Auto
        haskey(ctx.values, :dm31) || return nothing
        h = ctx.values[:dm31] >= 0 ? :NH : :IH
    end
    i = findfirst(b -> b.kind === :MassHierarchy && b.name === h, m.blocks)
    return i
end

# Scope 0 is the model's globals; scope i > 0 the locals of block i
function resolve(ctx::ModelContext, s::Symbol, scope::Int)
    m = ctx.model
    s = canonical_model_name(s)
    if haskey(ctx.dynamic, s)
        return ctx.dynamic[s]
    elseif s in DYNAMIC_SYMBOLS
        throw(DependsOnEnergy())
    end
    scopes = scope == 0 ? Int[] : [scope]
    h = selected_hierarchy(ctx)
    h === nothing || h in scopes || push!(scopes, h)
    for sc in scopes
        b = m.blocks[sc]
        haskey(b.locals, s) && return cached(ctx, (sc, s)) do
            model_eval(b.locals[s], ctx, sc)
        end
    end
    # Globals are evaluated where they are used, so that a block can
    # redefine their ingredients (R12 = ... changes PMNS in that block)
    haskey(m.globals, s) && return cached(ctx, (scope, s)) do
        model_eval(m.globals[s], ctx, scope)
    end
    haskey(ctx.values, s) && return ctx.values[s]
    s === :I && return im
    s === :N && return nflavours(m)
    haskey(MODEL_CONSTANTS, s) && return MODEL_CONSTANTS[s]
    s in (:U, :MixingMatrix) && return block_value(ctx, :MixingMatrix)
    s in (:M2, :MassSquared) && return block_value(ctx, :MassSquared)
    s in (:M, :MassMatrix) && return block_value(ctx, :MassMatrix)
    s === :Hvac && return vacuum_hamiltonian(ctx)
    i = findfirst(b -> b.kind === :AdditionalTerms && b.name === s, m.blocks)
    i === nothing || return block_value(ctx, m.blocks[i])
    throw(ArgumentError("$s is not defined in model $(m.name)"))
end

function cached(f, ctx::ModelContext, key)
    haskey(ctx.cache, key) && return ctx.cache[key]
    ctx.parent !== nothing && haskey(ctx.parent.cache, key) && return ctx.parent.cache[key]
    key in ctx.active && throw(ArgumentError("circular definition of $(key[2])"))
    push!(ctx.active, key)
    try
        v = f()
        ctx.cache[key] = v
        return v
    finally
        delete!(ctx.active, key)
    end
end

function block_value(ctx::ModelContext, kind::Symbol)
    m = ctx.model
    i = findfirst(b -> b.kind === kind && !b.off, m.blocks)
    i === nothing && throw(ArgumentError("model $(m.name) has no $kind"))
    return block_value(ctx, m.blocks[i])
end

function block_value(ctx::ModelContext, b::ModelBlock)
    i = findfirst(x -> x === b, ctx.model.blocks)
    return cached(ctx, (i, Symbol("#value"))) do
        v = model_eval(b.definition, ctx, i)
        b.kind === :MassSquared && v isa AbstractVector ? diagm(ComplexF64.(v)) : v
    end
end

has_block(m::OscillationModel, kind) = any(b -> b.kind === kind && !b.off, m.blocks)

# (M M† or U M2 U†) in eV², independent of energy
function vacuum_mass_squared(ctx::ModelContext)
    m = ctx.model
    if has_block(m, :MassMatrix)
        M = as_matrix(block_value(ctx, :MassMatrix), m)
        return M * M'
    end
    U = as_matrix(block_value(ctx, :MixingMatrix), m)
    M2 = as_matrix(block_value(ctx, :MassSquared), m)
    return U * M2 * U'
end

vacuum_hamiltonian(ctx::ModelContext) = vacuum_mass_squared(ctx) / (2 * dynamic(ctx, :E))

dynamic(ctx, s) = haskey(ctx.dynamic, s) ? ctx.dynamic[s] : throw(DependsOnEnergy())

function as_matrix(x, m::OscillationModel)
    n = nflavours(m)
    x isa Number && return Matrix{ComplexF64}(x * I, n, n)
    size(x) == (n, n) ||
        throw(DimensionMismatch("a matrix of model $(m.name) is $(size(x)), not $n×$n"))
    return ComplexF64.(x)
end

# Numerical value of a model expression
function model_eval(e, ctx::ModelContext, scope::Int)
    e isa MNumber && return e
    e isa String && return e
    e isa Symbol && return resolve(ctx, e, scope)
    e isa MExpr || return e
    h = e.head
    ev(x) = model_eval(x, ctx, scope)
    if h === :List
        items = map(ev, e.args)
        if !isempty(items) && all(x -> x isa AbstractVector, items) && allequal(length.(items))
            return permutedims(reduce(hcat, items))
        end
        return [x for x in items]
    end
    h === :Plus && return reduce(model_plus, map(ev, e.args))
    h === :Times && return reduce(model_times, map(ev, e.args))
    h === :Dot && return reduce(model_dot, map(ev, e.args))
    if h === :Power
        a, b = ev(e.args[1]), ev(e.args[2])
        a isa AbstractMatrix && return a^Int(b)
        return a isa Real && a < 0 && !(b isa Integer) ? complex(a)^b : a^b
    end
    h isa Symbol || throw(ArgumentError("cannot evaluate $(e.head)"))
    args = map(ev, e.args)
    x = isempty(args) ? nothing : args[1]
    h === :Sqrt && return x isa AbstractMatrix ? sqrt(x) : sqrt(complex_if_negative(x))
    h === :Exp && return x isa AbstractMatrix ? exp.(x) : exp(x)
    h === :MatrixExp && return exp(ComplexF64.(x))
    h === :Log && return log(complex_if_negative(x))
    h === :Sin && return sin(x)
    h === :Cos && return cos(x)
    h === :Tan && return tan(x)
    h === :ArcSin && return asin(x)
    h === :ArcCos && return acos(x)
    h === :ArcTan && return length(args) == 2 ? atan(args[2], args[1]) : atan(x)
    h === :Abs && return x isa AbstractArray ? abs.(x) : abs(x)
    h === :Arg && return angle(x)
    h === :Re && return real.(x)
    h === :Im && return imag.(x)
    h === :Sign && return sign(x)
    (h === :Conj || h === :Conjugate) && return conj.(x)
    h in (:Dagger, :ConjugateTranspose, :HC) && return collect(x')
    h === :Transpose && return permutedims(x)
    h === :Inverse && return inv(x)
    h === :Det && return det(x)
    h === :Tr && return tr(x)
    h === :MatrixPower && return x^Int(args[2])
    h === :Commutator && return x * args[2] - args[2] * x
    (h === :Diag || h === :DiagonalMatrix) && return diagm(ComplexF64.(collect(x)))
    h === :IdentityMatrix && return Matrix{ComplexF64}(
        I,
        isempty(args) ? nflavours(ctx.model) : Int(x),
        isempty(args) ? nflavours(ctx.model) : Int(x),
    )
    h === :ZeroMatrix && return zeros(ComplexF64, nflavours(ctx.model), nflavours(ctx.model))
    h === :Max && return maximum(real, args)
    h === :Min && return minimum(real, args)
    h === :Eigenvalues && return eigvals(Hermitian(x))
    h === :Eigenvectors && return eigvecs(Hermitian(x))
    if h === :Rotation
        i, j, θ = Int(args[1]), Int(args[2]), args[3]
        δ = length(args) > 3 ? args[4] : 0.0
        R = Matrix{ComplexF64}(I, nflavours(ctx.model), nflavours(ctx.model))
        R[i, i] = R[j, j] = cos(θ)
        R[i, j] = sin(θ) * cis(-δ)
        R[j, i] = -sin(θ) * cis(δ)
        return R
    end
    if h === :PhaseMatrix
        P = Matrix{ComplexF64}(I, nflavours(ctx.model), nflavours(ctx.model))
        P[Int(args[1]), Int(args[1])] = cis(args[2])
        return P
    end
    if h === :Range
        a, b = real(args[1]), real(args[2])
        step = length(args) > 2 ? real(args[3]) : 1.0
        return collect(a:step:(b + step * 1e-9))
    end
    throw(ArgumentError("unknown function $h"))
end

model_plus(a, b) =
    a isa Number && b isa AbstractArray ? a .+ b :
    b isa Number && a isa AbstractArray ? a .+ b : a + b
# `*` multiplies matrices like `.`: Hamiltonians are written with matrix products
model_times(a, b) =
    a isa AbstractMatrix && b isa AbstractArray ? a * b :
    a isa Number || b isa Number ? a * b : a .* b
model_dot(a, b) = a isa AbstractVector && b isa AbstractVector ? transpose(a) * b : a * b

# ---------------------------------------------------------------------------
#  Parameter values
# ---------------------------------------------------------------------------

function OscParams(m::OscillationModel; kwargs...)
    vals = model_values(m, Dict{Symbol, Float64}())
    for (k, v) in kwargs
        (k === :ρ || k === :rho) && continue
        set_model_value!(m, vals, k, v)
    end
    return model_oscparams(m, vals, get(kwargs, :ρ, get(kwargs, :rho, 1.0)))
end

function OscParams(m::OscillationModel, base::OscParams; kwargs...)
    vals = model_values(m, Dict{Symbol, Float64}())
    for (name, key, _, _) in STANDARD_PARAMETERS[1:6]
        vals[name] = getproperty(base, key)
    end
    for (i, k) in enumerate(base.extra_names)
        p = findfirst(x -> x.key === k, m.parameters)
        p === nothing || (vals[m.parameters[p].name] = base.values[6 + i])
    end
    for (k, v) in kwargs
        (k === :ρ || k === :rho) && continue
        set_model_value!(m, vals, k, v)
    end
    ρ = get(
        kwargs,
        :ρ,
        get(kwargs, :rho, length(base.densities) == 1 ? base.densities[1] : base.densities),
    )
    return model_oscparams(m, vals, ρ)
end

model_values(m::OscillationModel, overrides) =
    Dict{Symbol, Float64}(p.name => get(overrides, p.name, p.value) for p in m.parameters)

function set_model_value!(m::OscillationModel, vals, k::Symbol, v)
    i = findfirst(p -> p.name === canonical_model_name(k) || p.key === canonical(k), m.parameters)
    i === nothing && throw(
        ArgumentError(
            "$k is not a parameter of model $(m.name): $(join((p.name for p in m.parameters), ", "))",
        ),
    )
    vals[m.parameters[i].name] = Float64(v)
end

function model_oscparams(m::OscillationModel, vals, ρ)
    std = [vals[n] for (n, _, _, _) in STANDARD_PARAMETERS[1:6]]
    extras = m.parameters[7:end]
    return OscParams(
        [std; [vals[p.name] for p in extras]],
        [p.key for p in extras],
        Float64.(vcat(ρ)),
    )
end

# Model values (by model name) from an OscParams vector
function values_from_vector(m::OscillationModel, v::AbstractVector)
    vals = Dict{Symbol, Float64}()
    for (i, p) in enumerate(m.parameters)
        vals[p.name] = v[i]
    end
    return vals
end

# ---------------------------------------------------------------------------
#  Hamiltonian and probabilities
# ---------------------------------------------------------------------------

mutable struct ModelEngine
    model::OscillationModel
    cached_values::Vector{Float64}
    static::ModelContext
    K::Matrix{ComplexF64}                         # vacuum mass² matrix
    terms::Vector{Tuple{ModelBlock, Any}}         # (block, value or nothing if energy dependent)
    deco::Union{Nothing, ModelBlock}
    deco_static::Any
end

function ModelEngine(m::OscillationModel)
    eng = ModelEngine(
        m,
        Float64[],
        ModelContext(m, Dict{Symbol, Float64}()),
        zeros(ComplexF64, 0, 0),
        [],
        nothing,
        nothing,
    )
    i = findfirst(b -> b.kind === :Decoherence && !b.off, m.blocks)
    eng.deco = i === nothing ? nothing : m.blocks[i]
    return eng
end

function update!(eng::ModelEngine, v::AbstractVector)
    eng.cached_values == v && return eng
    m = eng.model
    eng.cached_values = collect(Float64, v)
    vals = values_from_vector(m, v)
    eng.static = ModelContext(m, vals)
    if has_block(m, :Hamiltonian)
        eng.K = zeros(ComplexF64, nflavours(m), nflavours(m))
        eng.terms = [(b, nothing) for b in m.blocks if b.kind === :Hamiltonian && !b.off]
    else
        eng.K = vacuum_mass_squared(eng.static)
        eng.terms = [
            (b, try_static(eng.static, b)) for
            b in m.blocks if b.kind === :AdditionalTerms && !b.off
        ]
    end
    eng.deco_static = eng.deco === nothing ? nothing : try_static(eng.static, eng.deco)
    return eng
end

function try_static(ctx::ModelContext, b::ModelBlock)
    try
        return block_value(ctx, b)
    catch err
        err isa DependsOnEnergy && return nothing
        rethrow()
    end
end

# The Hamiltonian in eV for energy E (eV) and matter potential A (eV)
function hamiltonian_matrix(eng::ModelEngine, E, A, anti::Bool; L = nothing, ρ = 0.0)
    m = eng.model
    n = nflavours(m)
    dyn = Dict{Symbol, Any}(:E => E, :A => anti ? -A : A, :CP => anti ? -1 : 1, :rho => ρ)
    L === nothing || (dyn[:L] = L)
    ctx = ModelContext(m, eng.static.values, dyn; parent = eng.static)
    H = eng.K / (2E)
    for (b, val) in eng.terms
        H = H + as_matrix(val === nothing ? block_value(ctx, b) : val, m)
    end
    return anti ? conj(H) : H, ctx
end

function damping_matrix(eng::ModelEngine, ctx::ModelContext)
    eng.deco === nothing && return nothing
    Γ = eng.deco_static === nothing ? block_value(ctx, eng.deco) : eng.deco_static
    n = nflavours(eng.model)
    Γ isa Number && return [i == j ? 0.0 : Float64(real(Γ)) for i = 1:n, j = 1:n]
    return Float64.(real.(Γ))
end

function model_probabilities!(P, eng::ModelEngine, v, E_GeV, anti, lengths, densities)
    update!(eng, v)
    m = eng.model
    n = nflavours(m)
    E = E_GeV * 1e9
    ρ_states = nothing
    S = Matrix{ComplexF64}(I, n, n)
    if eng.deco === nothing
        for k in eachindex(lengths)
            A = densities[k] * GLB_V_FACTOR * GLB_Ne_MANTLE
            L = km_to_ev(lengths[k])
            H, _ = hamiltonian_matrix(eng, E, A, anti; L, ρ = densities[k])
            F = eigen(Hermitian(H))
            S = F.vectors * Diagonal(cis.(-L .* F.values)) * F.vectors' * S
        end
        for α = 1:3, β = 1:3
            P[α, β] = abs2(S[β, α])
        end
        return P
    end
    ρ_states = [
        begin
            r = zeros(ComplexF64, n, n)
            r[α, α] = 1
            r
        end for α = 1:3
    ]
    for k in eachindex(lengths)
        A = densities[k] * GLB_V_FACTOR * GLB_Ne_MANTLE
        L = km_to_ev(lengths[k])
        H, ctx = hamiltonian_matrix(eng, E, A, anti; L, ρ = densities[k])
        F = eigen(Hermitian(H))
        Γ = damping_matrix(eng, ctx)
        D = [
            cis(-(F.values[i] - F.values[j]) * L) * exp(-Γ[max(i, j), min(i, j)] * L) for
            i = 1:n, j = 1:n
        ]
        for α = 1:3
            r = F.vectors' * ρ_states[α] * F.vectors
            ρ_states[α] = F.vectors * (r .* D) * F.vectors'
        end
    end
    for α = 1:3, β = 1:3
        P[α, β] = real(ρ_states[α][β, β])
    end
    return P
end

function use_model!(m::OscillationModel; setup::Setup = current_setup())
    eng = ModelEngine(m)
    P! =
        (P, params, E, anti, lengths, densities) ->
            model_probabilities!(P, eng, params.values, E, anti, lengths, densities)
    probability_engine!(P!, [p.key for p in m.parameters[7:end]]; setup)
    return m
end

"""
    hamiltonian(model, params; E, ρ=0, antineutrino=false) -> Matrix{ComplexF64}

The effective Hamiltonian of a model in eV, for neutrino energy `E` in GeV
in matter of density `ρ` (g/cm³).

```julia
hamiltonian(model, OscParams(model; eta_etau = 0.1); E = 2.5, ρ = 2.848)   # model, parameters; energy in GeV, density in g/cm³
```
"""
function hamiltonian(
    m::OscillationModel,
    p::OscParams;
    E::Real,
    ρ::Real = 0.0,
    antineutrino::Bool = false,
)
    eng = ModelEngine(m)
    update!(eng, p.values)
    H, _ = hamiltonian_matrix(eng, E * 1e9, ρ * GLB_V_FACTOR * GLB_Ne_MANTLE, antineutrino; ρ)
    return H
end

# Eigenvalues and eigenvectors of the Hamiltonian, the matter eigenstates
# numbered like the vacuum mass eigenstates they evolve from. Levels do not
# cross as the density grows, so matter state i has the rank (by eigenvalue)
# of vacuum state i; vacuum states are identified by their overlap with the
# columns of the mixing matrix (by their rank if the model has none).
function matter_eigensystem(
    m::OscillationModel,
    p::OscParams;
    E::Real,
    ρ::Real = 0.0,
    antineutrino::Bool = false,
)
    F = eigen(Hermitian(hamiltonian(m, p; E, ρ, antineutrino)))
    has_block(m, :MixingMatrix) && !has_block(m, :MassMatrix) && !has_block(m, :Hamiltonian) ||
        return F.values, F.vectors
    V0 = eigvecs(Hermitian(hamiltonian(m, p; E, ρ = 0.0, antineutrino)))
    eng = ModelEngine(m)
    update!(eng, p.values)
    U = as_matrix(block_value(eng.static, :MixingMatrix), m)
    antineutrino && (U = conj(U))
    overlap = abs2.(U' * V0)                  # vacuum state i ↔ eigenvalue rank j
    n = size(overlap, 1)
    rank = zeros(Int, n)
    free = trues(n)
    for (i, j) in Tuple.(sort(vec(CartesianIndices(overlap)); by = c -> -overlap[c]))
        rank[i] == 0 && free[j] || continue
        rank[i] = j
        free[j] = false
    end
    return F.values[rank], F.vectors[:, rank]
end

"""
    effective_mixing(model, params; E, ρ=0, antineutrino=false) -> Matrix{ComplexF64}

The effective mixing matrix in matter: the eigenvectors of the
Hamiltonian as columns, numbered like the vacuum mass eigenstates they
evolve from as the density increases (levels do not cross, so θ̃12 → 90°
above the solar resonance for neutrinos).

```julia
effective_mixing(model, OscParams(model); E = 2.5, ρ = 2.848)   # model, parameters; energy in GeV, density in g/cm³
```
"""
effective_mixing(m::OscillationModel, p::OscParams; kwargs...) =
    matter_eigensystem(m, p; kwargs...)[2]

"""
    effective_masses(model, params; E, ρ=0, antineutrino=false) -> Vector{Float64}

The effective masses squared in eV²: 2E times the eigenvalues of the
Hamiltonian, numbered like [`effective_mixing`](@ref).

```julia
effective_masses(model, OscParams(model); E = 2.5, ρ = 2.848, antineutrino = true)   # model, parameters; GeV, g/cm³: m̃² in eV²
```
"""
effective_masses(m::OscillationModel, p::OscParams; E::Real, kwargs...) =
    2E * 1e9 .* matter_eigensystem(m, p; E, kwargs...)[1]

"""
    effective_parameters(model, params; E, ρ=0, antineutrino=false) -> OscParams

The mixing angles, phase and mass splittings in matter (three flavours),
from [`effective_mixing`](@ref) and [`effective_masses`](@ref).

```julia
q = effective_parameters(model, OscParams(model); E = 2.5, ρ = 2.848)   # θ̃12, θ̃13, θ̃23, δ̃ and Δm̃² in matter
```
"""
function effective_parameters(m::OscillationModel, p::OscParams; E::Real, kwargs...)
    nflavours(m) == 3 || throw(ArgumentError("effective parameters are defined for three flavours"))
    λ, V = matter_eigensystem(m, p; E, kwargs...)
    θ12, θ13, θ23, δ = mixing_angles(V)
    m2 = 2E * 1e9 .* λ
    return OscParams(
        θ12 = θ12,
        θ13 = θ13,
        θ23 = θ23,
        δ = δ,
        Δm21 = m2[2] - m2[1],
        Δm31 = m2[3] - m2[1],
    )
end

# Angles and Dirac phase of a 3×3 unitary matrix in the PDG parametrisation
function mixing_angles(U::AbstractMatrix)
    θ13 = asin(min(1.0, abs(U[1, 3])))
    θ12 = atan(abs(U[1, 2]), abs(U[1, 1]))
    θ23 = atan(abs(U[2, 3]), abs(U[3, 3]))
    s12, c12, s13, c13, s23, c23 = sin(θ12), cos(θ12), sin(θ13), cos(θ13), sin(θ23), cos(θ23)
    J = imag(U[1, 1] * U[2, 2] * conj(U[1, 2]) * conj(U[2, 1]))
    nrm = s12 * c12 * s23 * c23 * s13
    nrm > 1e-12 || return θ12, θ13, θ23, 0.0
    cosδ = (abs2(U[2, 1]) - s12^2 * c23^2 - c12^2 * s23^2 * s13^2) / (2nrm)
    sinδ = J / (nrm * c13^2)
    return θ12, θ13, θ23, mod(atan(sinδ, cosδ), 2π)
end

"""
    model_probability(model, params, from, to; E, L, ρ=0, antineutrino=false)

Oscillation probability of a model, also for more than three flavours
(`from`, `to` are flavour names of the model or indices), without loading
it as the probability engine.

```julia
model_probability(model, OscParams(model; eta_etau = 0.1), :mu, :e; E = 2.0, L = 1300, ρ = 2.848)   # model, parameters, from, to; GeV, km, g/cm³
sterile = load_model("examples/models/Sterile")   # four flavours: e, mu, tau, s
model_probability(sterile, OscParams(sterile), :mu, :s; E = 0.5:0.1:5, L = 1300)   # flavours by name, energies as a range
```
"""
function model_probability(
    m::OscillationModel,
    p::OscParams,
    from,
    to;
    E,
    L::Real,
    ρ::Real = 0.0,
    antineutrino::Bool = false,
)
    eng = ModelEngine(m)
    function idx(f)
        f isa Integer && return Int(f)
        i = findfirst(==(Symbol(replace(String(f), "μ" => "mu", "τ" => "tau"))), m.flavours)
        return i === nothing ? flavour(f) : i
    end
    a, b = idx(from), idx(to)
    n = nflavours(m)
    one(Eg) = begin
        update!(eng, p.values)
        Ee = Eg * 1e9
        H, ctx = hamiltonian_matrix(
            eng,
            Ee,
            ρ * GLB_V_FACTOR * GLB_Ne_MANTLE,
            antineutrino;
            L = km_to_ev(L),
            ρ,
        )
        F = eigen(Hermitian(H))
        Lev = km_to_ev(L)
        Γ = damping_matrix(eng, ctx)
        if Γ === nothing
            amp = sum(F.vectors[b, i] * cis(-F.values[i] * Lev) * conj(F.vectors[a, i]) for i = 1:n)
            return abs2(amp)
        end
        r0 = zeros(ComplexF64, n, n)
        r0[a, a] = 1
        D = [
            cis(-(F.values[i] - F.values[j]) * Lev) * exp(-Γ[max(i, j), min(i, j)] * Lev)
            for i = 1:n, j = 1:n
        ]
        r = F.vectors * ((F.vectors' * r0 * F.vectors) .* D) * F.vectors'
        return real(r[b, b])
    end
    return E isa Number ? one(Float64(E)) : [one(Float64(x)) for x in E]
end
