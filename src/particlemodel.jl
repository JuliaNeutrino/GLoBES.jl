# Particle-physics models in Mathematica-style model files.
#
# A model is read from its files (`Model.m`, and optionally `particles.m`,
# `parameters.m`, `SPheno.m`) or built in Julia with the same blocks. GLoBES.jl
# then checks the Lagrangian for gauge and global invariance, computes the
# tree-level scalar spectrum after electroweak symmetry breaking, prints the
# model in LaTeX and, with a neutrino mass mechanism (see `mechanisms.jl`),
# predicts the oscillation parameters.

# ---------------------------------------------------------------------------
#  Model blocks
# ---------------------------------------------------------------------------

"""
    GaugeGroup(boson, group, name, coupling; N=1, expand=false, charge=1)

A factor of the gauge group: `group` is `:U1` or `:SU` (with `N`), as in
the model file's `Gauge[[i]] = {B, U[1], hypercharge, g1, False, 1}`.
"""
struct GaugeGroup
    boson::Symbol
    group::Symbol          # :U1 or :SU
    N::Int
    name::Symbol
    coupling::Symbol
    expand::Bool
    global_charges::Vector{Any}
end

"""
    GlobalSymmetry(kind, name; N=2)

A global symmetry: `kind` is `:Z` (a Z_N parity, charges are the N-th roots of
unity written as ±1 for Z₂, or integers mod N) or `:U1`.
"""
struct GlobalSymmetry
    kind::Symbol
    N::Int
    name::Symbol
end

"""
    ModelField

A gauge multiplet of fermions (left-handed Weyl spinors) or
complex scalars: `generations` copies, the SU(2) `components`, and one charge
per gauge group (the U(1) charge, or the dimension of the SU(N)
representation, negative for the conjugate) followed by the global charges.
"""
struct ModelField
    name::Symbol
    kind::Symbol                     # :fermion or :scalar
    generations::Int
    components::Vector{Symbol}
    conjugated::Bool                 # components given as conj[...] (e.g. d = conj[dR])
    charges::Vector{Any}
end

"""
    LagrangianTerm(name, expr; hc=false)

A named piece of the Lagrangian (`LagFer`, `LagH`, ... in the model file), a parsed
expression, and whether its Hermitian conjugate is added.
"""
struct LagrangianTerm
    name::Symbol
    expr::Any
    hc::Bool
end

"""
    BSMModel

A particle-physics model in the structure of its model file: symmetries, fields, the
Lagrangian and the definitions of the eigenstates after symmetry breaking,
with optional LaTeX names and LesHouches positions of parameters and
particles. Create one with [`load_particle_model`](@ref) or [`BSMModel(name)`](@ref)
and the `add_*!` functions.
"""
mutable struct BSMModel
    name::String
    latex_name::String
    authors::String
    date::String
    global_symmetries::Vector{GlobalSymmetry}
    gauge::Vector{GaugeGroup}
    fields::Vector{ModelField}
    states::Vector{Symbol}
    lagrangian::Vector{LagrangianTerm}
    definitions::Dict{Tuple{Symbol, Symbol}, Any}
    parameter_info::Dict{Symbol, Dict{Symbol, Any}}
    particle_info::Dict{Symbol, Dict{Symbol, Any}}
    spheno::Dict{Symbol, Any}
    assignments::Dict{Symbol, Any}
    source::String
end

"""
    BSMModel(name; latex=name, authors="", date="")

An empty model; add blocks with [`add_gauge!`](@ref), [`add_global!`](@ref),
[`add_fermion!`](@ref), [`add_scalar!`](@ref), [`add_lagrangian!`](@ref),
[`add_vevs!`](@ref) and [`add_parameter!`](@ref).

```julia
m = BSMModel("Scotogenic"; authors = "A. Physicist")   # name, then LaTeX name, authors and date
```
"""
BSMModel(name::AbstractString; latex = name, authors = "", date = "") = BSMModel(
    String(name),
    String(latex),
    String(authors),
    String(date),
    GlobalSymmetry[],
    GaugeGroup[],
    ModelField[],
    [:GaugeES, :EWSB],
    LagrangianTerm[],
    Dict{Tuple{Symbol, Symbol}, Any}(),
    Dict{Symbol, Dict{Symbol, Any}}(),
    Dict{Symbol, Dict{Symbol, Any}}(),
    Dict{Symbol, Any}(),
    Dict{Symbol, Any}(),
    "",
)

fermions(m::BSMModel) = filter(f -> f.kind === :fermion, m.fields)
scalars(m::BSMModel) = filter(f -> f.kind === :scalar, m.fields)

function Base.show(io::IO, m::BSMModel)
    print(
        io,
        "BSMModel(\"",
        m.name,
        "\": ",
        length(fermions(m)),
        " fermion and ",
        length(scalars(m)),
        " scalar multiplets, ",
        length(m.lagrangian),
        " Lagrangian terms)",
    )
end

function Base.show(io::IO, ::MIME"text/plain", m::BSMModel)
    println(io, "BSMModel \"", m.name, "\"", isempty(m.authors) ? "" : " by $(m.authors)")
    group =
        join((g.group === :U1 ? "U(1)_$(g.name)" : "SU($(g.N))_$(g.name)" for g in m.gauge), " × ")
    println(io, "  gauge group: ", group)
    isempty(m.global_symmetries) || println(
        io,
        "  global:      ",
        join((s.kind === :Z ? "Z$(s.N)" : "U(1)" for s in m.global_symmetries), " × "),
    )
    for f in m.fields
        @printf(
            io,
            "  %-8s %-7s ×%d  (%s)\n",
            f.name,
            f.kind,
            f.generations,
            join(map(charge_string, f.charges), ", ")
        )
    end
    print(io, "  Lagrangian:  ", join((t.name for t in m.lagrangian), ", "))
end

charge_string(q::Rational) =
    denominator(q) == 1 ? string(numerator(q)) : "$(numerator(q))/$(denominator(q))"
charge_string(q) = string(q)

# ---------------------------------------------------------------------------
#  Building a model in Julia
# ---------------------------------------------------------------------------

"""
    add_gauge!(model, boson, group, name, coupling; expand=false, charge=1)

Add a gauge group: `group` is `:U1`, or `(:SU, N)`. As in the model file,
`Gauge[[i]] = {boson, group, name, coupling, expand, charge}`.

```julia
add_gauge!(m, :B, :U1, :hypercharge, :g1)   # model, gauge boson, group, name of the group, coupling
add_gauge!(m, :WB, (:SU, 2), :left, :g2; expand = true)   # SU(2), with its indices expanded in expressions
add_gauge!(m, :G, (:SU, 3), :color, :g3)   # SU(3)
```
"""
function add_gauge!(m::BSMModel, boson, group, name, coupling; expand = false, charge = 1)
    kind, N = group === :U1 ? (:U1, 1) : (Symbol(group[1]), Int(group[2]))
    push!(
        m.gauge,
        GaugeGroup(Symbol(boson), kind, N, Symbol(name), Symbol(coupling), expand, Any[charge...]),
    )
    return m
end

"""
    add_global!(model, name; Z=2)
    add_global!(model, name; U1=true)

Add a global Z_N or U(1) symmetry.

```julia
add_global!(m, :Z2; Z = 2)   # model, name, a Z_N symmetry with N = 2
add_global!(m, :Lepton; U1 = true)   # a global U(1)
```
"""
function add_global!(m::BSMModel, name; Z = 2, U1 = false)
    push!(
        m.global_symmetries,
        U1 ? GlobalSymmetry(:U1, 1, Symbol(name)) : GlobalSymmetry(:Z, Z, Symbol(name)),
    )
    return m
end

"""
    add_fermion!(model, name, generations, components, charges...; conjugated=false)
    add_scalar!(model, name, generations, components, charges...)

Add a multiplet of left-handed Weyl fermions or complex scalars, with one
charge per gauge group in the order they were added, then one per global
symmetry. `conjugated=true` marks right-handed fields entered as their
conjugates, like `conj[dR]` in the model file.

```julia
add_fermion!(m, :l, 3, (:vL, :eL), -1 // 2, 2, 1, 1)   # model, name, generations, SU(2) components, hypercharge, SU(2) and SU(3) dimensions, Z2 charge
add_fermion!(m, :e, 3, :eR, 1, 1, 1, 1; conjugated = true)   # a right-handed field, entered as conj[eR]
add_fermion!(m, :n, 3, :nR, 0, 1, 1, -1; conjugated = true)   # the Z2-odd singlet fermions
add_scalar!(m, :H, 1, (:Hp, :H0), 1 // 2, 2, 1, 1)   # the Higgs doublet, entries as for fermions
add_scalar!(m, :Et, 1, (:etp, :et0), 1 // 2, 2, 1, -1)   # the Z2-odd inert doublet
```
"""
add_fermion!(m::BSMModel, name, generations, components, charges...; conjugated = false) =
    add_field!(m, :fermion, name, generations, components, charges, conjugated)
add_scalar!(m::BSMModel, name, generations, components, charges...; conjugated = false) =
    add_field!(m, :scalar, name, generations, components, charges, conjugated)
@doc (@doc add_fermion!) add_scalar!

function add_field!(m, kind, name, generations, components, charges, conjugated)
    comps = components isa Symbol ? [components] : collect(Symbol, components)
    push!(
        m.fields,
        ModelField(Symbol(name), kind, Int(generations), comps, conjugated, collect(Any, charges)),
    )
    return m
end

"""
    add_lagrangian!(model, name, expr; hc=false)

Add a piece of the Lagrangian, written as in the model file:

```julia
add_lagrangian!(m, :LagFer, "Ye conj[H].e.l + Yn Et.n.l"; hc = true)   # model, name, terms as in the model file, add the Hermitian conjugate
add_lagrangian!(m, :LagNV, "Mn/2 n.n"; hc = true)   # Majorana masses of the singlet fermions
add_lagrangian!(m, :LagH, "-(-mH2 conj[H].H + 1/2 lambda1 conj[H].H.conj[H].H)")   # self-conjugate: no Hermitian conjugate
```
"""
function add_lagrangian!(m::BSMModel, name, expr; hc = false)
    e = expr isa AbstractString ? only(parse_mathematica(expr)) : expr
    push!(m.lagrangian, LagrangianTerm(Symbol(name), e, hc))
    m.assignments[Symbol(name)] = e
    return m
end

"""
    add_vevs!(model, (component, vev, (pseudoscalar, coefficient), (scalar, coefficient)), ...)

Decompose neutral scalars after electroweak symmetry breaking, like
`DEFINITION[EWSB][VEVs]` in the model file. Coefficients are numbers or expressions:

```julia
add_vevs!(m, (:H0, :v, (:Ah, "I/Sqrt[2]"), (:hh, "1/Sqrt[2]")), (:et0, 0, (:etI, "I/Sqrt[2]"), (:etR, "1/Sqrt[2]")))   # model, then (neutral component, vev, (CP-odd state, coefficient), (CP-even state, coefficient)) for each
```
"""
function add_vevs!(m::BSMModel, entries...)
    to_m(x) = x isa AbstractString ? only(parse_mathematica(x)) : x
    list = MExpr(
        :List,
        Any[
            MExpr(
                :List,
                Any[
                    Symbol(e[1]),
                    e[2] == 0 ? MExpr(:List, Any[0, 0]) :
                    MExpr(:List, Any[Symbol(e[2]), MExpr(:Power, MExpr(:Sqrt, 2), -1)]),
                    MExpr(:List, Any[Symbol(e[3][1]), to_m(e[3][2])]),
                    MExpr(:List, Any[Symbol(e[4][1]), to_m(e[4][2])]),
                ],
            ) for e in entries
        ],
    )
    m.definitions[(:EWSB, :VEVs)] = list
    return m
end

"""
    add_parameter!(model, name; latex, description, real, leshouches)

LaTeX name, description and LesHouches position (`(:HDM, 1)` for an entry,
`:YN` for a matrix block) of a parameter, like a line of
`parameters.m`.

```julia
add_parameter!(m, :Yn; description = "Yukawa of the singlet fermions", leshouches = :YN)   # model, name; description, LesHouches block
add_parameter!(m, :lambda5; leshouches = (:HDM, 5))   # an entry of a LesHouches block
```
"""
function add_parameter!(m::BSMModel, name; kwargs...)
    info = get!(m.parameter_info, Symbol(name), Dict{Symbol, Any}())
    for (k, v) in kwargs
        info[k == :latex ? :LaTeX :
             k == :description ? :Description :
             k == :real ? :Real : k == :leshouches ? :LesHouches : k] = v
    end
    return m
end

# ---------------------------------------------------------------------------
#  Reading model files
# ---------------------------------------------------------------------------

"""
    load_particle_model(dir_or_file; particles, parameters, spheno) -> BSMModel

Read a particle model: the main model file and, if they exist next to it,
`particles.m`, `parameters.m` and `SPheno.m`. `dir_or_file` is the model
directory (which must contain `<name>.m`, or exactly one other model file)
or the main file itself.

```julia
scoto = load_particle_model("Scotogenic/")   # model directory with Scotogenic.m and particles.m, parameters.m, SPheno.m
scoto = load_particle_model("Scotogenic/Scotogenic.m"; spheno = "Scotogenic/SPheno.m")   # the model file, the others by name
```
"""
function load_particle_model(
    path::AbstractString;
    particles = nothing,
    parameters = nothing,
    spheno = nothing,
)
    path = expanduser(path)
    main = isdir(path) ? find_model_file(path) : path
    dir = dirname(main)
    m = BSMModel(splitext(basename(main))[1])
    m.source = main
    read_particle_file!(m, main)
    for (file, default) in
        ((particles, "particles.m"), (parameters, "parameters.m"), (spheno, "SPheno.m"))
        f = something(file, joinpath(dir, default))
        isfile(f) && read_particle_file!(m, f)
    end
    return m
end

function find_model_file(dir)
    candidates = filter(readdir(dir)) do f
        endswith(f, ".m") && !(lowercase(f) in ("particles.m", "parameters.m", "spheno.m"))
    end
    name = basename(rstrip(dir, '/')) * ".m"
    name in candidates && return joinpath(dir, name)
    length(candidates) == 1 && return joinpath(dir, only(candidates))
    throw(ArgumentError("cannot tell which model file in $dir to read: $(join(candidates, ", "))"))
end

"""
    read_particle_file!(model, file)

Read the statements of one model file into `model`. Unknown statements are
kept in `model.assignments`.
"""
function read_particle_file!(m::BSMModel, file::AbstractString)
    read_particle_string!(m, read(file, String))
    return m
end

function read_particle_string!(m::BSMModel, src::AbstractString)
    lag_input = nothing
    for st in parse_mathematica(src)
        is_call(st, :Set) || is_call(st, :SetDelayed) || continue
        lhs, rhs = st.args
        if lhs isa Symbol
            name = String(lhs)
            if startswith(name, "Model`")
                key = name[7:end]
                key == "Name" && (m.name = rhs)
                key == "NameLaTeX" && (m.latex_name = rhs)
                key == "Authors" && (m.authors = rhs)
                key == "Date" && (m.date = rhs)
            elseif lhs === :NameOfStates
                m.states = Symbol[s for s in rhs.args]
            elseif lhs === :ParameterDefinitions
                for entry in rhs.args
                    merge!(
                        get!(m.parameter_info, entry.args[1], Dict{Symbol, Any}()),
                        options(entry.args[2]),
                    )
                end
            elseif lhs === :WeylFermionAndIndermediate || lhs === :WeylFermionAndIntermediate
                for entry in rhs.args
                    merge!(
                        get!(m.particle_info, entry.args[1], Dict{Symbol, Any}()),
                        options(entry.args[2]),
                    )
                end
            elseif lhs in (
                :OnlyLowEnergySPheno,
                :MINPAR,
                :ParametersToSolveTadpoles,
                :BoundaryLowScaleInput,
                :ListDecayParticles,
                :ListDecayParticles3B,
                :EXTPAR,
                :BoundaryEWSBScale,
            )
                m.spheno[lhs] = rhs
            else
                m.assignments[lhs] = rhs
            end
        elseif is_call(lhs, :Part)
            block, i = lhs.args[1], lhs.args[2]
            if block === :Gauge
                set_index!(m.gauge, i, gauge_entry(rhs))
            elseif block === :Global
                set_index!(m.global_symmetries, i, global_entry(rhs))
            elseif block === :FermionFields
                set_index!(m.fields, field_slot(m, :fermion, i), field_entry(:fermion, rhs))
            elseif block === :ScalarFields
                set_index!(m.fields, field_slot(m, :scalar, i), field_entry(:scalar, rhs))
            else
                m.assignments[Symbol(block, "_", i)] = rhs
            end
        elseif lhs isa MExpr && lhs.head isa MExpr && lhs.head.head === :DEFINITION
            state, key = lhs.head.args[1], lhs.args[1]
            m.definitions[(state, key)] = rhs
            key === :LagrangianInput && (lag_input = rhs)
        elseif is_call(lhs, :ParticleDefinitions)
            for entry in rhs.args
                info = get!(m.particle_info, entry.args[1], Dict{Symbol, Any}())
                for (k, v) in options(entry.args[2])
                    haskey(info, k) || (info[k] = v)        # the first definition wins
                end
                get!(info, :states, Symbol[])
                push!(info[:states], lhs.args[1])
            end
        elseif is_call(lhs, :DEFINITION)
            m.spheno[Symbol(lhs.args[1])] = rhs
        end
    end
    if lag_input !== nothing
        for entry in lag_input.args
            name = entry.args[1]
            opts = options(entry.args[2])
            haskey(m.assignments, name) || throw(
                ArgumentError(
                    "Lagrangian piece $name is listed in LagrangianInput but never defined",
                ),
            )
            filter!(t -> t.name !== name, m.lagrangian)
            push!(
                m.lagrangian,
                LagrangianTerm(name, m.assignments[name], get(opts, :AddHC, false) === :True),
            )
        end
    end
    return m
end

# Rules {A->x, B->y} as a Dict; option names like PDG.IX are joined with dots
function options(list)
    d = Dict{Symbol, Any}()
    for r in (list isa MExpr && list.head === :List ? list.args : Any[list])
        is_call(r, :Rule) || continue
        key = r.args[1] isa Symbol ? r.args[1] : Symbol(join(string.(r.args[1].args), "."))
        d[key] = r.args[2]
    end
    return d
end

function set_index!(v::Vector, i::Integer, x)
    while length(v) < i
        push!(v, x)
    end
    v[i] = x
end

# Fermions and scalars share m.fields; the model file numbers them separately
function field_slot(m, kind, i)
    k = 0
    for (j, f) in enumerate(m.fields)
        f.kind === kind && (k += 1) == i && return j
    end
    return length(m.fields) + 1
end

function gauge_entry(rhs)
    boson, group, name, coupling, expand, charges... = rhs.args
    kind, N = group.head === :U ? (:U1, 1) : (:SU, Int(group.args[1]))
    GaugeGroup(boson, kind, N, name, coupling, expand === :True, collect(Any, charges))
end

function global_entry(rhs)
    kind, name = rhs.args
    kind.head === :Z ? GlobalSymmetry(:Z, Int(kind.args[1]), name) : GlobalSymmetry(:U1, 1, name)
end

function field_entry(kind, rhs)
    name, gen, comps, charges... = rhs.args
    conjugated = is_call(comps, :conj)
    list =
        conjugated ? Symbol[comps.args[1]] :
        comps isa Symbol ? Symbol[comps] : Symbol[c for c in comps.args]
    ModelField(name, kind, Int(gen), list, conjugated, collect(Any, charges))
end

# ---------------------------------------------------------------------------
#  Monomials of the Lagrangian
# ---------------------------------------------------------------------------

"""
    Monomial

One term `coefficient × parameters × fields` of an expanded Lagrangian piece.
`fields` lists `(name, conjugated)` in the order of the `Dot` product.
"""
struct Monomial
    coefficient::Any
    parameters::Vector{Any}
    fields::Vector{Tuple{Symbol, Bool}}
end

# The multiplet a field symbol belongs to (multiplet or component name)
function field_lookup(m::BSMModel, s::Symbol)
    for f in m.fields
        (f.name === s || s in f.components) && return f
    end
    return nothing
end

"""
    expand_lagrangian(model, expr) -> Vector{Monomial}

Expand a Lagrangian expression into monomials.
"""
function expand_lagrangian(m::BSMModel, expr)
    mono(c, p, f) = [Monomial(c, p, f)]
    function ex(e)
        e isa MNumber && return mono(e, Any[], Tuple{Symbol, Bool}[])
        if e isa Symbol
            field_lookup(m, e) === nothing || return mono(1, Any[], [(e, false)])
            return mono(1, Any[e], Tuple{Symbol, Bool}[])
        end
        if is_call(e, :conj) && e.args[1] isa Symbol && field_lookup(m, e.args[1]) !== nothing
            return mono(1, Any[], [(e.args[1], true)])
        end
        if is_call(e, :Plus)
            return reduce(vcat, (ex(a) for a in e.args))
        end
        if is_call(e, :Times) || is_call(e, :Dot)
            out = mono(1, Any[], Tuple{Symbol, Bool}[])
            for a in e.args
                out = [
                    Monomial(
                        x.coefficient * y.coefficient,
                        [x.parameters; y.parameters],
                        [x.fields; y.fields],
                    ) for x in out for y in ex(a)
                ]
            end
            return out
        end
        if is_call(e, :Power) && e.args[2] isa Integer && e.args[2] > 0
            return ex(MExpr(:Times, fill(e.args[1], e.args[2])))
        end
        return mono(1, Any[e], Tuple{Symbol, Bool}[])      # Sqrt[2], functions, ...
    end
    return ex(expr)
end

# ---------------------------------------------------------------------------
#  Invariance checks
# ---------------------------------------------------------------------------

"""
    ModelCheck

Result of [`check_model`](@ref): `problems` lists the terms that break a
symmetry, `warnings` other remarks (non-renormalisable terms, parameters
without a definition, ...).
"""
struct ModelCheck
    model::String
    terms::Int
    problems::Vector{String}
    warnings::Vector{String}
end

passed(c::ModelCheck) = isempty(c.problems)

function Base.show(io::IO, ::MIME"text/plain", c::ModelCheck)
    println(
        io,
        "Model check of \"",
        c.model,
        "\": ",
        c.terms,
        " terms, ",
        passed(c) ? "all invariant" : "$(length(c.problems)) problems",
    )
    for p in c.problems
        println(io, "  ✗ ", p)
    end
    for w in c.warnings
        println(io, "  ! ", w)
    end
end

"""
    check_model(model) -> ModelCheck

Check every term of the Lagrangian: the sum of U(1) charges vanishes, the
SU(N) representations can form a singlet (their N-alities add up to zero
mod N), the global Z_N or U(1) charges are conserved, each term has 0 or 2
fermions and mass dimension at most 4.

```julia
c = check_model(scoto)   # c.problems: terms that break a gauge or global symmetry; c.warnings
```
"""
function check_model(m::BSMModel)
    problems, warnings = String[], String[]
    n_terms = 0
    known = Set(keys(m.parameter_info))
    unknown = Set{Symbol}()
    for t in m.lagrangian
        for (k, mono) in enumerate(expand_lagrangian(m, t.expr))
            n_terms += 1
            isempty(mono.fields) && continue
            where_ = "$(t.name), term $k ($(monomial_string(mono)))"
            for (g, group) in enumerate(m.gauge)
                if group.group === :U1
                    total = sum(field_charge(m, f, g, c) for (f, c) in mono.fields)
                    total == 0 ||
                        push!(problems, "$where_: $(group.name) charge $(charge_string(total)) ≠ 0")
                else
                    nalities =
                        [n_ality(field_charge(m, f, g, c), group.N) for (f, c) in mono.fields]
                    nality = any(isnothing, nalities) ? nothing : sum(nalities)
                    if nality === nothing
                        push!(
                            warnings,
                            "$where_: cannot check $(group.name) for these representations",
                        )
                    elseif mod(nality, group.N) != 0
                        push!(
                            problems,
                            "$where_: no $(group.name) singlet (SU($(group.N)) N-ality $nality)",
                        )
                    end
                end
            end
            for (s, sym) in enumerate(m.global_symmetries)
                charges = [field_charge(m, f, length(m.gauge) + s, c) for (f, c) in mono.fields]
                if sym.kind === :Z && all(q -> q in (1, -1), charges) && sym.N == 2
                    prod(charges) == 1 || push!(problems, "$where_: $(sym.name) parity is odd")
                else
                    total = sum(charges)
                    ok = sym.kind === :Z ? mod(total, sym.N) == 0 : total == 0
                    ok || push!(
                        problems,
                        "$where_: $(sym.name) charge $(charge_string(total)) not conserved",
                    )
                end
            end
            nf = count(((f, c),) -> field_lookup(m, f).kind === :fermion, mono.fields)
            ns = length(mono.fields) - nf
            nf in (0, 2) || push!(problems, "$where_: $nf fermions (a Lorentz scalar needs 0 or 2)")
            dim = ns + 3nf / 2
            dim <= 4 || push!(warnings, "$where_: mass dimension $(dim) > 4 (not renormalisable)")
            for p in mono.parameters
                p isa Symbol && !(p in known) && !(p in (:I, :v)) && push!(unknown, p)
            end
        end
    end
    isempty(m.parameter_info) ||
        isempty(unknown) ||
        push!(
            warnings,
            "parameters not in parameters.m (no LaTeX name or LesHouches entry): " *
            join(sort!(collect(unknown)), ", "),
        )
    return ModelCheck(m.name, n_terms, problems, warnings)
end

# The charge of a field under symmetry k, conjugated if c. The model file gives
# the charges of the multiplet; a component of a multiplet entered as conj[...]
# (d = conj[dR]) has the conjugate charges.
function field_charge(m, name, k, c)
    f = field_lookup(m, name)
    q = f.charges[k]
    f.name !== name && f.conjugated && (q = conjugate_charge(m, k, q))
    return c ? conjugate_charge(m, k, q) : q
end

# U(1) charges and SU(N) dimensions change sign; Z₂ parities ±1 stay
function conjugate_charge(m, k, q)
    k <= length(m.gauge) && return -q
    sym = m.global_symmetries[k - length(m.gauge)]
    return sym.kind === :Z && sym.N == 2 && q in (1, -1) ? q : -q
end

# N-ality of an SU(N) representation given by its signed dimension
function n_ality(d, N)
    d = Int(d)
    abs(d) == 1 && return 0
    abs(d) == N && return d > 0 ? 1 : N - 1
    abs(d) == N^2 - 1 && return 0
    N == 3 && abs(d) == 6 && return d > 0 ? 2 : 1
    N == 2 && return isodd(abs(d)) ? 0 : 1           # spin j: 2j+1
    return nothing
end

function monomial_string(mono::Monomial)
    parts = String[]
    mono.coefficient == 1 || push!(parts, charge_string(mono.coefficient))
    append!(parts, string.(mono.parameters))
    append!(parts, (c ? "conj[$f]" : string(f) for (f, c) in mono.fields))
    join(parts, " ")
end
