# The control file of a model directory (Run.m): what to calculate, plot,
# fit and write to LaTeX, written like SARAH's SPheno.m. Each entry is
# {name, {options}} and the entries run in the order of the file.
#
#   Calculate[[1]] = {EffectiveMasses, {E -> 2.5, rho -> 2.848}};
#   Plot[[1]] = {Pmue, {Type -> Probability, From -> mu, To -> e, E -> Range[0.5, 20, 0.05], L -> 1300,
#                       Curves -> {"SM" -> {}, "NSI" -> {eta_etau -> 0.1}}, File -> "pmue.pdf"}};
#   Fit[[1]] = {EtaScan, {Experiment -> "DUNE.glb", Scan -> {eta_etau -> Range[-0.04, 0.04, 0.004]},
#                         Free -> {theta23, deltaCP, dm31}, Priors -> NuFit, Output -> "eta.dat"}};
#   LaTeX[[1]] = {"results.tex", {Title -> "...", Parameters -> On, MixingMatrix -> Off, Text -> {"..."}}};
#
# Plain numbers are in GeV (E), km (L), g/cm³ (rho) and radians (angles);
# with units (2.5 GeV, 1300 km, 198 Degree) they are converted.

"""
    ModelRun

What [`run_model`](@ref) did: the `model`, the results of `Calculate`
entries (`calculations`), of `Fit` entries (`fits`), the figure files
(`figures`) and the LaTeX file (`tex`, or `nothing`).
"""
struct ModelRun
    model::OscillationModel
    calculations::Vector{Any}       # (title, value, LaTeX)
    fits::Vector{Any}               # (name, NamedTuple of results)
    figures::Vector{Tuple{String, String}}
    tex::Union{Nothing, String}
end

function Base.show(io::IO, ::MIME"text/plain", r::ModelRun)
    println(
        io,
        "Run of model \"",
        r.model.name,
        "\": ",
        length(r.calculations),
        " calculations, ",
        length(r.fits),
        " fits, ",
        length(r.figures),
        " figures",
    )
    r.tex === nothing || print(io, "  LaTeX: ", r.tex)
end

"""
    run_model(dir; run="Run.m", output=dir, plots=true) -> ModelRun

Load the model in `dir` (see [`load_model`](@ref)) and execute its control
file: its `Calculate[[i]]`, `Plot[[i]]`, `Fit[[i]]`, `LaTeX[[i]]` entries
and `Print[...]` statements in the order of the file.
Output files go to `output`. Plot entries need `using Plots` (with
`plots=false` they are skipped).
"""
function run_model(dir::AbstractString; run = nothing, output = nothing, plots::Bool = true)
    m = load_model(dir)
    file = something(run, get(m.files, :run, nothing), joinpath(dir, "Run.m"))
    isfile(file) || throw(ArgumentError("no control file $file"))
    out = something(output, dirname(file))
    mkpath(out)
    statements = parse_mathematica(read(file, String))
    state = RunState(
        m,
        out,
        Setup(),
        Any[],
        Any[],
        Tuple{String, String}[],
        nothing,
        plots,
        Dict{Symbol, Any}(),
    )
    for st in statements
        execute!(state, st)
    end
    return ModelRun(m, state.calculations, state.fits, state.figures, state.tex)
end

mutable struct RunState
    model::OscillationModel
    output::String
    setup::Setup
    calculations::Vector{Any}
    fits::Vector{Any}
    figures::Vector{Tuple{String, String}}
    tex::Union{Nothing, String}
    plots::Bool
    settings::Dict{Symbol, Any}      # top-level defaults: E, L, rho, Hierarchy, ...
end

const RUN_ENTRIES = (:Calculate, :Plot, :Fit, :LaTeX)

function execute!(s::RunState, st)
    if is_call(st, :Print)
        println(join(string.(run_value(s, a, Dict()) for a in st.args)))
        return
    end
    is_call(st, :Set) ||
        throw(ArgumentError("unknown statement in control file: $(sprint(show_mexpr, st))"))
    lhs, rhs = st.args
    if lhs isa Symbol
        lhs === :Hierarchy ? (s.model.hierarchy = rhs) : (s.settings[lhs] = rhs)
        return
    end
    kind = is_call(lhs, :Part) ? lhs.args[1] : nothing
    kind in RUN_ENTRIES || throw(
        ArgumentError(
            "unknown entry $(sprint(show_mexpr, lhs)) in control file; use $(join(("$k[[i]]" for k in RUN_ENTRIES), ", "))",
        ),
    )
    is_call(rhs, :List) && length(rhs.args) == 2 && is_call(rhs.args[2], :List) ||
        throw(ArgumentError("write $kind[[i]] = {name, {option -> value, ...}}"))
    arg, list = rhs.args
    opts = Dict{Symbol, Any}()
    for r in list.args
        is_call(r, :Rule) ||
            throw(ArgumentError("options of $kind[[i]] are rules, not $(sprint(show_mexpr, r))"))
        opts[r.args[1]] = r.args[2]
    end
    if kind === :Calculate
        calculate!(s, arg, opts)
    elseif kind === :Plot
        s.plots && run_plot!(s, arg, opts)
    elseif kind === :Fit
        run_fit!(s, arg, opts)
    else
        run_latex!(s, arg, opts)
    end
end

# Option value with a fallback to the top-level settings
option(s::RunState, opts, k, default = nothing) = get(opts, k, get(s.settings, k, default))

# A value in the control file: units convert to GeV (E), km (L) or radians
function run_value(s::RunState, e, opts; unit = nothing)
    e === nothing && return nothing
    e isa String && return e
    ctx = ModelContext(s.model, model_values(s.model, Dict{Symbol, Float64}()))
    v = model_eval(e, ctx, 0)
    v isa AbstractArray && (v = real.(v) == v ? Float64.(real.(v)) : v)
    unit === nothing && return v
    return mentions_units(e) ? v ./ MODEL_CONSTANTS[unit] : v
end

function mentions_units(e)
    e isa Symbol && return e in (:eV, :keV, :MeV, :GeV, :TeV, :km, :m, :cm)
    e isa MExpr && return any(mentions_units, e.args)
    return false
end

# Parameters: the model defaults with the overrides of a rule list
function run_params(s::RunState, overrides)
    m = s.model
    kw = Dict{Symbol, Float64}()
    if overrides !== nothing
        for r in (is_call(overrides, :List) ? overrides.args : Any[overrides])
            is_call(r, :Rule) || continue
            kw[r.args[1]] = Float64(real(run_value(s, r.args[2], Dict())))
        end
    end
    return OscParams(m; kw...)
end

flavour_arg(x) = x isa String ? Symbol(x) : x === :mu ? :μ : x === :tau ? :τ : x

# ---------------------------------------------------------------------------
#  Calculate
# ---------------------------------------------------------------------------

function calculate!(s::RunState, what, opts)
    m = s.model
    p = run_params(s, get(opts, :Parameters, nothing))
    E = run_value(s, option(s, opts, :E, 1.0), opts; unit = :GeV)
    ρ = Float64(real(run_value(s, option(s, opts, :rho, 0.0), opts)))
    anti = option(s, opts, :Antineutrino, :False) === :True
    names = latex_names(m)
    where_ = "E = $(E) GeV, ρ = $(ρ) g/cm³" * (anti ? ", antineutrinos" : "")
    tex_where =
        "E = $(latex_number(first(E)))\\,\\mathrm{GeV}, \\rho = $(latex_number(ρ))\\,\\mathrm{g/cm^3}" *
        (anti ? ", \\bar\\nu" : "")
    title, value, tex_ = if what === :EffectiveMixingMatrix
        U = effective_mixing(m, p; E, ρ, antineutrino = anti)
        "Effective mixing matrix ($where_)", U, "\\tilde U = " * latex(U)
    elseif what === :EffectiveMasses
        m2 = effective_masses(m, p; E, ρ, antineutrino = anti)
        "Effective masses squared in eV² ($where_)",
        m2,
        "\\tilde m^2_i = " * latex(m2) * "\\,\\mathrm{eV}^2"
    elseif what === :EffectiveParameters
        q = effective_parameters(m, p; E, ρ, antineutrino = anti)
        "Effective parameters ($where_)",
        q,
        "\\sin^2\\tilde\\theta_{12} = $(latex_number(sin(q.θ12)^2)),\\ \\sin^2\\tilde\\theta_{13} = $(latex_number(sin(q.θ13)^2)),\\ \\sin^2\\tilde\\theta_{23} = $(latex_number(sin(q.θ23)^2)),\\ \\tilde\\delta = $(latex_number(rad2deg(q.δ)))^\\circ"
    elseif what === :Hamiltonian
        H = hamiltonian(m, p; E, ρ, antineutrino = anti)
        "Hamiltonian in eV ($where_)", H, "H = " * latex(H) * "\\,\\mathrm{eV}"
    elseif what === :MixingMatrix
        U = block_value(
            ModelContext(
                m,
                model_values(m, Dict(k.name => p.values[i] for (i, k) in enumerate(m.parameters))),
            ),
            :MixingMatrix,
        )
        "Vacuum mixing matrix", U, "U = " * latex(U)
    elseif what === :Probability
        from, to = flavour_arg(option(s, opts, :From, :mu)), flavour_arg(option(s, opts, :To, :e))
        L = Float64(real(run_value(s, option(s, opts, :L, 1300), opts; unit = :km)))
        P = model_probability(m, p, from, to; E, L, ρ, antineutrino = anti)
        "P($from → $to) at L = $L km ($where_)",
        P,
        "P(\\nu_{$from} \\to \\nu_{$to}) = " *
        (P isa Number ? latex_number(P) : latex(P)) *
        " \\quad (L = $(latex_number(L))\\,\\mathrm{km})"
    else
        throw(
            ArgumentError(
                "unknown calculation $what; use EffectiveMixingMatrix, EffectiveMasses, EffectiveParameters, Hamiltonian, MixingMatrix or Probability",
            ),
        )
    end
    println(title, ":")
    show(stdout, MIME"text/plain"(), value)
    println("\n")
    push!(s.calculations, (title, value, "\\[ " * tex_ * " \\quad (" * tex_where * ") \\]"))
end

# ---------------------------------------------------------------------------
#  Fits
# ---------------------------------------------------------------------------

function key_of(m::OscillationModel, name::Symbol)
    name = canonical_model_name(name)
    i = findfirst(p -> p.name === name || p.key === canonical(name), m.parameters)
    i === nothing && throw(ArgumentError("$name is not a parameter of model $(m.name)"))
    return m.parameters[i].key
end

function run_fit!(s::RunState, name, opts)
    m = s.model
    setup = Setup()
    exp_file = get(opts, :Experiment, nothing)
    exp_file === nothing && throw(ArgumentError("Fit[$name] needs an Experiment"))
    files = exp_file isa String ? [exp_file] : [x for x in exp_file.args]
    with_setup(setup) do
        for f in files
            path = isabspath(f) ? f : joinpath(m.dir, f)
            cd(() -> Experiment(basename(path)), dirname(path))
        end
    end
    use_model!(m; setup)
    truth = run_params(s, get(opts, :Truth, nothing))
    test = haskey(opts, :Test) ? run_params(s, opts[:Test]) : truth
    free = haskey(opts, :Free) ? Tuple(key_of(m, f) for f in opts[:Free].args) : ()
    prior, priors = nothing, NamedTuple()
    if get(opts, :Priors, nothing) === :NuFit
        t = nufit_chi2(m.nufit.version; sk = m.nufit.sk)
        prior =
            nufit_prior(t; use = Tuple(f for f in free if f in (:θ12, :θ13, :θ23, :Δm21, :Δm31)))
    elseif is_call(get(opts, :Priors, nothing), :List)
        priors = NamedTuple(
            key_of(m, r.args[1]) => Float64(real(run_value(s, r.args[2], opts))) for
            r in opts[:Priors].args
        )
    end
    type = get(opts, :Type, :Scan)
    result = if type === :Scan
        scan = only(opts[:Scan].args)
        key = key_of(m, scan.args[1])
        values = Float64.(real.(run_value(s, scan.args[2], opts)))
        χ² = chi2_profile(truth, key => values; test, free, priors, prior, setup)
        (parameter = scan.args[1], x = values, chi2 = χ² .- minimum(χ²))
    elseif type in (:MassOrdering, :CPViolation)
        δs = Float64.(
            real.(run_value(s, get(opts, :DeltaTrue, MExpr(:Range, Any[-π, π, π / 10])), opts)),
        )
        f = type === :MassOrdering ? mass_ordering_sensitivity : cp_violation_sensitivity
        kw =
            type === :MassOrdering ? (free = isempty(free) ? (:θ23, :δ, :Δm31) : free,) :
            (free = isempty(free) ? (:θ23, :Δm31) : free,)
        σ = f(truth, δs; kw..., priors, prior, setup)
        (parameter = :deltaCP, x = δs, sigma = σ)
    else
        throw(ArgumentError("unknown fit type $type; use Scan, MassOrdering or CPViolation"))
    end
    push!(s.fits, (String(name), result))
    out = get(opts, :Output, nothing)
    if out !== nothing
        open(joinpath(s.output, out), "w") do io
            y = haskey(result, :chi2) ? result.chi2 : result.sigma
            println(
                io,
                "# ",
                name,
                ": ",
                result.parameter,
                "  ",
                haskey(result, :chi2) ? "Δχ²" : "σ",
            )
            for (x, v) in zip(result.x, y)
                @printf(io, "%.8g %.8g\n", x, v)
            end
        end
    end
    println("Fit[$name] done", out === nothing ? "" : ", written to $out")
end

# ---------------------------------------------------------------------------
#  Plots
# ---------------------------------------------------------------------------

function run_plot!(s::RunState, name, opts)
    isempty(methods(plot_probability)) && throw(ArgumentError(PLOTS_MESSAGE))
    m = s.model
    type = get(opts, :Type, :Probability)
    file = String(get(opts, :File, "$(name).pdf"))
    title = String(get(opts, :Title, ""))
    labels(list) =
        [(r.args[1] isa String ? r.args[1] : string(r.args[1])) => r.args[2] for r in list.args]
    plt = if type === :Probability
        setup = Setup()
        use_model!(m; setup)
        E = run_value(s, option(s, opts, :E, MExpr(:Range, Any[0.5, 10, 0.05])), opts; unit = :GeV)
        L = Float64(real(run_value(s, option(s, opts, :L, 1300), opts; unit = :km)))
        ρ = Float64(real(run_value(s, option(s, opts, :rho, 0.0), opts)))
        anti = option(s, opts, :Antineutrino, :False) === :True
        curves =
            haskey(opts, :Curves) ?
            [l => run_params(s, v) for (l, v) in labels(opts[:Curves])] :
            ["" => run_params(s, nothing)]
        bands = Pair{String, Any}[]
        if haskey(opts, :Band)
            for (l, v) in labels(opts[:Band])
                rules = v.args
                scan = only(
                    r for r in rules if is_call(r, :Rule) &&
                        !(r.args[2] isa Number) &&
                        is_call(r.args[2], :Range)
                )
                fixed = MExpr(:List, Any[r for r in rules if r !== scan])
                push!(
                    bands,
                    l => (
                        run_params(s, fixed),
                        key_of(m, scan.args[1]) =>
                            Float64.(real.(run_value(s, scan.args[2], opts))),
                    ),
                )
            end
        end
        with_setup(setup) do
            plot_probability(
                curves...;
                from = flavour_arg(option(s, opts, :From, :mu)),
                to = flavour_arg(option(s, opts, :To, :e)),
                E,
                L,
                ρ,
                antineutrino = anti,
                bands,
                title,
                setup,
            )
        end
    elseif type === :Fit
        fits = [
            f for f in s.fits if Symbol(f[1]) in
            (haskey(opts, :Fits) ? opts[:Fits].args : [Symbol(f[1]) for f in s.fits])
        ]
        isempty(fits) && throw(ArgumentError("Plot[$name]: no fits to plot"))
        x = fits[1][2].x
        ys = [f[1] => (haskey(f[2], :chi2) ? f[2].chi2 : f[2].sigma) for f in fits]
        plot_curves(
            x,
            ys...;
            xlabel = String(get(opts, :XLabel, string(fits[1][2].parameter))),
            ylabel = String(get(opts, :YLabel, haskey(fits[1][2], :chi2) ? "Δχ²" : "σ")),
            title,
        )
    else
        throw(ArgumentError("unknown plot type $type; use Probability or Fit"))
    end
    path = joinpath(s.output, file)
    savefig_(plt, path)
    push!(s.figures, (file, String(get(opts, :Caption, title))))
    println("Plot[$name] written to $file")
end

# Plots.savefig, provided by the extension
function savefig_ end

# ---------------------------------------------------------------------------
#  LaTeX
# ---------------------------------------------------------------------------

function run_latex!(s::RunState, file, opts)
    m = s.model
    file = file === nothing ? "$(m.name).tex" : String(file)
    on(k, default = true) = get(opts, k, default ? :On : :Off) in (:On, :True)
    r = TexReport(
        String(get(opts, :Title, "The $(m.name) model"));
        author = String(get(opts, :Author, m.authors)),
    )
    sections = Dict(k => on(k) for k in MODEL_SECTIONS)
    add!(r, latex(m; sections...))
    if on(:Results) && !isempty(s.calculations)
        section!(r, "Results")
        for (title, _, tex_) in s.calculations
            add!(r, escape_latex(title) * ":\n" * tex_ * "\n")
        end
    end
    if on(:Fits) && !isempty(s.fits)
        section!(r, "Fits")
        for (name, res) in s.fits
            y = haskey(res, :chi2) ? res.chi2 : res.sigma
            i = argmin(haskey(res, :chi2) ? y : -y)
            add!(
                r,
                "\\paragraph{$(escape_latex(name))} $(length(res.x)) points of \$$(guess_symbol_latex(res.parameter))\$; " *
                (
                    haskey(res, :chi2) ?
                    "\$\\Delta\\chi^2 = 0\$ at \$$(latex_number(res.x[i]))\$." :
                    "largest significance \$$(latex_number(y[i]))\\sigma\$ at \$$(latex_number(res.x[i]))\$."
                ) *
                "\n",
            )
        end
    end
    if on(:Figures) && !isempty(s.figures)
        section!(r, "Figures")
        for (f, caption) in s.figures
            figure!(r, f; caption = escape_latex(caption))
        end
    end
    text = get(opts, :Text, MExpr(:List, Any[]))
    for x in (is_call(text, :List) ? text.args : Any[text])
        if x isa String
            add!(r, x * "\n")
        elseif is_call(x, :Section)
            section!(r, String(x.args[1]))
        else
            throw(ArgumentError("Text of LaTeX[[i]] holds strings and Section[\"title\"]"))
        end
    end
    path = joinpath(s.output, file)
    write_tex(path, r; compile = get(opts, :Compile, :False) in (:On, :True))
    s.tex = path
    println("LaTeX written to $file")
end
