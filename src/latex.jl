# LaTeX output: expressions of model files, models, parameters, NuFIT
# tables, fit results, and whole reports with sections that can be switched
# on and off.

# ---------------------------------------------------------------------------
#  Expressions
# ---------------------------------------------------------------------------

const LATEX_FUNCTIONS = Dict(
    :Sin => "\\sin",
    :Cos => "\\cos",
    :Tan => "\\tan",
    :ArcSin => "\\arcsin",
    :ArcCos => "\\arccos",
    :ArcTan => "\\arctan",
    :Log => "\\ln",
    :Det => "\\det",
    :Tr => "\\operatorname{tr}",
    :Arg => "\\arg",
    :Re => "\\operatorname{Re}",
    :Im => "\\operatorname{Im}",
)

const GREEK_LETTERS = Set([
    "alpha",
    "beta",
    "gamma",
    "delta",
    "epsilon",
    "zeta",
    "eta",
    "theta",
    "iota",
    "kappa",
    "lambda",
    "mu",
    "nu",
    "xi",
    "pi",
    "rho",
    "sigma",
    "tau",
    "upsilon",
    "phi",
    "chi",
    "psi",
    "omega",
    "Gamma",
    "Delta",
    "Theta",
    "Lambda",
    "Xi",
    "Pi",
    "Sigma",
    "Phi",
    "Psi",
    "Omega",
])

"""
    latex_expression(expr; names=Dict()) -> String

LaTeX for a model-file expression. `names` maps symbols to their LaTeX
(parameter names, fields); other symbols are guessed (`theta12` → θ₁₂,
`eta_etau` → η_eτ).

```julia
latex_expression("Sqrt[Abs[dm31]] eta_etau")   # an expression in the Wolfram language
```
"""
latex_expression(e; names = Dict{Symbol, String}()) = tex(e, names, 0)

# Precedence: 0 sum, 1 product, 2 power base
function tex(e, names, prec)
    e isa Integer && return string(e)
    e isa Rational && return "\\frac{$(numerator(e))}{$(denominator(e))}"
    e isa AbstractFloat && return latex_number(e)
    e isa String && return "\\text{$(escape_latex(e))}"
    e isa Symbol && return symbol_latex(e, names)
    e isa MExpr || return string(e)
    h, a = e.head, e.args
    sub(x, p = 0) = tex(x, names, p)
    if h === :Plus
        out = sub(a[1])
        for x in a[2:end]
            neg, y = split_sign(x)
            out *= neg ? " - " * sub(y, 1) : " + " * sub(x)
        end
        return prec > 0 ? "\\left($out\\right)" : out
    elseif h === :Times
        neg, y = split_sign(e)
        if neg
            s = "-" * sub(y, 1)
            return prec > 0 ? "\\left($s\\right)" : s
        end
        num, den = Any[], Any[]
        for x in a
            if is_call(x, :Power) && x.args[2] isa Number && x.args[2] < 0
                push!(den, x.args[2] == -1 ? x.args[1] : MExpr(:Power, Any[x.args[1], -x.args[2]]))
            elseif x isa Rational
                numerator(x) == 1 || push!(num, numerator(x))
                push!(den, denominator(x))
            else
                push!(num, x)
            end
        end
        n = isempty(num) ? "1" : join((sub(x, 1) for x in num), " ")
        s = isempty(den) ? n : "\\frac{$n}{$(join((sub(x, 1) for x in den), " "))}"
        return prec > 1 ? "\\left($s\\right)" : s
    elseif h === :Dot
        s = join((sub(x, 1) for x in a), " \\, ")
        return prec > 1 ? "\\left($s\\right)" : s
    elseif h === :Power
        base, ex = a
        ex == 1 // 2 && return "\\sqrt{$(sub(base))}"
        ex == -1 && return "\\frac{1}{$(sub(base))}"
        return "{$(sub(base, 2))}^{$(sub(ex))}"
    elseif h === :List
        if !isempty(a) && all(x -> is_call(x, :List), a)
            rows = (join((sub(x) for x in r.args), " & ") for r in a)
            return "\\begin{pmatrix} " * join(rows, " \\\\ ") * " \\end{pmatrix}"
        end
        return "\\left(" * join((sub(x) for x in a), ", ") * "\\right)"
    elseif h === :Sqrt
        return "\\sqrt{$(sub(a[1]))}"
    elseif h === :Exp
        return "e^{$(sub(a[1]))}"
    elseif h === :Abs
        return "\\left|$(sub(a[1]))\\right|"
    elseif h in (:Dagger, :ConjugateTranspose, :HC, :conj)
        return "{$(sub(a[1], 2))}^\\dagger"
    elseif h in (:Conj, :Conjugate)
        return "{$(sub(a[1], 2))}^*"
    elseif h === :Transpose
        return "{$(sub(a[1], 2))}^T"
    elseif h === :Inverse
        return "{$(sub(a[1], 2))}^{-1}"
    elseif h in (:Diag, :DiagonalMatrix)
        items = is_call(a[1], :List) ? a[1].args : a
        return "\\operatorname{diag}\\left(" * join((sub(x) for x in items), ", ") * "\\right)"
    elseif h === :Rotation
        δ = length(a) > 3 ? ", " * sub(a[4]) : ""
        return "R_{$(a[1])$(a[2])}\\left($(sub(a[3]))$δ\\right)"
    elseif h === :PhaseMatrix
        return "P_{$(a[1])}\\left($(sub(a[2]))\\right)"
    elseif h === :IdentityMatrix
        return "\\mathbb{1}"
    elseif h === :Rule
        return sub(a[1]) * " \\to " * sub(a[2])
    elseif h === :Set
        return sub(a[1]) * " = " * sub(a[2])
    elseif h isa Symbol && haskey(LATEX_FUNCTIONS, h)
        return LATEX_FUNCTIONS[h] * "\\left(" * join((sub(x) for x in a), ", ") * "\\right)"
    end
    head = h isa Symbol ? "\\operatorname{$(escape_latex(String(h)))}" : sub(h)
    return head * "\\left(" * join((sub(x) for x in a), ", ") * "\\right)"
end

# (negative?, expression without the sign)
function split_sign(x)
    if is_call(x, :Times) && x.args[1] isa Number && x.args[1] < 0
        c = -x.args[1]
        rest = x.args[2:end]
        return true,
        c == 1 ? (length(rest) == 1 ? rest[1] : MExpr(:Times, rest)) : MExpr(:Times, Any[c; rest])
    end
    x isa Number && x < 0 && return true, -x
    return false, x
end

function symbol_latex(s::Symbol, names)
    haskey(names, s) && return names[s]
    s === :I && return "i"
    s === :Pi && return "\\pi"
    s === :E && return "E"
    s === :Degree && return "^\\circ"
    s in (:GeV, :MeV, :keV, :eV, :TeV, :km) && return "\\,\\mathrm{$s}"
    return guess_symbol_latex(s)
end

function guess_symbol_latex(s::Symbol)
    str = String(s)
    str in GREEK_LETTERS && return "\\" * str
    m = match(r"^([A-Za-z]+?)(\d+)$", str)
    if m !== nothing
        base = m.captures[1]
        return (base in GREEK_LETTERS ? "\\" * base : base) * "_{" * m.captures[2] * "}"
    end
    occursin('_', str) && return guess_latex(s)
    length(str) == 1 && return str
    any(g -> startswith(str, g) && length(str) > length(g), GREEK_LETTERS) && return guess_latex(s)
    return "\\mathrm{$str}"
end

escape_latex(s::AbstractString) =
    replace(s, r"([&%$#_{}])" => s"\\\1", "~" => "\\textasciitilde{}", "^" => "\\textasciicircum{}")

function latex_number(x::Real; digits = 4)
    x == 0 && return "0"
    isinteger(x) && abs(x) < 1e6 && return string(Int(x))
    e = floor(Int, log10(abs(x)))
    if -3 <= e <= 4
        return @sprintf("%.*g", digits, x)
    end
    mant = x / 10.0^e
    m = @sprintf("%.*g", digits, mant)
    return m == "1" ? "10^{$e}" : "$m \\times 10^{$e}"
end

latex_number(z::Complex; digits = 4) =
    imag(z) == 0 ? latex_number(real(z); digits) :
    real(z) == 0 ? latex_number(imag(z); digits) * "\\,i" :
    latex_number(real(z); digits) *
    (imag(z) < 0 ? " - " : " + ") *
    latex_number(abs(imag(z)); digits) *
    "\\,i"

"""
    latex(x) -> String

LaTeX for a GLoBES.jl object: a matrix or number, `OscParams`, `NuFit`,
`FitResult`, `OscillationModel` or a model-file expression.

```julia
latex(nufit("6.1"))   # NuFIT tables, fit results, OscParams, matrices and numbers
latex([1 2; 3 4])   # a matrix
```
"""
latex(x::Number; digits = 4) = latex_number(x; digits)
latex(e::MExpr; names = Dict{Symbol, String}()) = latex_expression(e; names)
latex(s::Symbol) = guess_symbol_latex(s)

function latex(M::AbstractMatrix; digits = 4, threshold = 1e-12)
    clean(z) = abs(z) < threshold ? zero(z) : z
    rows = (
        join((latex_number(clean(M[i, j]); digits) for j in axes(M, 2)), " & ") for i in axes(M, 1)
    )
    return "\\begin{pmatrix} " * join(rows, " \\\\ ") * " \\end{pmatrix}"
end
latex(v::AbstractVector{<:Number}; digits = 4) =
    "\\left(" * join((latex_number(x; digits) for x in v), ", ") * "\\right)"

function latex(p::OscParams)
    rows = [
        ("\\sin^2\\theta_{12}", latex_number(sin(p.θ12)^2)),
        ("\\sin^2\\theta_{13}", latex_number(sin(p.θ13)^2)),
        ("\\sin^2\\theta_{23}", latex_number(sin(p.θ23)^2)),
        ("\\delta_{\\mathrm{CP}}", latex_number(rad2deg(p.δ)) * "^\\circ"),
        ("\\Delta m^2_{21}", latex_number(p.Δm21) * "\\,\\mathrm{eV}^2"),
        ("\\Delta m^2_{31}", latex_number(p.Δm31) * "\\,\\mathrm{eV}^2"),
    ]
    for (i, n) in enumerate(p.extra_names)
        push!(rows, (guess_symbol_latex(n), latex_number(p.values[6 + i])))
    end
    return latex_table(["Parameter", "Value"], [["\$$a\$", "\$$b\$"] for (a, b) in rows])
end

function latex(nf::NuFit)
    v(t, scale = "") =
        "\$" *
        latex_number(t.best) *
        "^{+" *
        latex_number(t.plus) *
        "}_{-" *
        latex_number(t.minus) *
        "}\$"
    r(t) = "\$" * latex_number(t.lo) * " \\to " * latex_number(t.hi) * "\$"
    rows = [
        ["\$\\sin^2\\theta_{12}\$", v(nf.s12), r(nf.s12)],
        ["\$\\theta_{12}/^\\circ\$", v(nf.t12), r(nf.t12)],
        ["\$\\sin^2\\theta_{23}\$", v(nf.s23), r(nf.s23)],
        ["\$\\theta_{23}/^\\circ\$", v(nf.t23), r(nf.t23)],
        ["\$\\sin^2\\theta_{13}\$", v(nf.s13), r(nf.s13)],
        ["\$\\theta_{13}/^\\circ\$", v(nf.t13), r(nf.t13)],
        ["\$\\delta_{\\mathrm{CP}}/^\\circ\$", v(nf.dcp), r(nf.dcp)],
        ["\$\\Delta m^2_{21}/10^{-5}\\,\\mathrm{eV}^2\$", v(nf.dm21), r(nf.dm21)],
        [
            nf.ordering === :normal ? "\$\\Delta m^2_{31}/10^{-3}\\,\\mathrm{eV}^2\$" :
            "\$\\Delta m^2_{32}/10^{-3}\\,\\mathrm{eV}^2\$",
            v(nf.dm3l),
            r(nf.dm3l),
        ],
    ]
    caption = "NuFIT $(nf.version), $(nf.ordering) ordering (\$\\Delta\\chi^2 = $(nf.chi2)\$), $(dataset(nf))."
    return latex_table(["", "best fit \$\\pm 1\\sigma\$", "\$3\\sigma\$ range"], rows; caption)
end

function latex(r::FitResult)
    head = "\\chi^2_{\\min} = $(latex_number(r.chi2))"
    free =
        isempty(r.free) ? "" :
        " \\quad \\text{(free: \$" * join(guess_symbol_latex.(r.free), ", ") * "\$)}"
    return "\\[ $head$free \\]\n" * latex(r.params)
end

"""
    latex_table(header, rows; caption="", label="") -> String

A booktabs table; cells are LaTeX strings.

```julia
latex_table(["Parameter", "Value"], [["θ12", "33.7°"], ["θ13", "8.6°"]]; caption = "Best fit")   # header, rows, caption
```
"""
function latex_table(header, rows; caption = "", label = "")
    io = IOBuffer()
    println(io, "\\begin{table}[htbp]\n\\centering")
    println(io, "\\begin{tabular}{", "l"^length(header), "}\n\\toprule")
    println(io, join(header, " & "), " \\\\\n\\midrule")
    for r in rows
        println(io, join(r, " & "), " \\\\")
    end
    println(io, "\\bottomrule\n\\end{tabular}")
    isempty(caption) || println(io, "\\caption{", caption, "}")
    isempty(label) || println(io, "\\label{", label, "}")
    print(io, "\\end{table}\n")
    return String(take!(io))
end

# ---------------------------------------------------------------------------
#  Oscillation models
# ---------------------------------------------------------------------------

# LaTeX names of the symbols of a model
function latex_names(m::OscillationModel)
    names = Dict{Symbol, String}(p.name => p.latex for p in m.parameters)
    names[:deltaCP] = "\\delta_{\\mathrm{CP}}"
    names[:U] = "U"
    names[:M2] = "M^2"
    names[:M] = "M"
    names[:A] = "A"
    names[:PMNS] = "U_{\\mathrm{PMNS}}"
    names[:Ud] = "U_\\delta"
    names[:Gamma] = "\\Gamma"
    for b in m.blocks
        b.kind === :AdditionalTerms && (names[b.name] = "H_{\\mathrm{$(b.name)}}")
    end
    return names
end

const MODEL_SECTIONS = (
    :Model,
    :Parameters,
    :NuFit,
    :Definitions,
    :MixingMatrix,
    :MassSquared,
    :MassHierarchy,
    :MassMatrix,
    :AdditionalTerms,
    :Decoherence,
    :Hamiltonian,
)

"""
    latex(model::OscillationModel; sections...) -> String

The model in LaTeX: one subsection per block, the composed Hamiltonian and
the parameters. Switch parts off with `Parameters=false` etc.; the parts
are $(join(MODEL_SECTIONS, ", ")).

```julia
latex(load_model("examples/models/SNSI"); Parameters = false, NuFit = false)   # model; parts switched off by name
```
"""
function latex(m::OscillationModel; kwargs...)
    on(k) = get(kwargs, k, true)
    names = latex_names(m)
    tx(e) = latex_expression(e; names)
    io = IOBuffer()
    if on(:Model)
        println(
            io,
            "\\section{The ",
            escape_latex(isempty(m.latex_name) ? m.name : m.latex_name),
            " model}",
        )
        isempty(m.authors) || println(io, "Model by ", escape_latex(m.authors), ".\n")
        println(
            io,
            "Flavours: \$",
            join((flavour_latex(f) for f in m.flavours), ", "),
            "\$ (",
            nflavours(m),
            " states). Energies, masses and the Hamiltonian are in eV; \$A = \\sqrt{2} G_F N_e\$ is the matter potential.\n",
        )
    end
    if on(:Hamiltonian)
        println(
            io,
            "\\subsection*{Effective Hamiltonian}\n\\begin{equation}\n  H = ",
            hamiltonian_latex(m, names),
            "\n\\end{equation}",
        )
        println(io, "For antineutrinos \$H \\to H^*\$ with \$A \\to -A\$.\n")
    end
    if on(:Definitions)
        defs = [k for k in m.global_order if any(b -> !b.off && block_uses(b, k, m), m.blocks)]
        if !isempty(defs)
            println(io, "\\subsection*{Definitions}\n\\begin{align}")
            println(
                io,
                join(
                    ("  " * symbol_latex(k, names) * " &= " * tx(m.globals[k]) for k in defs),
                    " \\\\\n",
                ),
            )
            println(io, "\\end{align}")
        end
    end
    for kind in
        (:MixingMatrix, :MassSquared, :MassHierarchy, :MassMatrix, :AdditionalTerms, :Decoherence)
        on(kind) || continue
        for b in m.blocks
            b.kind === kind && !b.off || continue
            kind === :MassHierarchy &&
                b.name !== active_hierarchy_name(m) &&
                m.hierarchy !== :Auto &&
                continue
            println(io, "\\subsection*{", block_title(b), b.default ? " (default)" : "", "}")
            lines = String[]
            main = main_symbol(kind)
            if b.definition !== nothing
                lhs =
                    kind === :AdditionalTerms ? names[b.name] :
                    kind === :Decoherence ? "\\Gamma_{ij}" :
                    symbol_latex(something(main, :X), names)
                push!(lines, lhs * " &= " * tx(b.definition))
            end
            for k in b.order
                k === main && continue
                push!(lines, symbol_latex(k, names) * " &= " * tx(b.locals[k]))
            end
            isempty(lines) ||
                println(io, "\\begin{align}\n  ", join(lines, " \\\\\n  "), "\n\\end{align}")
        end
    end
    on(:Parameters) && print(io, latex_parameters(m))
    on(:NuFit) && print(io, latex(m.nufit))
    return String(take!(io))
end

flavour_latex(f::Symbol) =
    f === :mu ? "\\nu_\\mu" : f === :tau ? "\\nu_\\tau" : f === :e ? "\\nu_e" : "\\nu_{$f}"

active_hierarchy_name(m) = m.hierarchy

function block_uses(b::ModelBlock, k::Symbol, m)
    found = false
    walk(e) = e === k ? (found = true) : e isa MExpr && (walk(e.head); foreach(walk, e.args))
    b.definition === nothing || walk(b.definition)
    foreach(walk, values(b.locals))
    # uses through other globals
    if !found
        for g in m.global_order
            g === k && continue
            uses_g = false
            w2(e) = e === g ? (uses_g = true) : e isa MExpr && (w2(e.head); foreach(w2, e.args))
            b.definition === nothing || w2(b.definition)
            if uses_g
                inner(e) =
                    e === k ? (found = true) :
                    e isa MExpr && (inner(e.head); foreach(inner, e.args))
                inner(m.globals[g])
            end
        end
    end
    return found
end

function block_title(b::ModelBlock)
    title = Dict(
        :MixingMatrix => "Mixing matrix",
        :MassSquared => "Masses squared",
        :MassHierarchy => "Mass hierarchy",
        :MassMatrix => "Mass matrix",
        :AdditionalTerms => "Additional term",
        :Decoherence => "Decoherence",
        :Hamiltonian => "Hamiltonian",
    )[b.kind]
    return b.kind in LIST_BLOCKS ? title * " " * String(b.name) : title
end

function hamiltonian_latex(m::OscillationModel, names)
    hb = findfirst(b -> b.kind === :Hamiltonian && !b.off, m.blocks)
    hb === nothing || return latex_expression(m.blocks[hb].definition; names)
    vac =
        has_block(m, :MassMatrix) ? "\\frac{1}{2E} M M^\\dagger" : "\\frac{1}{2E} U M^2 U^\\dagger"
    terms = [names[b.name] for b in m.blocks if b.kind === :AdditionalTerms && !b.off]
    return join([vac; terms], " + ")
end

function latex_parameters(m::OscillationModel)
    rows = Vector{String}[]
    for p in m.parameters
        unit = isempty(p.unit) ? "" : p.unit
        shown = p.value
        scale = get(MODEL_CONSTANTS, Symbol(unit), 1.0)
        unit in ("°",) && (shown = rad2deg(p.value))
        unit in ("GeV", "MeV", "keV", "TeV") && (shown = p.value / scale)
        src =
            p.source === :nufit ? "NuFIT $(m.nufit.version)" :
            p.source === :file ? "parameters.m" : "default"
        range =
            p.range === nothing ? "" :
            "\$[" * latex_number(p.range[1]) * ", " * latex_number(p.range[2]) * "]\$"
        unit_tex = unit == "°" ? "\$^\\circ\$" : unit == "eV²" ? "eV\$^2\$" : escape_latex(unit)
        push!(rows, ["\$$(p.latex)\$", "\$" * latex_number(shown) * "\$", unit_tex, range, src])
    end
    return latex_table(
        ["Parameter", "Value", "Unit", "Range", "Source"],
        rows;
        caption = "Parameters of the $(escape_latex(m.name)) model.",
    )
end

# ---------------------------------------------------------------------------
#  Reports
# ---------------------------------------------------------------------------

"""
    TexReport(title; author="", date="\\today")

A LaTeX document assembled from parts: text, sections, objects
([`latex`](@ref) of models, parameters, fit results, NuFIT tables), tables
and figures. Write it with [`write_tex`](@ref).

```julia
r = TexReport("SNSI at DUNE"; author = "A. Physicist")   # title, author, date
add!(r, load_model("examples/models/SNSI"); Parameters = false)   # report, object, options of latex
section!(r, "Results")   # report, section title
figure!(r, "pmue.pdf"; caption = "P(ν_μ → ν_e) at 1300 km")   # report, file, caption
write_tex("report.tex", r; compile = true)   # file, report, run pdflatex
```
"""
mutable struct TexReport
    title::String
    author::String
    date::String
    parts::Vector{String}
    packages::Vector{String}
end

TexReport(title::AbstractString = "GLoBES.jl results"; author = "", date = "\\today") = TexReport(
    String(title),
    String(author),
    String(date),
    String[],
    ["amsmath", "amssymb", "booktabs", "graphicx", "hyperref"],
)

"""
    add!(report, x; options...)

Append LaTeX to a report: a string is added as it is; other objects as
`latex(x; options...)`.

```julia
add!(r, "Text of a paragraph.")   # report, LaTeX text
add!(r, nufit("6.1"))   # report, any object that latex knows
```
"""
add!(r::TexReport, s::AbstractString) = (push!(r.parts, String(s)); r)
add!(r::TexReport, x; kwargs...) = (push!(r.parts, latex(x; kwargs...)); r)

"""
    section!(report, title; level=1)

```julia
section!(r, "Results"; level = 2)   # report, title, 1 for a section and 2 for a subsection
```
"""
section!(r::TexReport, title::AbstractString; level::Integer = 1) =
    add!(r, "\\" * ("sub"^(level - 1)) * "section{" * title * "}\n")

"""
    figure!(report, file; caption="", width="0.8\\textwidth")

```julia
figure!(r, "pmue.pdf"; caption = "Appearance probability")   # report, file, caption
```
"""
function figure!(
    r::TexReport,
    file::AbstractString;
    caption = "",
    width = "0.8\\textwidth",
    label = "",
)
    s = "\\begin{figure}[htbp]\n\\centering\n\\includegraphics[width=$width]{$file}\n"
    isempty(caption) || (s *= "\\caption{$caption}\n")
    isempty(label) || (s *= "\\label{$label}\n")
    return add!(r, s * "\\end{figure}\n")
end

"""
    table!(report, header, rows; caption="")

A table of strings or numbers.

```julia
table!(r, ["δ_CP", "σ"], [["0", "0.1"], ["π/2", "5.2"]]; caption = "CP violation")   # report, header, rows, caption
```
"""
table!(r::TexReport, header, rows; caption = "") = add!(
    r,
    latex_table(
        header,
        [[x isa Number ? "\$" * latex_number(x) * "\$" : string(x) for x in row] for row in rows];
        caption,
    ),
)

"""
    write_tex(file, report; compile=false) -> file

Write the report as a complete LaTeX document. `compile=true` runs
`latexmk` or `pdflatex` if one is installed.

```julia
write_tex("report.tex", r; compile = true)   # file, report, run pdflatex
```
"""
function write_tex(file::AbstractString, r::TexReport; compile::Bool = false)
    open(file, "w") do io
        println(io, "% Written by GLoBES.jl")
        println(io, "\\documentclass[11pt]{article}")
        println(io, "\\usepackage[margin=2.5cm]{geometry}")
        for p in r.packages
            println(io, "\\usepackage{", p, "}")
        end
        println(io, "\\title{", r.title, "}")
        isempty(r.author) || println(io, "\\author{", r.author, "}")
        println(io, "\\date{", r.date, "}")
        println(io, "\\begin{document}\n\\maketitle\n")
        for p in r.parts
            println(io, p)
        end
        println(io, "\\end{document}")
    end
    compile && compile_tex(file)
    return file
end

"""
    compile_tex(file) -> Bool

Run `latexmk -pdf` (or `pdflatex` twice) on a LaTeX file in its directory.
Returns `false` if no LaTeX is installed.

```julia
compile_tex("report.tex")   # file: true if pdflatex succeeded
```
"""
function compile_tex(file::AbstractString)
    dir, name = dirname(abspath(file)), basename(file)
    cmd = if Sys.which("latexmk") !== nothing
        `latexmk -pdf -interaction=nonstopmode -halt-on-error $name`
    elseif Sys.which("pdflatex") !== nothing
        `pdflatex -interaction=nonstopmode -halt-on-error $name`
    else
        @warn "no LaTeX installation found; $file was written but not compiled"
        return false
    end
    cd(dir) do
        run(pipeline(cmd; stdout = devnull, stderr = devnull))
        cmd.exec[1] == "pdflatex" && run(pipeline(cmd; stdout = devnull, stderr = devnull))
    end
    return true
end
