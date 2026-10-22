# AEDL ("Abstract Experiment Definition Language") reader.
#
# Replaces glb_lexer.l / glb_parser.y with a hand-written lexer and a
# recursive-descent parser. The semantics follow the bison grammar: groups
# `rule(#name)< ... >` etc., `$global` and `@local` directives, user variables,
# `%lists`, `include "file"`, `#NEXT#` separators and the built-in functions.

struct AEDLError <: Exception
    msg::String
end
Base.showerror(io::IO, e::AEDLError) = print(io, "AEDL error: ", e.msg)

struct Token
    kind::Symbol      # :num :str :ident :dollar :at :hash :pct :op :rulemult :nl :next :eof
    text::String
    val::Float64
    file::String
    line::Int
    dir::String       # directory of the AEDL file the token comes from
end

const MAX_INCLUDE_DEPTH = 10

isidstart(c) = ('a' <= c <= 'z') || ('A' <= c <= 'Z')
isidchar(c) = isidstart(c) || isdigit(c) || c == '_'

function lex_number(s, i, n)
    j = i
    while j <= n && isdigit(s[j])
        j += 1
    end
    if j <= n && s[j] == '.'
        j += 1
        while j <= n && isdigit(s[j])
            j += 1
        end
    end
    if j <= n && (s[j] == 'e' || s[j] == 'E')
        k = j + 1
        k <= n && (s[k] == '+' || s[k] == '-') && (k += 1)
        if k <= n && isdigit(s[k])
            while k <= n && isdigit(s[k])
                k += 1
            end
            j = k
        end
    end
    return j
end

"""
    aedl_tokenize(path; depth=0) -> Vector{Token}

Tokenize an AEDL file, splicing in `include "file"` directives.
"""
function aedl_tokenize(path::AbstractString; depth::Int = 0, display_name::AbstractString = path)
    depth > MAX_INCLUDE_DEPTH && throw(AEDLError("Includes nested too deeply"))
    s = collect(read(path, String))
    n = length(s)
    toks = Token[]
    line = Ref(1)    # a Ref, since `tok` captures it
    i = 1
    dir = dirname(abspath(path))
    tok(kind, text, val = 0.0) = push!(toks, Token(kind, text, val, display_name, line[], dir))
    while i <= n
        c = s[i]
        if c == '\n'
            tok(:nl, "\n")
            line[] += 1
            i += 1
        elseif c == ' ' || c == '\t' || c == '\r'
            i += 1
        elseif c == '/' && i < n && s[i + 1] == '/'
            while i <= n && s[i] != '\n'
                i += 1
            end
        elseif c == '/' && i < n && s[i + 1] == '*'
            i += 2
            while i <= n && !(s[i] == '*' && i < n && s[i + 1] == '/')
                s[i] == '\n' && (line[] += 1)
                i += 1
            end
            i += 2
        elseif c == '%' && i + 7 <= n && String(s[i:(i + 7)]) == "%!GLoBES"
            i += 8
        elseif c == '"'
            j = i + 1
            while j <= n && s[j] != '"'
                s[j] == '\n' && throw(AEDLError("Unterminated string in $display_name:$(line[])"))
                j += 1
            end
            j > n && throw(AEDLError("Unterminated string in $display_name:$(line[])"))
            str = String(s[(i + 1):(j - 1)])
            i = j + 1
            if !isempty(toks) && toks[end].kind == :ident && toks[end].text == "include"
                pop!(toks)
                incpath = glb_find_file(str, dirname(path))
                incpath === nothing && throw(
                    AEDLError("Could not open include file \"$str\" ($display_name:$(line[]))"),
                )
                inc = aedl_tokenize(incpath; depth = depth + 1, display_name = str)
                inc[end].kind == :eof && pop!(inc)
                append!(toks, inc)
            else
                tok(:str, str)
            end
        elseif isdigit(c) || (c == '.' && i < n && isdigit(s[i + 1]))
            j = lex_number(s, i, n)
            txt = String(s[i:(j - 1)])
            tok(:num, txt, parse(Float64, txt))
            i = j
        elseif isidstart(c)
            j = i
            while j <= n && isidchar(s[j])
                j += 1
            end
            tok(:ident, String(s[i:(j - 1)]))
            i = j
        elseif (c == '$' || c == '@' || c == '%') && i < n && isidstart(s[i + 1])
            j = i + 1
            while j <= n && isidchar(s[j])
                j += 1
            end
            tok(c == '$' ? :dollar : c == '@' ? :at : :pct, String(s[i:(j - 1)]))
            i = j
        elseif c == '@'
            tok(:rulemult, "@")
            i += 1
        elseif c == '#'
            if i + 5 <= n && String(s[i:(i + 5)]) == "#NEXT#"
                tok(:next, "#NEXT#")
                i += 6
            elseif i < n && isidstart(s[i + 1])
                j = i + 1
                while j <= n && isidchar(s[j])
                    j += 1
                end
                tok(:hash, String(s[i:(j - 1)]))
                i = j
            else
                throw(AEDLError("Unexpected '#' in $display_name:$(line[])"))
            end
        else
            tok(:op, string(c))
            i += 1
        end
    end
    tok(:eof, "")
    return toks
end

# ---------------------------------------------------------------------------
#  Parser state
# ---------------------------------------------------------------------------

const GROUP_NAMES = ("rule", "channel", "energy", "cross", "flux", "nuflux")

mutable struct AEDLParser
    toks::Vector{Token}
    pos::Int
    exp::ExperimentData
    finished::Vector{ExperimentData}
    context::String
    loc_count::Int
    vars::Dict{String, Float64}
    lists::Dict{String, Vector{Float64}}
    names::Dict{String, Tuple{String, Int}}     # name -> (context, value)
    counters::Dict{String, Int}
    flt::Flux
    xsc::XSec
    ibf::SmearData
    filename::String
end

# User-predefined AEDL variables and lists (glbDefineAEDLVariable/List)
const aedl_presym_vars = Dict{String, Float64}()
const aedl_presym_lists = Dict{String, Vector{Float64}}()

function new_parser(toks, filename)
    p = AEDLParser(
        toks,
        1,
        ExperimentData(),
        ExperimentData[],
        "global",
        -1,
        Dict{String, Float64}(),
        Dict{String, Vector{Float64}}(),
        Dict{String, Tuple{String, Int}}(),
        Dict{String, Int}(),
        Flux(),
        XSec(),
        SmearData(),
        filename,
    )
    reset_parser_tables!(p)
    return p
end

function reset_parser_tables!(p::AEDLParser)
    p.vars = copy(aedl_presym_vars)
    p.lists = Dict(k => copy(v) for (k, v) in aedl_presym_lists)
    p.names = Dict{String, Tuple{String, Int}}()
    for k in keys(SIGMA_FUNCTIONS)
        p.names[k] = ("energy", -1)
    end
    p.counters = Dict("rule" => 0, "channel" => 0, "energy" => 0, "cross" => 0, "flux" => 0)
    p.exp = ExperimentData()
    p.exp.filename = p.filename
    p.context = "global"
    p.loc_count = -1
end

peek(p::AEDLParser) = p.toks[p.pos]
peekk(p::AEDLParser, k = 1) = p.toks[min(p.pos + k, length(p.toks))]
advance!(p::AEDLParser) = (t = p.toks[p.pos]; p.pos < length(p.toks) && (p.pos += 1); t)

function perr(p::AEDLParser, msg, t::Token = peek(p))
    throw(AEDLError("$msg in file [$(t.file)]: in line $(t.line)"))
end

isop(t::Token, c) = t.kind == :op && t.text == c

function expect_op!(p::AEDLParser, c)
    t = peek(p)
    isop(t, c) || perr(p, "Expected '$c' but found '$(t.text)'")
    advance!(p)
end

skipnl!(p::AEDLParser) =
    while peek(p).kind == :nl
        advance!(p)
    end

# ---------------------------------------------------------------------------
#  Expressions and lists
# ---------------------------------------------------------------------------

const SCALAR_FUNCTIONS = Dict{String, Function}(
    "sin" => sin,
    "cos" => cos,
    "tan" => tan,
    "asin" => asin,
    "acos" => acos,
    "atan" => atan,
    "log" => log,
    "exp" => exp,
    "log10" => log10,
    "sqrt" => sqrt,
)
const LIST_FUNCTIONS =
    ("reverse", "copy", "echo", "echon", "line", "interpolation", "bincenter", "samplingbincenter")

aedl_isfunction(name) = haskey(SCALAR_FUNCTIONS, name) || name in LIST_FUNCTIONS

function parse_expr(p::AEDLParser)
    v = parse_term(p)
    while true
        t = peek(p)
        if isop(t, "+") || isop(t, "-")
            advance!(p)
            skipnl!(p)
            r = parse_term(p)
            v = scalar(p, v)
            r = scalar(p, r)
            v = t.text == "+" ? v + r : v - r
        else
            return v
        end
    end
end

function parse_term(p::AEDLParser)
    v = parse_unary(p)
    while true
        t = peek(p)
        if isop(t, "*") || isop(t, "/")
            advance!(p)
            skipnl!(p)
            r = scalar(p, parse_unary(p))
            v = scalar(p, v)
            v = t.text == "*" ? v * r : v / r
        else
            return v
        end
    end
end

function parse_unary(p::AEDLParser)
    if isop(peek(p), "-")
        advance!(p)
        return -scalar(p, parse_unary(p))
    elseif isop(peek(p), "+")
        advance!(p)
        return parse_unary(p)
    end
    return parse_power(p)
end

function parse_power(p::AEDLParser)
    v = parse_primary(p)
    if isop(peek(p), "^")
        advance!(p)
        e = scalar(p, parse_unary(p))
        return scalar(p, v)^e
    end
    return v
end

scalar(p::AEDLParser, v::Float64) = v
scalar(p::AEDLParser, v) = perr(p, "List used where a number is expected")

function parse_list_literal(p::AEDLParser)
    expect_op!(p, "{")
    out = Float64[]
    skipnl!(p)
    while !isop(peek(p), "}")
        peek(p).kind == :eof && perr(p, "Unterminated list")
        push!(out, scalar(p, parse_expr(p)))
        skipnl!(p)
        if isop(peek(p), ",")
            advance!(p)
            skipnl!(p)
        end
    end
    advance!(p)
    return out
end

function parse_primary(p::AEDLParser)
    t = peek(p)
    if t.kind == :num
        advance!(p)
        return t.val
    elseif isop(t, "(")
        advance!(p)
        skipnl!(p)
        v = parse_expr(p)
        skipnl!(p)
        expect_op!(p, ")")
        return v
    elseif isop(t, "{")
        return parse_list_literal(p)
    elseif t.kind == :pct
        advance!(p)
        haskey(p.lists, t.text) || perr(p, "Undefined list $(t.text)")
        return p.lists[t.text]
    elseif t.kind == :hash
        advance!(p)
        haskey(p.names, t.text) || perr(p, "Unknown name $(t.text)")
        return Float64(p.names[t.text][2])
    elseif t.kind == :ident
        if isop(peekk(p), "(") && aedl_isfunction(t.text)
            return parse_function_call(p)
        end
        advance!(p)
        haskey(p.vars, t.text) || perr(p, "Undefined variable $(t.text)")
        return p.vars[t.text]
    end
    perr(p, "Syntax error near '$(t.text)'")
end

function parse_function_call(p::AEDLParser)
    name = advance!(p).text
    expect_op!(p, "(")
    skipnl!(p)
    if isop(peek(p), ")")
        advance!(p)
        name == "bincenter" && return aedl_bincenter(p)
        name == "samplingbincenter" && return aedl_samplingbincenter(p)
        perr(p, "Function $name needs an argument")
    end
    args = Any[parse_expr(p)]
    skipnl!(p)
    while isop(peek(p), ",")
        advance!(p)
        skipnl!(p)
        push!(args, parse_expr(p))
        skipnl!(p)
    end
    expect_op!(p, ")")
    if name == "interpolation"
        length(args) == 4 || perr(p, "interpolation needs 4 arguments")
        return aedl_interpolation(p, args[1], args[2], floor(Int, scalar(p, args[3])), args[4])
    end
    length(args) == 1 || perr(p, "Function $name takes one argument")
    a = args[1]
    if haskey(SCALAR_FUNCTIONS, name)
        f = SCALAR_FUNCTIONS[name]
        return a isa Float64 ? f(a) : map(f, a)
    elseif name == "reverse"
        return a isa Float64 ? a : reverse(a)
    elseif name == "copy"
        return copy(a)
    elseif name == "echo"
        foreach(x -> print(x, " "), a isa Float64 ? [a] : a)
        return a
    elseif name == "echon"
        foreach(println, a isa Float64 ? [a] : a)
        return a
    elseif name == "line"
        for x in (a isa Float64 ? [a] : a)
            print("\n"^floor(Int, x))
        end
        return a
    end
    perr(p, "Improper use of special function $name")
end

function aedl_bincenter(p::AEDLParser)
    p.exp.numofbins < 0 && perr(p, "Cannot compute bincenter. Binning not set up properly")
    sd = SmearData()
    set_up_smear_data!(sd, p.exp)
    return copy(sd.bincenter)
end

function aedl_samplingbincenter(p::AEDLParser)
    p.exp.simbins < 0 &&
        perr(p, "Cannot compute samplingbincenter. Sampling-binning not set up properly")
    sd = SmearData()
    set_up_smear_data!(sd, p.exp)
    return copy(sd.simbincenter)
end

# Linear (flag 1) or natural cubic spline (flag 2) interpolation, as GSL does
function aedl_interpolation(p::AEDLParser, xs, ys, flag, where)
    (xs isa Vector && ys isa Vector) || perr(p, "interpolation expects lists")
    length(xs) == length(ys) || perr(p, "Xval and Yval in interpolation are not of the same length")
    w = where isa Float64 ? [where] : where
    if flag == 1
        return [interp_linear(p, xs, ys, x) for x in w]
    elseif flag == 2
        m = natural_spline_moments(xs, ys)
        return [interp_cspline(p, xs, ys, m, x) for x in w]
    end
    perr(p, "Invalid interpolation scheme flag")
end

function interp_index(p, xs, x)
    (x < xs[1] || x > xs[end]) && perr(p, "interpolation: point $x outside of data range")
    i = searchsortedlast(xs, x)
    return clamp(i, 1, length(xs) - 1)
end

function interp_linear(p, xs, ys, x)
    i = interp_index(p, xs, x)
    dx = xs[i + 1] - xs[i]
    return ys[i] + (ys[i + 1] - ys[i]) * (x - xs[i]) / dx
end

function natural_spline_moments(xs, ys)
    n = length(xs)
    M = zeros(n)
    n < 3 && return M
    h = diff(xs)
    a = zeros(n - 2)
    b = zeros(n - 2)
    c = zeros(n - 2)
    r = zeros(n - 2)
    for i = 2:(n - 1)
        k = i - 1
        a[k] = h[i - 1]
        b[k] = 2 * (h[i - 1] + h[i])
        c[k] = h[i]
        r[k] = 6 * ((ys[i + 1] - ys[i]) / h[i] - (ys[i] - ys[i - 1]) / h[i - 1])
    end
    # Thomas algorithm
    for k = 2:(n - 2)
        w = a[k] / b[k - 1]
        b[k] -= w * c[k - 1]
        r[k] -= w * r[k - 1]
    end
    sol = zeros(n - 2)
    sol[end] = r[end] / b[end]
    for k = (n - 3):-1:1
        sol[k] = (r[k] - c[k] * sol[k + 1]) / b[k]
    end
    M[2:(n - 1)] .= sol
    return M
end

function interp_cspline(p, xs, ys, M, x)
    i = interp_index(p, xs, x)
    h = xs[i + 1] - xs[i]
    A = (xs[i + 1] - x) / h
    B = (x - xs[i]) / h
    return A * ys[i] + B * ys[i + 1] + ((A^3 - A) * M[i] + (B^3 - B) * M[i + 1]) * h^2 / 6
end

# A right-hand side: number, list, `a : b` pair, or comma-separated list
function parse_rhs(p::AEDLParser)
    v = parse_expr(p)
    if v isa Float64
        if isop(peek(p), ":")
            advance!(p)
            return (v, scalar(p, parse_expr(p)))
        elseif isop(peek(p), ",")
            out = [v]
            while isop(peek(p), ",")
                advance!(p)
                skipnl!(p)
                (peek(p).kind == :nl || peek(p).kind == :eof) && break
                push!(out, scalar(p, parse_expr(p)))
            end
            return out
        end
    end
    return v
end

function expect_string!(p::AEDLParser)
    t = peek(p)
    t.kind == :str || perr(p, "Expected a quoted string")
    advance!(p)
    return t.text
end

# A data file named in an AEDL file: like include files, it is looked up in
# the current directory, then next to the AEDL file, then in GLB_PATH.
function expect_data_file!(p::AEDLParser)
    t = peek(p)
    name = expect_string!(p)
    path = glb_find_file(name, t.dir)
    return path === nothing ? name : path
end

# ---------------------------------------------------------------------------
#  Statements
# ---------------------------------------------------------------------------

function parse_aedl!(p::AEDLParser)
    while true
        skipnl!(p)
        t = peek(p)
        if t.kind == :eof
            finish_experiment!(p)
            return p.finished
        elseif t.kind == :next
            advance!(p)
            finish_experiment!(p)
            reset_parser_tables!(p)
        else
            parse_statement!(p)
        end
    end
end

function finish_experiment!(p::AEDLParser)
    p.exp.names = [(k, v[1], v[2]) for (k, v) in p.names if v[2] != -1]
    push!(p.finished, p.exp)
end

function end_of_statement!(p::AEDLParser)
    t = peek(p)
    if t.kind == :nl || t.kind == :eof || t.kind == :next || isop(t, ">") || isop(t, ";")
        isop(t, ";") && advance!(p)
        return
    end
    # Several directives on one line are fine as long as a new one starts here
    (t.kind == :at || t.kind == :dollar) && return
    perr(p, "Unexpected '$(t.text)'")
end

function parse_statement!(p::AEDLParser)
    t = peek(p)
    if t.kind == :ident && t.text in GROUP_NAMES && isop(peekk(p), "(") && p.context == "global"
        parse_group!(p)
    elseif t.kind == :dollar
        advance!(p)
        expect_op!(p, "=")
        if t.text == "\$version"
            p.exp.version = expect_string!(p)
        else
            set_global!(p, t, parse_rhs(p))
        end
        end_of_statement!(p)
    elseif t.kind == :at
        p.context == "global" && perr(p, "Directive $(t.text) outside of a group")
        advance!(p)
        expect_op!(p, "=")
        parse_local!(p, t)
        end_of_statement!(p)
    elseif t.kind == :pct
        advance!(p)
        expect_op!(p, "=")
        v = parse_rhs(p)
        p.lists[t.text] = v isa Float64 ? [v] : v isa Tuple ? collect(v) : copy(v)
        end_of_statement!(p)
    elseif t.kind == :ident && isop(peekk(p), "=") && !aedl_isfunction(t.text)
        advance!(p)
        advance!(p)
        p.vars[t.text] = scalar(p, parse_expr(p))
        end_of_statement!(p)
    else
        parse_expr(p)   # e.g. echo(...)
        end_of_statement!(p)
    end
end

function parse_group!(p::AEDLParser)
    gid = advance!(p).text
    expect_op!(p, "(")
    nt = peek(p)
    nt.kind == :hash || perr(p, "Expected #name in $gid definition")
    advance!(p)
    expect_op!(p, ")")
    counter = gid == "nuflux" ? "flux" : gid
    ctx, val = get(p.names, nt.text, (gid, -1))
    ctx != gid && perr(p, "Redefinition of an automatic variable $(nt.text)")
    if val == -1
        p.counters[counter] += 1
        val = p.counters[counter]
    end
    p.names[nt.text] = (gid, val)
    p.loc_count = val
    p.context = gid
    group_start!(p, gid)
    skipnl!(p)
    expect_op!(p, "<")
    while true
        skipnl!(p)
        t = peek(p)
        t.kind == :eof && perr(p, "Unterminated $gid group")
        if isop(t, ">")
            advance!(p)
            break
        end
        parse_statement!(p)
    end
    group_end!(p, gid)
    p.context = "global"
end

function ensure_rule!(p::AEDLParser, k)
    while length(p.exp.rules) < k
        push!(p.exp.rules, Rule())
    end
    return p.exp.rules[k]
end

function ensure_channel!(p::AEDLParser, k)
    e = p.exp
    while length(e.channels) < k
        push!(e.channels, Channel(0, 0, 0, 0, 0, 0))
        push!(e.user_pre_smearing_channel, nothing)
        push!(e.user_post_smearing_channel, nothing)
        push!(e.user_pre_smearing_background, nothing)
        push!(e.user_post_smearing_background, nothing)
    end
    return e.channels[k]
end

function ensure_smear!(p::AEDLParser, k)
    e = p.exp
    while length(e.smear_data) < k
        push!(e.smear_data, SmearData())
        push!(e.smear, Vector{Float64}[])
        push!(e.lowrange, Int[])
        push!(e.uprange, Int[])
        push!(e.user_smear_given, false)
    end
end

function group_start!(p::AEDLParser, gid)
    if gid == "rule"
        r = ensure_rule!(p, p.loc_count)
        r.errordim_sys_on = -1
        r.errordim_sys_off = -1
    elseif gid == "channel"
        ensure_channel!(p, p.loc_count)
    elseif gid == "energy"
        ensure_smear!(p, p.loc_count)
        p.ibf = SmearData()
    elseif gid == "flux" || gid == "nuflux"
        p.flt = Flux()
    elseif gid == "cross"
        p.xsc = XSec()
    end
end

function group_end!(p::AEDLParser, gid)
    e = p.exp
    k = p.loc_count
    if gid == "energy"
        e.smear_data[k] = p.ibf
    elseif gid == "flux" || gid == "nuflux"
        if gid == "flux"
            glb_error(
                "The 'flux' directive is deprecated.\nThe flux normalization may not be what you expect.\nPlease, consult the manual!",
            )
            p.flt.builtin <= 0 && (p.flt.builtin = GLB_OLD_NORM)
        end
        while length(e.fluxes) < k
            push!(e.fluxes, Flux())
        end
        e.fluxes[k] = p.flt
    elseif gid == "cross"
        while length(e.xsecs) < k
            push!(e.xsecs, XSec())
        end
        e.xsecs[k] = p.xsc
    elseif gid == "rule"
        r = e.rules[k]
        r.sys_on_string === nothing &&
            r.errordim_sys_on >= 0 &&
            (r.sys_on_string = convert_errordim(r.errordim_sys_on))
        r.sys_off_string === nothing &&
            r.errordim_sys_off >= 0 &&
            (r.sys_off_string = convert_errordim(r.errordim_sys_off))
    end
end

function check_range(p, name, v, lo, hi)
    (lo <= v <= hi) || perr(p, "Value for $name out of range")
    return v
end

aslist(p, v::Float64) = [v]
aslist(p, v::Vector{Float64}) = copy(v)
aslist(p, v) = perr(p, "Expected a list")
asnum(p, v::Float64) = v
asnum(p, v) = perr(p, "Expected a number")

function set_counter!(p, name, field, v, lo, hi)
    check_range(p, name, v, lo, hi)
    old = getfield(p.exp, field)
    if !(old == -1 || old == Int(v))
        glb_warning("Given length does not match actual length")
        return
    end
    setfield!(p.exp, field, Int(v))
end

function set_list_with_length!(p, field, lenfield, v)
    l = aslist(p, v)
    old = getfield(p.exp, lenfield)
    if old == -1
        setfield!(p.exp, lenfield, length(l))
    elseif old != length(l)
        glb_warning("Length mismatch or list length changed")
    end
    setfield!(p.exp, field, l)
end

const GMAX = 1e100

function set_global!(p::AEDLParser, t::Token, v)
    e = p.exp
    name = t.text
    if name == "\$parent_energy"
        e.emax = check_range(p, name, asnum(p, v), 0, GMAX)
    elseif name == "\$target_mass"
        e.targetmass = check_range(p, name, asnum(p, v), 0, GMAX)
    elseif name == "\$sampling_points"
        set_counter!(p, name, :simbins, asnum(p, v), 0, 500)
    elseif name == "\$sampling_min"
        e.simtresh = check_range(p, name, asnum(p, v), 0, GMAX)
    elseif name == "\$sampling_max"
        e.simbeam = check_range(p, name, asnum(p, v), 0, GMAX)
    elseif name == "\$binsize"
        set_list_with_length!(p, :binsize, :numofbins, v)
    elseif name == "\$bins"
        set_counter!(p, name, :numofbins, asnum(p, v), 0, 500)
    elseif name == "\$emin"
        e.emin = check_range(p, name, asnum(p, v), 0, GMAX)
    elseif name == "\$emax"
        e.emax = check_range(p, name, asnum(p, v), 0, GMAX)
    elseif name == "\$baseline"
        e.baseline = check_range(p, name, asnum(p, v), 0, 2 * GLB_EARTH_RADIUS)
    elseif name == "\$profiletype"
        e.density_profile_type = Int(check_range(p, name, asnum(p, v), 1, 3))
    elseif name == "\$sampling_stepsize"
        set_list_with_length!(p, :simbinsize, :simbins, v)
    elseif name == "\$densitytab"
        set_list_with_length!(p, :densitytab, :psteps, v)
    elseif name == "\$lengthtab"
        l = aslist(p, v)
        foreach(x -> check_range(p, name, x, 0, 2 * GLB_EARTH_RADIUS), l)
        set_list_with_length!(p, :lengthtab, :psteps, l)
    elseif name == "\$densitysteps"
        set_counter!(p, name, :psteps, asnum(p, v), 1, GMAX)
    elseif name == "\$filter_state"
        e.filter_state = Int(check_range(p, name, asnum(p, v), GLB_OFF, GLB_ON))
    elseif name == "\$filter_value"
        e.filter_value = check_range(p, name, asnum(p, v), 0, GMAX)
    else
        perr(p, "Unknown identifier $name", t)
    end
end

function flavour_of(p::AEDLParser, t::Token)
    t.kind == :ident || perr(p, "Expected a flavour")
    s = t.text
    nosc = 0
    if startswith(s, "NOSC_")
        nosc = 10
        s = s[6:end]
    end
    occursin(r"^[eE][lectronLECTRON]*$", s) && return 1 + nosc
    occursin(r"^[mM][uonUON]*$", s) && return 2 + nosc
    occursin(r"^[tT][auAU]*$", s) && return 3 + nosc
    perr(p, "Unknown flavour $(t.text)", t)
end

function name_value(p::AEDLParser)
    t = advance!(p)
    t.kind == :hash || perr(p, "Expected #name", t)
    haskey(p.names, t.text) || perr(p, "Unknown name $(t.text)", t)
    return p.names[t.text][2]
end

function parse_rule_parts(p::AEDLParser)
    coeffs = Float64[]
    chans = Int[]
    while true
        c = scalar(p, parse_expr(p))
        peek(p).kind == :rulemult || perr(p, "Expected '@' in rule")
        advance!(p)
        ch = scalar(p, parse_expr(p))
        push!(coeffs, c)
        push!(chans, Int(ch))
        isop(peek(p), ":") || break
        advance!(p)
        skipnl!(p)
    end
    return coeffs, chans
end

function parse_energy_matrix!(p::AEDLParser)
    e = p.exp
    k = p.loc_count
    rows = Vector{Float64}[]
    while true
        skipnl!(p)
        push!(rows, aslist(p, parse_expr(p)))
        skipnl!(p)
        if isop(peek(p), ":")
            advance!(p)
        elseif isop(peek(p), ";")
            advance!(p)
            break
        else
            perr(p, "Expected ':' or ';' in @energy")
        end
    end
    mat = Vector{Vector{Float64}}(undef, length(rows))
    low = zeros(Int, length(rows))
    up = zeros(Int, length(rows))
    for (l, r) in enumerate(rows)
        length(r) < 2 && perr(p, "In @energy number $k: sublist $(l-1) is too short")
        v1 = Int(r[1])
        v2 = Int(r[2])
        ok = v1 >= 0 && v2 >= v1 && v2 - v1 == length(r) - 3 && (e.simbins == -1 || v2 <= e.simbins)
        ok || perr(p, "Value for ranges in smear out of range")
        low[l] = v1 + 1
        up[l] = v2 + 1
        mat[l] = r[3:end]
    end
    e.smear[k] = mat
    e.lowrange[k] = low
    e.uprange[k] = up
    e.user_smear_given[k] = true
end

function parse_local!(p::AEDLParser, t::Token)
    e = p.exp
    k = p.loc_count
    ctx = p.context
    name = t.text
    if ctx == "rule"
        r = e.rules[k]
        if name == "@signal"
            r.signal_coeff, r.signal_channels = parse_rule_parts(p)
        elseif name == "@background"
            r.bg_coeff, r.bg_channels = parse_rule_parts(p)
        elseif name == "@sys_on_function"
            r.sys_on_string = expect_string!(p)
        elseif name == "@sys_off_function"
            r.sys_off_string = expect_string!(p)
        elseif name in ("@signalerror", "@backgrounderror", "@backgroundcenter", "@energy_window")
            v = parse_rhs(p)
            v isa Tuple || perr(p, "$name expects a pair a : b")
            check_range(p, name, v[1], 0, name == "@energy_window" ? GMAX : 100)
            target =
                name == "@signalerror" ? r.signal_errors :
                name == "@backgrounderror" ? r.bg_errors :
                name == "@backgroundcenter" ? r.bg_centers : r.energy_window
            target[1] = v[1]
            target[2] = v[2]
        elseif name == "@sys_on_errors"
            r.sys_on_errors = aslist(p, parse_rhs(p))
        elseif name == "@sys_off_errors"
            r.sys_off_errors = aslist(p, parse_rhs(p))
        elseif name == "@errordim_sys_on"
            r.errordim_sys_on = Int(check_range(p, name, asnum(p, parse_rhs(p)), 0, 20))
        elseif name == "@errordim_sys_off"
            r.errordim_sys_off = Int(check_range(p, name, asnum(p, parse_rhs(p)), 0, 20))
        else
            perr(p, "Unknown identifier $name", t)
        end
    elseif ctx == "channel"
        if name == "@channel"
            fl = name_value(p)
            expect_op!(p, ":")
            pm = advance!(p)
            cp = isop(pm, "+") ? 1 : isop(pm, "-") ? -1 : perr(p, "Expected + or -", pm)
            expect_op!(p, ":")
            ini = flavour_of(p, advance!(p))
            expect_op!(p, ":")
            fin = flavour_of(p, advance!(p))
            expect_op!(p, ":")
            xs = name_value(p)
            expect_op!(p, ":")
            en = name_value(p)
            e.channels[k] = Channel(fl, cp, ini, fin, xs, en)
        elseif name == "@pre_smearing_efficiencies"
            e.user_pre_smearing_channel[k] = aslist(p, parse_rhs(p))
        elseif name == "@pre_smearing_background"
            e.user_pre_smearing_background[k] = aslist(p, parse_rhs(p))
        elseif name == "@post_smearing_efficiencies"
            e.user_post_smearing_channel[k] = aslist(p, parse_rhs(p))
        elseif name == "@post_smearing_background"
            e.user_post_smearing_background[k] = aslist(p, parse_rhs(p))
        else
            perr(p, "Unknown identifier $name", t)
        end
    elseif ctx == "energy"
        sd = p.ibf
        if name == "@energy"
            parse_energy_matrix!(p)
        elseif name == "@type"
            sd.type = Int(check_range(p, name, asnum(p, parse_rhs(p)), 1, 2))
        elseif name == "@sigma_e"
            sd.sigma = aslist(p, parse_rhs(p))
            sd.num_of_params = length(sd.sigma)
        elseif name == "@sigma_function"
            nt = advance!(p)
            (nt.kind == :hash && haskey(SIGMA_FUNCTIONS, nt.text)) ||
                perr(p, "Unknown sigma function $(nt.text)", nt)
            sd.sig_f = SIGMA_FUNCTIONS[nt.text]
        else
            perr(p, "Unknown identifier $name", t)
        end
    elseif ctx == "cross"
        if name == "@cross_file"
            p.xsc.file_name = expect_data_file!(p)
        else
            perr(p, "Unknown identifier $name", t)
        end
    elseif ctx == "flux" || ctx == "nuflux"
        f = p.flt
        if name == "@flux_file" || name == "@flux_file_m"
            f.file_name = expect_data_file!(p)
        elseif name == "@builtin"
            f.builtin = Int(check_range(p, name, asnum(p, parse_rhs(p)), 1, 4))
        else
            v = check_range(p, name, asnum(p, parse_rhs(p)), 0, GMAX)
            if name == "@time"
                f.time = v
            elseif name == "@power"
                f.target_power = v
            elseif name == "@stored_muons" || name == "@stored_ions"
                f.stored_muons = v
            elseif name == "@parent_energy"
                f.parent_energy = v
            elseif name == "@norm"
                f.norm = v
            elseif name == "@gamma"
                f.gamma = v
            elseif name == "@end_point"
                f.end_point = v
            else
                perr(p, "Unknown identifier $name", t)
            end
        end
    end
end

"""
    parse_aedl_file(filename) -> Vector{ExperimentData}

Parse an AEDL file into (not yet initialized) experiments. Several
experiments can be separated by `#NEXT#`.
"""
function parse_aedl_file(path::AbstractString)
    open(path) do io
        hdr = String(read(io, min(8, filesize(path))))
        hdr == "%!GLoBES" || throw(AEDLError("Not a GLoBES file: $path"))
    end
    toks = aedl_tokenize(path; display_name = path)
    p = new_parser(toks, path)
    return parse_aedl!(p)
end

"""
    glbDefineAEDLVariable(name, value)

Predefine a variable that AEDL files read later can use, e.g. `BASELINE`.
"""
glbDefineAEDLVariable(name::AbstractString, value::Real) =
    (aedl_presym_vars[String(name)] = Float64(value); nothing)

"""
    glbDefineAEDLList(name, list)

Predefine an AEDL list; `name` must start with `%`.
"""
function glbDefineAEDLList(name::AbstractString, list::AbstractVector{<:Real})
    startswith(name, "%") || (glb_error("AEDL lists have to start with '%'"); return)
    aedl_presym_lists[String(name)] = Float64.(collect(list))
    nothing
end

"""
    glbClearAEDLVariables()

Remove all variables and lists defined with `glbDefineAEDLVariable` and
`glbDefineAEDLList`.
"""
function glbClearAEDLVariables()
    empty!(aedl_presym_vars)
    empty!(aedl_presym_lists)
    nothing
end
