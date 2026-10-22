# A reader for the subset of the Wolfram language used in SARAH model files:
# assignments, lists, rules, function calls f[x], parts x[[i]], strings,
# numbers, + - * / ^, the non-commutative product `.` and comments (* *).
#
# Expressions become `MExpr(head, args)`; symbols are Julia `Symbol`s, numbers
# `Int`, `Rational{Int}` or `Float64`, strings `String`. `a - b` is read as
# Plus[a, Times[-1, b]] and `a/b` as Times[a, Power[b, -1]], as in Mathematica.

"""
    MExpr(head, args)

A Wolfram-language expression `head[args...]` read from a SARAH model file.
"""
struct MExpr
    head::Any
    args::Vector{Any}
    MExpr(head, args::Vector{Any}) = new(head, args)
end
MExpr(head, args...) = MExpr(head, collect(Any, args))

Base.:(==)(a::MExpr, b::MExpr) = a.head == b.head && a.args == b.args
Base.hash(e::MExpr, h::UInt) = hash(e.args, hash(e.head, h))

function Base.show(io::IO, e::MExpr)
    print(io, e.head === :List ? "" : e.head, e.head === :List ? "{" : "[")
    join(io, (sprint(show_mexpr, a) for a in e.args), ", ")
    print(io, e.head === :List ? "}" : "]")
end
show_mexpr(io, x::Symbol) = print(io, x)
show_mexpr(io, x::Rational) = print(io, numerator(x), "/", denominator(x))
show_mexpr(io, x) = show(io, x)

is_call(x, head) = x isa MExpr && x.head === head

struct MathematicaSyntaxError <: Exception
    msg::String
    line::Int
end
Base.showerror(io::IO, e::MathematicaSyntaxError) =
    print(io, "syntax error in model file, line ", e.line, ": ", e.msg)

# ---------------------------------------------------------------------------
#  Tokens
# ---------------------------------------------------------------------------

struct MToken
    kind::Symbol          # :sym, :num, :str, :op, :eof
    value::Any
    line::Int
    newline::Bool         # a line break separates it from the previous token
    space::Bool           # whitespace separates it from the previous token
end

# Named characters that SARAH files use as operators or letters
const NAMED_CHARACTERS = Dict(
    "RightArrow" => "->",
    "Rule" => "->",
    "RuleDelayed" => ":>",
    "Times" => "*",
    "Alpha" => "α",
    "Beta" => "β",
    "Gamma" => "γ",
    "Delta" => "δ",
    "Epsilon" => "ϵ",
    "Eta" => "η",
    "Theta" => "θ",
    "Lambda" => "λ",
    "Mu" => "μ",
    "Nu" => "ν",
    "Pi" => "π",
    "Sigma" => "σ",
    "Tau" => "τ",
    "Phi" => "ϕ",
    "Chi" => "χ",
    "Psi" => "ψ",
    "Omega" => "ω",
    "Congruent" => "≡",
)

const OPERATORS = ("[[", ":=", ":>", "->", "==", "!=", "&&", "||", "@@", "/.", "//")

function tokenize(src::AbstractString)
    s = collect(src)
    n = length(s)
    tokens = MToken[]
    i, line, newline, space = 1, 1, false, false
    err(msg) = throw(MathematicaSyntaxError(msg, line))
    push(kind, value) =
        (push!(tokens, MToken(kind, value, line, newline, space)); newline =  space = false )
    peek(k = 0) = i + k <= n ? s[i + k] : '\0'
    while i <= n
        c = s[i]
        if c == '\n'
            line += 1
            newline = space = true
            i += 1
        elseif isspace(c)
            space = true
            i += 1
        elseif c == '(' && peek(1) == '*'        # comment, possibly nested
            space = true
            depth, i = 1, i + 2
            while depth > 0
                i > n && err("unterminated comment")
                if s[i] == '(' && peek(1) == '*'
                    depth += 1
                    i += 2
                elseif s[i] == '*' && peek(1) == ')'
                    depth -= 1
                    i += 2
                else
                    s[i] == '\n' && (line += 1; newline = true)
                    i += 1
                end
            end
        elseif c == '"'
            buf = IOBuffer()
            i += 1
            while true
                i > n && err("unterminated string")
                if s[i] == '\\' && i < n
                    nxt = s[i + 1]
                    write(buf, nxt == 'n' ? '\n' : nxt == 't' ? '\t' : nxt)
                    i += 2
                elseif s[i] == '"'
                    i += 1
                    break
                else
                    s[i] == '\n' && (line += 1)
                    write(buf, s[i])
                    i += 1
                end
            end
            push(:str, String(take!(buf)))
        elseif isdigit(c) || (c == '.' && isdigit(peek(1)) && !last_is_operand(tokens))
            j = i
            while j <= n && isdigit(s[j])
                j += 1
            end
            isfloat = false
            if j <= n && s[j] == '.' && (j + 1 > n || isdigit(s[j + 1]) || !isletter(s[j + 1]))
                # a decimal point, unless it is a Dot product like 2.x
                if j + 1 <= n && (isdigit(s[j + 1]) || !(isletter(s[j + 1]) || s[j + 1] in "([{"))
                    isfloat = true
                    j += 1
                    while j <= n && isdigit(s[j])
                        j += 1
                    end
                end
            end
            mantissa = String(s[i:(j - 1)])
            exponent = 0
            if j + 1 <= n && s[j] == '*' && s[j + 1] == '^'          # 1.5*^-3
                k = j + 2
                k <= n && s[k] in "+-" && (k += 1)
                while k <= n && isdigit(s[k])
                    k += 1
                end
                exponent = parse(Int, String(s[(j + 2):(k - 1)]))
                isfloat = true
                j = k
            elseif j <= n && s[j] in "eE" && j + 1 <= n && (isdigit(s[j + 1]) || s[j + 1] in "+-")
                k = j + 1
                s[k] in "+-" && (k += 1)
                while k <= n && isdigit(s[k])
                    k += 1
                end
                exponent = parse(Int, String(s[(j + 1):(k - 1)]))
                isfloat = true
                j = k
            end
            push(:num, isfloat ? parse(Float64, mantissa) * 10.0^exponent : parse(Int, mantissa))
            i = j
        elseif isletter(c) || c == '$' || c == '\\' && peek(1) == '['
            buf = IOBuffer()
            while i <= n
                if s[i] == '\\' && i < n && s[i + 1] == '['
                    k = findnext(==(']'), s, i)
                    k === nothing && err("unterminated named character")
                    name = String(s[(i + 2):(k - 1)])
                    rep = get(NAMED_CHARACTERS, name, name)
                    if rep in ("->", ":>", "*")
                        position(buf) == 0 || break          # an operator ends the symbol
                        push(:op, rep)
                        i = k + 1
                        @goto next_token
                    end
                    write(buf, rep)
                    i = k + 1
                elseif isletter(s[i]) || isdigit(s[i]) || s[i] in "\$`"
                    write(buf, s[i])
                    i += 1
                elseif s[i] == '_' && i < n && (isletter(s[i + 1]) || isdigit(s[i + 1]))
                    write(buf, '_')                       # eta_etau (not a Wolfram pattern)
                    i += 1
                elseif s[i] == ':' && i < n && s[i + 1] == ':'      # General::spell
                    write(buf, "::")
                    i += 2
                else
                    break
                end
            end
            push(:sym, Symbol(String(take!(buf))))
            @label next_token
        else
            op = nothing
            for o in OPERATORS
                if i + length(o) - 1 <= n && String(s[i:(i + length(o) - 1)]) == o
                    op = o
                    break
                end
            end
            op === nothing &&
                (c in "=+-*/^.,;{}[]()&@!_?|<>~'#%" || err("unexpected character '$c'"))
            op = something(op, string(c))
            push(:op, op)
            i += length(op)
        end
    end
    push!(tokens, MToken(:eof, nothing, line, true, true))
    return tokens
end

last_is_operand(tokens) =
    !isempty(tokens) &&
    (tokens[end].kind in (:sym, :num, :str) || tokens[end].value in ("]", ")", "}"))

# ---------------------------------------------------------------------------
#  Parser (precedence climbing)
# ---------------------------------------------------------------------------

mutable struct MParser
    tokens::Vector{MToken}
    pos::Int
    depth::Int            # open brackets: inside them, line breaks are spaces
end

peek_token(p::MParser) = p.tokens[p.pos]
next_token!(p::MParser) = (t = p.tokens[p.pos]; p.pos += 1; t)
is_op(t::MToken, op) = t.kind === :op && t.value == op
syntax_error(p::MParser, msg) = throw(MathematicaSyntaxError(msg, peek_token(p).line))

function expect!(p::MParser, op)
    t = next_token!(p)
    is_op(t, op) || throw(
        MathematicaSyntaxError(
            "expected '$op', found $(t.kind === :eof ? "end of file" : repr(t.value))",
            t.line,
        ),
    )
    return t
end

"""
    parse_mathematica(src) -> Vector{Any}

Parse the statements (separated by `;` or line breaks) of Wolfram-language
source text into `MExpr` trees.
"""
function parse_mathematica(src::AbstractString)
    p = MParser(tokenize(src), 1, 0)
    statements = Any[]
    while peek_token(p).kind !== :eof
        if is_op(peek_token(p), ";")
            next_token!(p)
            continue
        end
        push!(statements, parse_expr(p, 0))
        t = peek_token(p)
        t.kind === :eof ||
            is_op(t, ";") ||
            t.newline ||
            syntax_error(p, "unexpected $(repr(t.value))")
    end
    return statements
end

# Binding powers of infix operators: (left, right)
const INFIX = Dict(
    "=" => (10, 9),
    ":=" => (10, 9),
    "->" => (20, 19),
    ":>" => (20, 19),
    "/." => (15, 16),
    "||" => (25, 26),
    "&&" => (27, 28),
    "==" => (30, 31),
    "!=" => (30, 31),
    "+" => (40, 41),
    "-" => (40, 41),
    "*" => (50, 51),
    "/" => (50, 51),
    "." => (60, 61),
    "^" => (71, 70),
)
const IMPLICIT_TIMES = (50, 51)
const PREFIX_MINUS = 55
const POSTFIX = 80

starts_operand(t::MToken) =
    t.kind in (:sym, :num, :str) || (t.kind === :op && t.value in ("(", "{"))

function parse_expr(p::MParser, min_bp)
    lhs = parse_prefix(p)
    while true
        t = peek_token(p)
        # A line break ends the statement at the top level
        at_break = t.newline && p.depth == 0
        if t.kind === :op && (t.value == "[" || t.value == "[[") && !at_break
            POSTFIX < min_bp && break
            lhs = parse_postfix(p, lhs)
            continue
        end
        # Head[...]{ ... } with no space before the brace: a block of definitions
        if is_op(t, "{") && !t.space && (lhs isa Symbol || lhs isa MExpr && lhs.head isa Symbol)
            POSTFIX < min_bp && break
            next_token!(p)
            lhs = MExpr(:Block, Any[lhs; parse_block_items(p)])
            continue
        end
        if t.kind === :op && haskey(INFIX, t.value)
            lbp, rbp = INFIX[t.value]
            (lbp < min_bp || at_break) && break
            next_token!(p)
            rhs = parse_expr(p, rbp)
            lhs = combine(t.value, lhs, rhs)
        elseif starts_operand(t) && !at_break && !ends_item_statement(lhs)
            lbp, rbp = IMPLICIT_TIMES
            lbp < min_bp && break
            rhs = parse_expr(p, rbp)
            lhs = times(lhs, rhs)
        else
            break
        end
    end
    return lhs
end

function parse_prefix(p::MParser)
    t = next_token!(p)
    if t.kind in (:num, :str)
        return t.value
    elseif t.kind === :sym
        return t.value
    elseif is_op(t, "-")
        return times(-1, parse_expr(p, PREFIX_MINUS))
    elseif is_op(t, "+")
        return parse_expr(p, PREFIX_MINUS)
    elseif is_op(t, "(")
        p.depth += 1
        e = parse_expr(p, 0)
        while is_op(peek_token(p), ";")             # (a; b) evaluates to b
            next_token!(p)
            e = parse_expr(p, 0)
        end
        expect!(p, ")")
        p.depth -= 1
        return e
    elseif is_op(t, "{")
        p.depth += 1
        items = parse_sequence(p, "}")
        p.depth -= 1
        return MExpr(:List, items)
    else
        throw(
            MathematicaSyntaxError(
                "unexpected $(t.kind === :eof ? "end of file" : repr(t.value))",
                t.line,
            ),
        )
    end
end

# Items of a block, separated by commas, semicolons or line breaks
function parse_block_items(p::MParser)
    saved = p.depth
    p.depth = 0
    items = Any[]
    while true
        t = peek_token(p)
        if is_op(t, ",") || is_op(t, ";")
            next_token!(p)
        elseif is_op(t, "}")
            next_token!(p)
            break
        elseif t.kind === :eof
            syntax_error(p, "unterminated block")
        else
            push!(items, parse_expr(p, 0))
            u = peek_token(p)
            u.newline ||
                is_op(u, ",") ||
                is_op(u, ";") ||
                is_op(u, "}") ||
                ends_item_statement(items[end]) ||
                syntax_error(p, "unexpected $(repr(u.value)) in block")
        end
    end
    p.depth = saved
    return items
end

# Switches and print commands do not multiply what follows them, so that
# a block can list `X = Off Print["..."] Y = On` on one line
ends_item(e) =
    e === :On || e === :Off || is_call(e, :Print) || is_call(e, :printf) || is_call(e, :Section)
ends_item_statement(e) = ends_item(e) || is_call(e, :Set) && ends_item(e.args[2])

function parse_sequence(p::MParser, close)
    items = Any[]
    if is_op(peek_token(p), close)
        next_token!(p)
        return items
    end
    while true
        push!(items, parse_expr(p, 0))
        t = next_token!(p)
        is_op(t, close) && return items
        is_op(t, ",") || throw(MathematicaSyntaxError("expected ',' or '$close'", t.line))
    end
end

function parse_postfix(p::MParser, head)
    t = next_token!(p)
    p.depth += 1
    if t.value == "[["
        items = parse_sequence(p, "]")
        expect!(p, "]")
        p.depth -= 1
        return MExpr(:Part, Any[head; items])
    end
    items = parse_sequence(p, "]")
    p.depth -= 1
    return MExpr(head, items)
end

function combine(op, a, b)
    op == "+" && return plus(a, b)
    op == "-" && return plus(a, times(-1, b))
    op == "*" && return times(a, b)
    op == "/" && return times(a, power(b, -1))
    op == "^" && return power(a, b)
    op == "." && return MExpr(:Dot, Any[flatten_args(:Dot, a); flatten_args(:Dot, b)])
    head = Dict(
        "=" => :Set,
        ":=" => :SetDelayed,
        "->" => :Rule,
        ":>" => :RuleDelayed,
        "/." => :ReplaceAll,
        "==" => :Equal,
        "!=" => :Unequal,
        "&&" => :And,
        "||" => :Or,
    )[op]
    return MExpr(head, a, b)
end

flatten_args(head, x) = is_call(x, head) ? x.args : Any[x]

const MNumber = Union{Int, Rational{Int}, Float64}

function times(a, b)
    args = Any[flatten_args(:Times, a); flatten_args(:Times, b)]
    c = 1
    rest = Any[]
    for x in args
        x isa MNumber ? (c *= x) : push!(rest, x)
    end
    c isa Rational && denominator(c) == 1 && (c = numerator(c))
    isempty(rest) && return c
    c == 0 && return 0
    c == 1 && length(rest) == 1 && return rest[1]
    c == 1 || pushfirst!(rest, c)
    return MExpr(:Times, rest)
end

function plus(a, b)
    a isa MNumber && b isa MNumber && return a + b
    return MExpr(:Plus, Any[flatten_args(:Plus, a); flatten_args(:Plus, b)])
end

function power(a, b)
    if a isa Integer && b isa Integer
        return b >= 0 ? a^b : 1 // a^(-b)
    end
    a isa MNumber && b isa MNumber && return float(a)^b
    return MExpr(:Power, a, b)
end
