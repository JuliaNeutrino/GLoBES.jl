# Tree-level spectrum of a BSMModel after electroweak symmetry breaking, and
# parameter values in LesHouches files.

# ---------------------------------------------------------------------------
#  Numerical values of SARAH expressions
# ---------------------------------------------------------------------------

"""
    evaluate(expr, values) -> Number or Array

Numerical value of a SARAH expression, with the symbols taken from `values`
(a `Dict` or named tuple). `I`, `Pi`, `Sqrt`, `Exp`, `Log`, `Sin`, `Cos`,
`Conjugate` and `conj` are known.
"""
function evaluate(e, values)
    e isa MNumber && return e
    if e isa Symbol
        e === :I && return im
        e === :Pi && return π
        e === :True && return true
        e === :False && return false
        v = lookup(values, e)
        v === nothing && throw(ArgumentError("no value for parameter $e"))
        return v
    end
    e isa MExpr || return e
    args = [evaluate(a, values) for a in e.args]
    h = e.head
    h === :Plus && return sum(args)
    h === :Times && return prod(args)
    h === :Dot && return prod(args)
    h === :Power && return args[1]^args[2]
    h === :List && return args
    h === :Sqrt && return sqrt(complex_if_negative(args[1]))
    h === :Exp && return exp(args[1])
    h === :Log && return log(args[1])
    h === :Sin && return sin(args[1])
    h === :Cos && return cos(args[1])
    (h === :conj || h === :Conjugate) && return conj(args[1])
    h === :Abs && return abs(args[1])
    throw(ArgumentError("cannot evaluate $h[...]"))
end

complex_if_negative(x::Real) = x < 0 ? complex(x) : x
complex_if_negative(x) = x

lookup(values::AbstractDict, k::Symbol) = get(values, k, get(values, String(k), nothing))
lookup(values::NamedTuple, k::Symbol) = get(values, k, nothing)

# ---------------------------------------------------------------------------
#  Scalar potential
# ---------------------------------------------------------------------------

# How each scalar component depends on the real fields after EWSB:
# component = vev + Σ coefficient × real field
struct ScalarComponent
    name::Symbol
    vev::Any                                 # expression of the VEV term
    terms::Vector{Tuple{Symbol, Any}}        # (real field, coefficient expression)
end

function scalar_components(m::BSMModel)
    vevs = Dict{Symbol, Any}()
    if haskey(m.definitions, (:EWSB, :VEVs))
        for entry in m.definitions[(:EWSB, :VEVs)].args
            vevs[entry.args[1]] = entry.args[2:end]
        end
    end
    comps = ScalarComponent[]
    for f in scalars(m)
        f.generations == 1 || throw(
            ArgumentError(
                "scalar $(f.name) has $(f.generations) generations; only one is supported",
            ),
        )
        for c in f.components
            if haskey(vevs, c)
                (v, cv), rest = vevs[c][1].args, vevs[c][2:end]
                vev = v == 0 ? 0 : MExpr(:Times, Any[v, cv])
                push!(comps, ScalarComponent(c, vev, [(r.args[1], r.args[2]) for r in rest]))
            else
                s = MExpr(:Power, MExpr(:Sqrt, 2), -1)
                push!(
                    comps,
                    ScalarComponent(
                        c,
                        0,
                        [(Symbol(c, "_R"), s), (Symbol(c, "_I"), MExpr(:Times, Any[:I, s]))],
                    ),
                )
            end
        end
    end
    return comps
end

# The scalar potential V(x) = -L(x) as a function of the real fields x
struct ScalarPotential
    model::BSMModel
    components::Vector{ScalarComponent}
    real_fields::Vector{Symbol}
    monomials::Vector{Tuple{Monomial, Bool}}       # (monomial, add h.c.)
end

function ScalarPotential(m::BSMModel)
    comps = scalar_components(m)
    names = unique!([t[1] for c in comps for t in c.terms])
    monos = Tuple{Monomial, Bool}[]
    for t in m.lagrangian, mono in expand_lagrangian(m, t.expr)
        isempty(mono.fields) && continue
        all(((f, c),) -> field_lookup(m, f).kind === :scalar, mono.fields) || continue
        push!(monos, (mono, t.hc))
    end
    ScalarPotential(m, comps, names, monos)
end

function (V::ScalarPotential)(x::AbstractVector, values)
    m = V.model
    fieldval = Dict(zip(V.real_fields, x))
    compval = Dict{Symbol, ComplexF64}()
    for c in V.components
        z = complex(float(evaluate(c.vev, values)))
        for (r, coeff) in c.terms
            z += evaluate(coeff, values) * fieldval[r]
        end
        compval[c.name] = z
    end
    multiplet(f) = [compval[c] for c in field_lookup(m, f).components]
    L = zero(ComplexF64)
    for (mono, hc) in V.monomials
        val = complex(float(evaluate(mono.coefficient, values)))
        for p in mono.parameters
            val *= evaluate(p, values)
        end
        fs = mono.fields
        isodd(length(fs)) &&
            throw(ArgumentError("odd number of scalars in $(monomial_string(mono))"))
        for k = 1:2:length(fs)
            (a, ca), (b, cb) = fs[k], fs[k + 1]
            va, vb = multiplet(a), multiplet(b)
            ca && (va = conj.(va))
            cb && (vb = conj.(vb))
            val *= if length(va) == 1 || length(vb) == 1
                va[1] * vb[1]
            elseif ca != cb
                sum(va .* vb)                    # H†η
            else
                va[1] * vb[2] - va[2] * vb[1]    # ε contraction of two doublets
            end
        end
        L += hc ? val + conj(val) : val
    end
    return -real(L)
end

# Second derivatives of a polynomial of degree ≤ 4: Richardson extrapolation of
# central differences is exact up to rounding
function hessian(f, x0::Vector{Float64}, h)
    n = length(x0)
    H = zeros(n, n)
    function d2(i, j, h)
        e(k) = (v = zeros(n); v[k] = h; v)
        i == j && return (f(x0 + e(i)) - 2f(x0) + f(x0 - e(i))) / h^2
        return (
            f(x0 + e(i) + e(j)) - f(x0 + e(i) - e(j)) - f(x0 - e(i) + e(j)) + f(x0 - e(i) - e(j))
        ) / (4h^2)
    end
    for i = 1:n, j = i:n
        H[i, j] = H[j, i] = (4d2(i, j, h / 2) - d2(i, j, h)) / 3
    end
    return H
end

function gradient(f, x0::Vector{Float64}, h)
    n = length(x0)
    g = zeros(n)
    for i = 1:n
        e = zeros(n)
        e[i] = h
        g1 = (f(x0 + e) - f(x0 - e)) / (2h)
        g2 = (f(x0 + e / 2) - f(x0 - e / 2)) / h
        g[i] = (4g2 - g1) / 3
    end
    return g
end

# ---------------------------------------------------------------------------
#  Spectrum
# ---------------------------------------------------------------------------

"""
    Spectrum

Tree-level masses after electroweak symmetry breaking, from
[`spectrum`](@ref): `scalars` maps each mass eigenstate to its mass in GeV
(negative for a tachyon, `-sqrt(-m²)`), `fermions` the Majorana masses of
gauge-singlet fermions, `tadpoles` the parameters fixed by the minimum
conditions and `values` all parameter values used.
"""
struct Spectrum
    model::String
    scalars::Vector{Pair{Symbol, Float64}}
    fermions::Vector{Pair{Symbol, Vector{Float64}}}
    tadpoles::Dict{Symbol, Float64}
    values::Dict{Symbol, Any}
end

function Base.show(io::IO, ::MIME"text/plain", s::Spectrum)
    println(io, "Tree-level spectrum of \"", s.model, "\" (GeV)")
    for (k, v) in s.tadpoles
        @printf(io, "  %-10s = %.6g  (tadpole)\n", k, v)
    end
    for (n, mass) in s.scalars
        @printf(io, "  m(%s) = %.6g%s\n", n, mass, mass < 0 ? "  (tachyonic: m² < 0)" : "")
    end
    for (n, ms) in s.fermions
        println(io, "  m($n) = ", join((@sprintf("%.6g", x) for x in ms), ", "))
    end
end

"""
    spectrum(model, values; tadpoles=from SPheno.m) -> Spectrum

Tree-level masses of the scalars after electroweak symmetry breaking (from
the potential in the Lagrangian and `DEFINITION[EWSB][VEVs]`) and of the
gauge-singlet Majorana fermions. `values` holds the numerical parameters,
by their SARAH names. The parameters in `tadpoles` (by default SARAH's
`ParametersToSolveTadpoles`) are solved from the minimum conditions.

```julia
values = Dict(:v => 246.22, :lambda1 => 0.26, :lambda2 => 0.1, :lambda3 => 0.2, :lambda4 => -0.1, :lambda5 => 1e-9, :mEt2 => 1e6, :Mn => Diagonal([1e3, 2e3, 3e3]))   # parameters by SARAH name, in GeV
s = spectrum(scoto, values; tadpoles = [:mH2])   # model, values, parameters fixed by the minimum conditions
s.scalars, s.fermions, s.tadpoles   # masses in GeV and the solved parameters
```
"""
function spectrum(m::BSMModel, values; tadpoles = nothing)
    vals = Dict{Symbol, Any}(Symbol(k) => v for (k, v) in pairs(values))
    haskey(vals, :v) || (vals[:v] = HIGGS_VEV)
    V = ScalarPotential(m)
    tp =
        tadpoles === nothing ?
        (
            haskey(m.spheno, :ParametersToSolveTadpoles) ?
            Symbol[s for s in m.spheno[:ParametersToSolveTadpoles].args] : Symbol[]
        ) : collect(Symbol, tadpoles)
    x0 = zeros(length(V.real_fields))
    h = max(1.0, 1e-3 * abs(float(get(vals, :v, 1.0))))
    solved = Dict{Symbol, Float64}()
    if !isempty(tp)
        # The minimum conditions are linear in the tadpole parameters
        f(p) = (for (k, x) in zip(tp, p)
            vals[k] = x
        end; gradient(x -> V(x, vals), x0, h))
        g0 = f(zeros(length(tp)))
        A = hcat((f([i == j ? 1.0 : 0.0 for j = 1:length(tp)]) - g0 for i = 1:length(tp))...)
        p = -(A \ g0)
        f(p)
        for (k, x) in zip(tp, p)
            solved[k] = x
        end
    end
    M2 = hessian(x -> V(x, vals), x0, h)
    λ, Q = eigen(Symmetric(M2))
    λ[abs.(λ) .< 1e-9 * maximum(abs, λ; init = 1.0)] .= 0      # Goldstone bosons
    scalars_out = Pair{Symbol, Float64}[]
    taken = Set{Symbol}()
    for i in sortperm(λ)
        k = argmax(abs.(Q[:, i]))
        name = V.real_fields[k]
        base = Symbol(replace(String(name), r"_[RI]$" => ""))
        # The two real fields of a charged scalar are one particle
        endswith(String(name), "_R") || endswith(String(name), "_I") ? (base in taken && continue) :
        nothing
        push!(taken, base)
        push!(scalars_out, base => (λ[i] >= 0 ? sqrt(λ[i]) : -sqrt(-λ[i])))
    end
    fermions_out = Pair{Symbol, Vector{Float64}}[]
    for t in m.lagrangian, mono in expand_lagrangian(m, t.expr)
        length(mono.fields) == 2 || continue
        (a, _), (b, _) = mono.fields
        fa = field_lookup(m, a)
        fa.kind === :fermion && a === b || continue
        any(f -> f.name === a, fermions(m)) || continue
        M = evaluate(MExpr(:Times, Any[2 * mono.coefficient; mono.parameters]), vals)
        push!(fermions_out, a => sort(svdvals(M isa Number ? fill(M, 1, 1) : Matrix(M))))
    end
    return Spectrum(m.name, scalars_out, fermions_out, solved, vals)
end

mass(s::Spectrum, name::Symbol) = last(s.scalars[findfirst(p -> p.first === name, s.scalars)])

# ---------------------------------------------------------------------------
#  LesHouches files
# ---------------------------------------------------------------------------

"""
    read_leshouches(file) -> Dict{String, Dict}

The blocks of a LesHouches (SLHA) file: block name ⇒ (indices ⇒ value), with
`IM<name>` blocks merged into the imaginary parts.

```julia
blocks = read_leshouches("LesHouches.in")   # file: blocks["MINPAR"], blocks["YN"], ...
```
"""
function read_leshouches(file::AbstractString)
    blocks = Dict{String, Dict{Tuple, Any}}()
    current = nothing
    for raw in eachline(file)
        line = strip(split(raw, '#'; limit = 2)[1])
        isempty(line) && continue
        words = split(line)
        if uppercase(words[1]) == "BLOCK"
            current = uppercase(words[2])
            blocks[current] = Dict{Tuple, Any}()
        elseif uppercase(words[1]) == "DECAY"
            current = nothing
        elseif current !== nothing
            nums = [tryparse(Float64, w) for w in words]
            any(isnothing, nums) && continue
            idx = Tuple(Int.(nums[1:(end - 1)]))
            blocks[current][idx] = nums[end]
        end
    end
    for name in collect(keys(blocks))
        startswith(name, "IM") && haskey(blocks, name[3:end]) || continue
        re = blocks[name[3:end]]
        for (k, v) in blocks[name]
            re[k] = get(re, k, 0.0) + im * v
        end
        delete!(blocks, name)
    end
    return blocks
end

"""
    leshouches_values(model, file) -> Dict{Symbol, Any}

Parameter values from a LesHouches file, at the positions that
`parameters.m` gives (`LesHouches -> {HDM, 1}` or `LesHouches -> YN`), and
from `MINPAR` through `BoundaryLowScaleInput` of `SPheno.m`.

```julia
values = leshouches_values(scoto, "LesHouches.in")   # model, LesHouches file: parameter values by SARAH name
```
"""
function leshouches_values(m::BSMModel, file::AbstractString)
    blocks = read_leshouches(file)
    vals = Dict{Symbol, Any}()
    for (p, info) in m.parameter_info
        haskey(info, :LesHouches) || continue
        pos = info[:LesHouches]
        if pos isa Symbol
            b = get(blocks, uppercase(String(pos)), nothing)
            b === nothing && continue
            n = maximum(k -> maximum(k), keys(b))
            M = zeros(ComplexF64, n, n)
            for (k, v) in b
                length(k) == 2 && (M[k...] = v)
            end
            vals[p] = all(iszero ∘ imag, M) ? real(M) : M
        elseif is_call(pos, :List)
            b = get(blocks, uppercase(String(pos.args[1])), nothing)
            b === nothing || haskey(b, (pos.args[2],)) && (vals[p] = b[(pos.args[2],)])
        end
    end
    if haskey(m.spheno, :MINPAR) && haskey(blocks, "MINPAR")
        inputs = Dict(
            e.args[2] => blocks["MINPAR"][(e.args[1],)] for
            e in m.spheno[:MINPAR].args if haskey(blocks["MINPAR"], (e.args[1],))
        )
        if haskey(m.spheno, :BoundaryLowScaleInput)
            for e in m.spheno[:BoundaryLowScaleInput].args
                haskey(inputs, e.args[2]) && (vals[e.args[1]] = inputs[e.args[2]])
            end
        end
    end
    return vals
end

"""
    write_leshouches(file, model, values)

Write parameter values to a LesHouches input file in the blocks of
`parameters.m` and `SPheno.m` (`MINPAR`), for SPheno or for
[`leshouches_values`](@ref).

```julia
write_leshouches("LesHouches.in", scoto, values)   # file, model, parameter values by SARAH name
```
"""
function write_leshouches(file::AbstractString, m::BSMModel, values)
    vals = Dict{Symbol, Any}(Symbol(k) => v for (k, v) in pairs(values))
    open(file, "w") do io
        println(io, "# LesHouches input for model ", m.name, ", written by GLoBES.jl")
        if haskey(m.spheno, :MINPAR) && haskey(m.spheno, :BoundaryLowScaleInput)
            to_param = Dict(e.args[2] => e.args[1] for e in m.spheno[:BoundaryLowScaleInput].args)
            println(io, "Block MINPAR")
            for e in m.spheno[:MINPAR].args
                p = get(to_param, e.args[2], nothing)
                p !== nothing &&
                    haskey(vals, p) &&
                    @printf(io, "%4d  %16.8E   # %s\n", e.args[1], real(vals[p]), e.args[2])
            end
        end
        for (p, info) in sort!(collect(m.parameter_info); by = first)
            pos = get(info, :LesHouches, nothing)
            pos isa Symbol && haskey(vals, p) || continue
            M = vals[p]
            for part in (real, imag)
                part === imag && all(iszero ∘ imag, M) && continue
                println(io, "Block ", part === imag ? "IM" : "", uppercase(String(pos)), "   # ", p)
                for i in axes(M, 1), j in axes(M, 2)
                    @printf(io, "%3d %3d  %16.8E\n", i, j, part(M[i, j]))
                end
            end
        end
    end
    return file
end
