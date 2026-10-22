# Oscillation physics beyond the standard three-flavour picture, assembled
# from blocks: scalar NSI, matter (vector) NSI and decoherence. A `NewPhysics`
# model turns its blocks into a probability engine with one named parameter
# per coupling, so the rest of GLoBES (rates, χ², fits) works unchanged.
#
# The Hamiltonian in flavour basis, for neutrinos of energy E, is
#
#     H = (M + δM)(M + δM)† / 2E + V [diag(1, 0, 0) + ε],
#
# with M = U diag(m₁, m₂, m₃) U† (or U diag(0, Δm²21, Δm²31) U† without
# scalar NSI), δM the scalar NSI and ε the matter NSI; antineutrinos have
# H* and V → -V. Decoherence damps the coherences between the eigenstates of
# H by exp(-Γ_ij (E/E₀)ⁿ L).

"""
    PropagationBlock

A modification of neutrino propagation: [`ScalarNSI`](@ref),
[`MatterNSI`](@ref) or [`Decoherence`](@ref). Combine blocks with
[`NewPhysics`](@ref).
"""
abstract type PropagationBlock end

const FLAVOUR_LABELS = ("e", "μ", "τ")
const ALL_PAIRS = ((1, 1), (1, 2), (1, 3), (2, 2), (2, 3), (3, 3))

pair_label(a, b) = FLAVOUR_LABELS[a] * FLAVOUR_LABELS[b]
ascii_label(a, b) = ("e", "mu", "tau")[a] * ("e", "mu", "tau")[b]

function parse_pairs(couplings)
    isempty(couplings) && return collect(ALL_PAIRS)
    map(collect(couplings)) do c
        s = replace(String(c), "mu" => "μ", "tau" => "τ")
        letters = collect(s)
        length(letters) == 2 && all(in(('e', 'μ', 'τ')), letters) ||
            throw(ArgumentError("unknown flavour pair $c; use :ee, :eμ, :eτ, :μμ, :μτ or :ττ"))
        a, b = (findfirst(==(string(l)), FLAVOUR_LABELS) for l in letters)
        a <= b ? (a, b) : (b, a)
    end
end

"""
    ScalarNSI(couplings...; scale=:Δm31)

Scalar non-standard interactions: a correction δM to the neutrino mass
matrix, δM_αβ = s |η_αβ| e^{iδ_αβ} (Hermitian), with s = √|Δm²31| (the
convention of Ge and Parke, arXiv:1812.08376, and Medhi et al.) or a fixed
scale in eV. Couplings are flavour pairs (`:ee`, `:eμ`, `:eτ`, `:μμ`, `:μτ`,
`:ττ`, or ASCII `:emu`, ...; all six by default). Parameters: `η_αβ` and,
off the diagonal, the phase `δ_αβ`.

Scalar NSI depend on the absolute neutrino masses: the parameter
`m_lightest` (eV, default 0) is the lightest mass.
"""
struct ScalarNSI <: PropagationBlock
    pairs::Vector{Tuple{Int, Int}}
    scale::Union{Symbol, Float64}
    ScalarNSI(couplings...; scale = :Δm31) =
        new(parse_pairs(couplings), scale isa Symbol ? scale : Float64(scale))
end

"""
    MatterNSI(couplings...)

Neutral-current (vector) non-standard interactions with matter, H ∋ V ε,
in units of the standard matter potential V = √2 G_F N_e. Parameters
`ε_αβ` and, off the diagonal, the phase `δε_αβ`.
"""
struct MatterNSI <: PropagationBlock
    pairs::Vector{Tuple{Int, Int}}
    MatterNSI(couplings...) = new(parse_pairs(couplings))
end

"""
    Decoherence(; n=0, E0=1.0, pairs=false)

Quantum decoherence: the coherence between the propagation eigenstates i, j
decays as exp(-Γ_ij (E/E₀)ⁿ L), with Γ in GeV and E₀ in GeV. One parameter
`Γ` for all pairs, or `Γ21`, `Γ31`, `Γ32` with `pairs=true` (the eigenstates
of the Hamiltonian in matter, numbered by increasing eigenvalue).
"""
Base.@kwdef struct Decoherence <: PropagationBlock
    n::Float64 = 0.0
    E0::Float64 = 1.0
    pairs::Bool = false
end

block_parameters(b::ScalarNSI) = [
    p for (a, c) in b.pairs for p in (
        a == c ? [Symbol("η_", pair_label(a, c))] :
        [Symbol("η_", pair_label(a, c)), Symbol("δ_", pair_label(a, c))]
    )
]
block_parameters(b::MatterNSI) = [
    p for (a, c) in b.pairs for p in (
        a == c ? [Symbol("ε_", pair_label(a, c))] :
        [Symbol("ε_", pair_label(a, c)), Symbol("δε_", pair_label(a, c))]
    )
]
block_parameters(b::Decoherence) = b.pairs ? [:Γ21, :Γ31, :Γ32] : [:Γ]

# ASCII names for all of these, accepted wherever parameters are named
for (a, c) in ALL_PAIRS
    PARAM_ALIASES[Symbol("eta_", ascii_label(a, c))] = Symbol("η_", pair_label(a, c))
    PARAM_ALIASES[Symbol("eps_", ascii_label(a, c))] = Symbol("ε_", pair_label(a, c))
    if a != c
        PARAM_ALIASES[Symbol("delta_", ascii_label(a, c))] = Symbol("δ_", pair_label(a, c))
        PARAM_ALIASES[Symbol("deltaeps_", ascii_label(a, c))] = Symbol("δε_", pair_label(a, c))
    end
end
for k in
    (:Gamma => :Γ, :Gamma21 => :Γ21, :Gamma31 => :Γ31, :Gamma32 => :Γ32, :mlightest => :m_lightest)
    PARAM_ALIASES[k.first] = k.second
end

"""
    NewPhysics(blocks...)

A model of neutrino propagation made of [`PropagationBlock`](@ref)s. Use it
as the probability engine with [`use_model!`](@ref) and make parameters with
`OscParams(model; ...)`:

```julia
model = NewPhysics(ScalarNSI(:eτ), Decoherence(n=0))
use_model!(model)
truth = OscParams(model, nufit_params("6.1"); η_eτ=0.1, Γ=1e-23)
probability(truth, :μ, :e; E=2.5, L=1300, ρ=2.848)
```
"""
struct NewPhysics
    blocks::Vector{PropagationBlock}
    names::Vector{Symbol}
end

function NewPhysics(blocks::PropagationBlock...)
    names = Symbol[]
    any(b -> b isa ScalarNSI, blocks) && push!(names, :m_lightest)
    for b in blocks
        append!(names, block_parameters(b))
    end
    allunique(names) || throw(ArgumentError("two blocks define the same parameters"))
    return NewPhysics(collect(PropagationBlock, blocks), names)
end

Base.show(io::IO, m::NewPhysics) =
    print(io, "NewPhysics(", join(map(block_string, m.blocks), ", "), ")")
block_string(b::ScalarNSI) = "ScalarNSI(" * join((pair_label(p...) for p in b.pairs), ", ") * ")"
block_string(b::MatterNSI) = "MatterNSI(" * join((pair_label(p...) for p in b.pairs), ", ") * ")"
block_string(b::Decoherence) = "Decoherence(n=$(b.n)$(b.pairs ? ", per pair" : ""))"

parameter_names(m::NewPhysics) = m.names

"""
    OscParams(model::NewPhysics; θ12, ..., new physics parameters...)
    OscParams(model::NewPhysics, base::OscParams; changes...)

Oscillation parameters for a [`NewPhysics`](@ref) model: the standard ones
and the model's parameters (all zero unless given).
"""
function OscParams(m::NewPhysics; kwargs...)
    kw = Dict(canonical(k) => v for (k, v) in kwargs)
    std = Dict(k => v for (k, v) in kw if k in PARAM_NAMES || k === :ρ)
    extra = NamedTuple(n => Float64(get(kw, n, 0.0)) for n in m.names)
    for k in keys(kw)
        k in PARAM_NAMES ||
            k === :ρ ||
            k in m.names ||
            throw(ArgumentError("$k is not a parameter of $m"))
    end
    return OscParams(; std..., extra)
end

function OscParams(m::NewPhysics, base::OscParams; kwargs...)
    p = OscParams(
        m;
        θ12 = base.θ12,
        θ13 = base.θ13,
        θ23 = base.θ23,
        δ = base.δ,
        Δm21 = base.Δm21,
        Δm31 = base.Δm31,
        ρ = length(base.densities) == 1 ? base.densities[1] : base.densities,
        (k => base.values[6 + i] for (i, k) in enumerate(base.extra_names) if k in m.names)...,
    )
    return isempty(kwargs) ? p : OscParams(p; kwargs...)
end

# ---------------------------------------------------------------------------
#  The probability engine
# ---------------------------------------------------------------------------

mutable struct NewPhysicsEngine
    model::NewPhysics
    index::Dict{Symbol, Int}
    cached::Vector{Float64}          # parameters of the cached vacuum part
    K::Matrix{ComplexF64}            # (M + δM)(M + δM)† in eV²
    ε::Matrix{ComplexF64}
    Γ::NTuple{3, Float64}            # Γ21, Γ31, Γ32 in GeV
    deco::Union{Nothing, Decoherence}
    H::Matrix{ComplexF64}
    Q::Matrix{ComplexF64}
    λ::Vector{Float64}
    S::Matrix{ComplexF64}
    T::Matrix{ComplexF64}
    ρ::Array{ComplexF64, 3}
end

function NewPhysicsEngine(m::NewPhysics)
    z() = zeros(ComplexF64, 3, 3)
    index = Dict(n => 6 + i for (i, n) in enumerate(m.names))
    deco = findfirst(b -> b isa Decoherence, m.blocks)
    NewPhysicsEngine(
        m,
        index,
        Float64[],
        z(),
        z(),
        (0.0, 0.0, 0.0),
        deco === nothing ? nothing : m.blocks[deco],
        z(),
        z(),
        zeros(3),
        z(),
        z(),
        zeros(ComplexF64, 3, 3, 3),
    )
end

function update_vacuum!(eng::NewPhysicsEngine, v::Vector{Float64})
    eng.cached == v && return
    eng.cached = copy(v)
    get_(n) = v[eng.index[n]]
    p = OscParams(v[1:6], Symbol[], [1.0])
    U = pmns_matrix(p)
    scalar = findfirst(b -> b isa ScalarNSI, eng.model.blocks)
    if scalar === nothing
        eng.K .= U * Diagonal([0.0, p.Δm21, p.Δm31]) * U'
    else
        b = eng.model.blocks[scalar]
        m = light_masses(p; lightest = get_(:m_lightest))
        s = b.scale === :Δm31 ? sqrt(abs(p.Δm31)) : b.scale
        δM = zeros(ComplexF64, 3, 3)
        for (a, c) in b.pairs
            η = get_(Symbol("η_", pair_label(a, c)))
            if a == c
                δM[a, a] = s * η
            else
                δM[a, c] = s * η * cis(get_(Symbol("δ_", pair_label(a, c))))
                δM[c, a] = conj(δM[a, c])
            end
        end
        Meff = U * Diagonal(m) * U' + δM
        eng.K .= Meff * Meff'
    end
    fill!(eng.ε, 0)
    for b in eng.model.blocks
        b isa MatterNSI || continue
        for (a, c) in b.pairs
            ε = get_(Symbol("ε_", pair_label(a, c)))
            if a == c
                eng.ε[a, a] = ε
            else
                eng.ε[a, c] = ε * cis(get_(Symbol("δε_", pair_label(a, c))))
                eng.ε[c, a] = conj(eng.ε[a, c])
            end
        end
    end
    if eng.deco !== nothing
        eng.Γ = eng.deco.pairs ? (get_(:Γ21), get_(:Γ31), get_(:Γ32)) : ntuple(_ -> get_(:Γ), 3)
    end
end

# Eigensystem of H = K/2E ± V(e⊗e + ε) for one layer (E in eV, V in eV)
function layer_eigensystem!(eng::NewPhysicsEngine, E, V, anti)
    H = eng.H
    for i = 1:3, j = 1:3
        vacuum = eng.K[i, j] / (2E)
        matter = V * ((i == j == 1) + eng.ε[i, j])
        H[i, j] = anti ? conj(vacuum) - conj(matter) : vacuum + matter
    end
    F = eigen!(Hermitian(H))
    eng.λ .= F.values                     # ascending
    eng.Q .= F.vectors
end

function new_physics_probabilities!(P, eng::NewPhysicsEngine, v, E_GeV, anti, lengths, densities)
    update_vacuum!(eng, v)
    E = E_GeV * 1e9
    if eng.deco === nothing
        S = eng.S
        fill!(S, 0)
        S[1, 1] = S[2, 2] = S[3, 3] = 1
        for k in eachindex(lengths)
            layer_eigensystem!(eng, E, densities[k] * GLB_V_FACTOR * GLB_Ne_MANTLE, anti)
            L = km_to_ev(lengths[k])
            Q = eng.Q
            eng.T .= Q * Diagonal(cis.(-L .* eng.λ)) * Q' * S
            S .= eng.T
        end
        for α = 1:3, β = 1:3
            P[α, β] = abs2(S[β, α])
        end
        return P
    end
    # Density matrices for initial flavour α, evolved layer by layer
    ρ = eng.ρ
    fill!(ρ, 0)
    for α = 1:3
        ρ[α, α, α] = 1
    end
    scale = (E_GeV / eng.deco.E0)^eng.deco.n * 1e9           # Γ in GeV → eV
    for k in eachindex(lengths)
        layer_eigensystem!(eng, E, densities[k] * GLB_V_FACTOR * GLB_Ne_MANTLE, anti)
        L = km_to_ev(lengths[k])
        Q, λ = eng.Q, eng.λ
        D = Matrix{ComplexF64}(undef, 3, 3)
        for i = 1:3, j = 1:3
            Γ = i == j ? 0.0 : eng.Γ[i + j == 3 ? 1 : i + j == 4 ? 2 : 3]
            D[i, j] = cis(-(λ[i] - λ[j]) * L) * exp(-Γ * scale * L)
        end
        for α = 1:3
            r = Q' * view(ρ, :, :, α) * Q
            ρ[:, :, α] .= Q * (r .* D) * Q'
        end
    end
    for α = 1:3, β = 1:3
        P[α, β] = real(ρ[β, β, α])
    end
    return P
end

"""
    use_model!(model::NewPhysics; setup=current_setup())

Use the oscillation probabilities of `model` from now on (in `setup`), for
all experiments. `standard_engine!` switches back (see [`probability_engine!`](@ref)).
"""
function use_model!(m::NewPhysics; setup::Setup = current_setup())
    eng = NewPhysicsEngine(m)
    P! =
        (P, params, E, anti, lengths, densities) ->
            new_physics_probabilities!(P, eng, params.values, E, anti, lengths, densities)
    probability_engine!(P!, m.names; setup)
    return m
end
