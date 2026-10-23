# Neutrino mass mechanisms: from the parameters of a model to the light
# neutrino mass matrix, the oscillation parameters, the Casas-Ibarra
# parametrisation of the Yukawa couplings and lepton flavour violation.
#
# Conventions: the Yukawa matrix Y has one row per heavy state N_k and one
# column per lepton flavour (the Lagrangian term `Yn Et.n.l` is Y[k,α] N_k η l_α), and the
# light neutrino mass matrix is Mν = Yᵀ Λ Y = U* diag(m) U†, with U the PMNS
# matrix including the Majorana phases.

const HIGGS_VEV = 246.21965           # GeV, (√2 G_F)^(-1/2)
const G_FERMI = 1.1663788e-5          # GeV⁻², PDG 2024
const ALPHA_EM = 1 / 137.035999177    # Thomson limit
const MEG_LIMIT = 1.5e-13             # BR(μ → eγ), 90% CL, MEG II (2025)

# BR(ℓα → ℓβ ν ν̄) (PDG 2024), the normalisation of ℓα → ℓβ γ
const LEPTONIC_BR = Dict((:μ, :e) => 1.0, (:τ, :e) => 0.1782, (:τ, :μ) => 0.1739)

"""
    NeutrinoMassMechanism

A way to generate the light neutrino masses, Mν = Yᵀ Λ Y with Λ diagonal:
[`ScotogenicLoop`](@ref) and [`TypeISeesaw`](@ref). Subtypes implement
`loop_factors(mechanism, values)`, the diagonal of Λ in GeV.
"""
abstract type NeutrinoMassMechanism end

"""
    ScotogenicLoop(; yukawa=:Yn, majorana=:Mn, mη2=:mEt2, λ3=:lambda3, λ4=:lambda4, λ5=:lambda5, v=:v)

Radiative neutrino masses of the scotogenic model (E. Ma, hep-ph/0601225):
one loop of the Z₂-odd singlet fermions N_k and the inert doublet η,

    (Mν)_αβ = Σ_k Y_kα Y_kβ M_k/(32π²) [m_R²/(m_R² − M_k²) ln(m_R²/M_k²) − (R → I)],

with m_R,I² = m_η² + (λ₃ + λ₄ ± λ₅) v²/2. The keywords are the names of the
parameters in the model values, by default those of the scotogenic model file.

```julia
mech = ScotogenicLoop()   # parameter names as in the scotogenic model file
mech = ScotogenicLoop(yukawa = :Y, majorana = :MN)   # other names for the Yukawa matrix and the Majorana masses
```
"""
Base.@kwdef struct ScotogenicLoop <: NeutrinoMassMechanism
    yukawa::Symbol = :Yn
    majorana::Symbol = :Mn
    mη2::Symbol = :mEt2
    λ3::Symbol = :lambda3
    λ4::Symbol = :lambda4
    λ5::Symbol = :lambda5
    v::Symbol = :v
end

"""
    TypeISeesaw(; yukawa=:Yv, majorana=:Mn, v=:v)

Type-I seesaw: Mν = -(v²/2) Yᵀ M⁻¹ Y for heavy Majorana masses M.

```julia
mech = TypeISeesaw()   # Yukawa matrix :Yv, Majorana masses :Mn, vev :v
neutrino_mass_matrix(mech, Dict(:Yv => [1e-6 0 0; 0 2e-6 0; 0 0 3e-6], :Mn => [1e10, 2e10, 3e10], :v => 246.22))   # mechanism, values in GeV: Mν in eV
```
"""
Base.@kwdef struct TypeISeesaw <: NeutrinoMassMechanism
    yukawa::Symbol = :Yv
    majorana::Symbol = :Mn
    v::Symbol = :v
end

function param(values, k::Symbol, default = nothing)
    v = lookup(values, k)
    v === nothing || return v
    default === nothing && throw(ArgumentError("no value for parameter $k"))
    return default
end

vev(mech, values) = Float64(param(values, mech.v, HIGGS_VEV))

# The heavy Majorana masses, which must be given in their mass basis
function heavy_masses(values, name)
    M = param(values, name)
    M isa Number && return [Float64(M)]
    M isa AbstractVector && return Float64.(M)
    isdiag(M) || throw(
        ArgumentError("$name must be diagonal (work in the mass basis of the heavy fermions)"),
    )
    return Float64.(real.(diag(M)))
end

"""
    scotogenic_masses(mechanism, values) -> (; mR, mI, mηp)

Masses in GeV of the neutral CP-even and CP-odd and of the charged inert
scalars: m_R,I² = m_η² + (λ₃ + λ₄ ± λ₅)v²/2, m_η+² = m_η² + λ₃v²/2.

```julia
scotogenic_masses(ScotogenicLoop(), values)   # mechanism, parameter values: (mR, mI, mηp) in GeV
```
"""
function scotogenic_masses(mech::ScotogenicLoop, values)
    v = vev(mech, values)
    m2, l3, l4, l5 = (Float64(param(values, k)) for k in (mech.mη2, mech.λ3, mech.λ4, mech.λ5))
    mR2 = m2 + (l3 + l4 + l5) * v^2 / 2
    mI2 = m2 + (l3 + l4 - l5) * v^2 / 2
    mp2 = m2 + l3 * v^2 / 2
    min(mR2, mI2, mp2) > 0 || throw(ArgumentError("inert scalar mass² < 0: the vacuum is unstable"))
    return (mR = sqrt(mR2), mI = sqrt(mI2), mηp = sqrt(mp2))
end

"""
    loop_factors(mechanism, values) -> Vector{Float64}

Λ_k in GeV, such that Mν = Yᵀ diag(Λ) Y.
"""
function loop_factors(mech::ScotogenicLoop, values)
    v = vev(mech, values)
    m2, l3, l4, l5 = (Float64(param(values, k)) for k in (mech.mη2, mech.λ3, mech.λ4, mech.λ5))
    mR2 = m2 + (l3 + l4 + l5) * v^2 / 2
    mI2 = m2 + (l3 + l4 - l5) * v^2 / 2
    return [
        M / (32π^2) * loop_difference(mR2, mI2, M^2) for M in heavy_masses(values, mech.majorana)
    ]
end

function loop_factors(mech::TypeISeesaw, values)
    v = vev(mech, values)
    return [-v^2 / (2M) for M in heavy_masses(values, mech.majorana)]
end

# f(x) = x/(x − M²) ln(x/M²); f(a) − f(b) without the cancellation of the
# two terms when λ₅ (a − b = λ₅v²) is tiny
function loop_difference(a, b, M2)
    f(x) = (t = x / M2; abs(t - 1) < 1e-4 ? 1 + (t - 1) / 2 - (t - 1)^2 / 6 : t * log(t) / (t - 1))
    abs(a - b) > 1e-4 * abs(a + b) && return f(a) - f(b)
    x = (a + b) / 2
    t = x / M2
    u = t - 1
    df = abs(u) < 1e-4 ? (1 / 2 - u / 3 + u^2 / 4) / M2 : (t - 1 - log(t)) / (t - 1)^2 / M2
    return df * (a - b)
end

"""
    neutrino_mass_matrix(mechanism, values) -> Matrix{ComplexF64}

The light neutrino mass matrix Mν = Yᵀ Λ Y in eV.

```julia
Mν = neutrino_mass_matrix(ScotogenicLoop(), values)   # mechanism, parameter values: Mν in eV
```
"""
function neutrino_mass_matrix(mech::NeutrinoMassMechanism, values)
    Y = ComplexF64.(param(values, mech.yukawa))
    Λ = loop_factors(mech, values)
    size(Y, 1) == length(Λ) || throw(
        DimensionMismatch("Yukawa matrix has $(size(Y, 1)) rows for $(length(Λ)) heavy states"),
    )
    return transpose(Y) * Diagonal(Λ) * Y * 1e9
end

# ---------------------------------------------------------------------------
#  From the mass matrix to the oscillation parameters
# ---------------------------------------------------------------------------

"""
    NeutrinoSpectrum

The light neutrinos of a model, from [`neutrino_spectrum`](@ref):
`params` (oscillation parameters), the masses `m` in eV, `ordering`
(`:normal` or `:inverted`), `U` (PMNS matrix with Majorana phases),
`majorana_phases` (α21, α31), `mββ` (the effective mass of neutrinoless
double-β decay), `mβ` (of β decay) and `Σm`, all in eV.
"""
struct NeutrinoSpectrum
    params::OscParams
    m::Vector{Float64}
    ordering::Symbol
    U::Matrix{ComplexF64}
    majorana_phases::Tuple{Float64, Float64}
    mββ::Float64
    mβ::Float64
    Σm::Float64
end

function Base.show(io::IO, ::MIME"text/plain", s::NeutrinoSpectrum)
    p = s.params
    println(io, "Light neutrinos (", s.ordering, " ordering)")
    @printf(io, "  m₁, m₂, m₃ = %.4g, %.4g, %.4g eV   Σm = %.4g eV\n", s.m..., s.Σm)
    @printf(
        io,
        "  sin²θ12 = %.4f   sin²θ13 = %.5f   sin²θ23 = %.4f   δ = %.1f°\n",
        sin(p.θ12)^2,
        sin(p.θ13)^2,
        sin(p.θ23)^2,
        rad2deg(p.δ)
    )
    @printf(io, "  Δm²21 = %.4g eV²   Δm²31 = %.4g eV²\n", p.Δm21, p.Δm31)
    @printf(
        io,
        "  α21 = %.1f°, α31 = %.1f°   mββ = %.4g eV   mβ = %.4g eV",
        rad2deg.(s.majorana_phases)...,
        s.mββ,
        s.mβ
    )
end

"""
    neutrino_spectrum(Mν) -> NeutrinoSpectrum
    neutrino_spectrum(mechanism, values) -> NeutrinoSpectrum

Masses, mixing angles and phases of a light neutrino mass matrix in eV
(Mν = U* diag(m) U†). The mass states are numbered as usual: m₁ < m₂ with
the smaller splitting between them, so that Δm²31 < 0 means inverted
ordering. Phases are in [0, 2π).

```julia
s = neutrino_spectrum(Mν)   # light neutrino mass matrix in eV
s = neutrino_spectrum(ScotogenicLoop(), values)   # mechanism, parameter values
s.params, s.m, s.ordering, s.mββ, s.Σm   # oscillation parameters, masses and observables in eV
```
"""
neutrino_spectrum(mech::NeutrinoMassMechanism, values) =
    neutrino_spectrum(neutrino_mass_matrix(mech, values))

function neutrino_spectrum(Mν::AbstractMatrix)
    Mν = ComplexF64.(Mν)
    issymmetric(Mν) ||
        isapprox(Mν, transpose(Mν); rtol = 1e-10) ||
        throw(ArgumentError("a Majorana mass matrix must be symmetric"))
    λ, V = eigen(Hermitian(Mν' * Mν))
    m = sqrt.(max.(λ, 0.0))                   # ascending
    # Order: 1,2 are the pair with the smaller splitting
    order = (m[2]^2 - m[1]^2) < (m[3]^2 - m[2]^2) ? [1, 2, 3] : [2, 3, 1]
    m, V = m[order], V[:, order]
    # Takagi: make Uᵀ Mν U real and positive
    D = transpose(V) * Mν * V
    U = V * Diagonal([cis(-angle(D[i, i]) / 2) for i = 1:3])
    θ13 = asin(min(1.0, abs(U[1, 3])))
    θ12 = atan(abs(U[1, 2]), abs(U[1, 1]))
    θ23 = atan(abs(U[2, 3]), abs(U[3, 3]))
    s12, c12, s13, c13, s23, c23 = sin(θ12), cos(θ12), sin(θ13), cos(θ13), sin(θ23), cos(θ23)
    # δ from the Jarlskog invariant and |U_μ1|
    J = imag(U[1, 1] * U[2, 2] * conj(U[1, 2]) * conj(U[2, 1]))
    norm_ = s12 * c12 * s23 * c23 * s13
    δ = 0.0
    if norm_ > 1e-12
        cosδ = (abs2(U[2, 1]) - s12^2 * c23^2 - c12^2 * s23^2 * s13^2) / (2norm_)
        sinδ = J / (norm_ * c13^2)
        δ = mod(atan(sinδ, cosδ), 2π)
    end
    # Majorana phases after removing the charged-lepton phase of row e
    βe = angle(U[1, 1])
    phase(x) = (y = mod(x, 2π); y > 2π - 1e-9 ? 0.0 : y)
    α21 = phase(2 * (angle(U[1, 2]) - βe))
    α31 = phase(2 * (angle(U[1, 3]) - βe + δ))
    params = OscParams(
        θ12 = θ12,
        θ13 = θ13,
        θ23 = θ23,
        δ = δ,
        Δm21 = m[2]^2 - m[1]^2,
        Δm31 = m[3]^2 - m[1]^2,
    )
    ordering = m[3] > m[1] ? :normal : :inverted
    mβ = sqrt(sum(abs2(U[1, i]) * m[i]^2 for i = 1:3))
    return NeutrinoSpectrum(params, m, ordering, U, (α21, α31), abs(Mν[1, 1]), mβ, sum(m))
end

"""
    pmns_matrix(params; majorana_phases=(0, 0)) -> Matrix{ComplexF64}

The PMNS matrix in the PDG parametrisation, times diag(1, e^{iα21/2}, e^{iα31/2}).

```julia
U = pmns_matrix(nufit_params("6.1"); majorana_phases = (0, π / 2))   # oscillation parameters, (α21, α31)
```
"""
function pmns_matrix(p::OscParams; majorana_phases = (0.0, 0.0))
    s12, c12 = sincos(p.θ12)
    s13, c13 = sincos(p.θ13)
    s23, c23 = sincos(p.θ23)
    e = cis(p.δ)
    U = ComplexF64[
        c12*c13 s12*c13 s13/e
        -s12*c23-c12*s23*s13*e c12*c23-s12*s23*s13*e s23*c13
        s12*s23-c12*c23*s13*e -c12*s23-s12*c23*s13*e c23*c13
    ]
    return U * Diagonal([1, cis(majorana_phases[1] / 2), cis(majorana_phases[2] / 2)])
end

"""
    light_masses(params; lightest=0) -> Vector{Float64}

m₁, m₂, m₃ in eV from Δm²21, Δm²31 and the lightest mass (m₁ for normal,
m₃ for inverted ordering).

```julia
light_masses(nufit_params("6.1"); lightest = 0.01)   # oscillation parameters, lightest mass in eV: [m1, m2, m3]
```
"""
function light_masses(p::OscParams; lightest::Real = 0.0)
    if p.Δm31 > 0
        m1 = Float64(lightest)
        return [m1, sqrt(m1^2 + p.Δm21), sqrt(m1^2 + p.Δm31)]
    end
    m3 = Float64(lightest)
    m1 = sqrt(m3^2 - p.Δm31)
    return [m1, sqrt(m1^2 + p.Δm21), m3]
end

"""
    orthogonal_matrix(z1, z2, z3) -> Matrix{ComplexF64}

A complex orthogonal matrix R (RᵀR = 1), the product of rotations by the
complex angles z1, z2, z3 in the 23, 13 and 12 planes.

```julia
R = orthogonal_matrix(0.1, 0.2 + 0.1im, 0)   # complex angles in the 23, 13 and 12 planes
```
"""
function orthogonal_matrix(z1, z2, z3)
    rot(i, j, z) = (
        R = Matrix{ComplexF64}(I, 3, 3);
        R[i, i] =  R[j, j] = cos(z) ;
        R[i, j] = sin(z);
        R[j, i] = -sin(z);
        R
    )
    return rot(2, 3, z1) * rot(1, 3, z2) * rot(1, 2, z3)
end

"""
    casas_ibarra(mechanism, values, params; lightest=0, R=I, majorana_phases=(0, 0)) -> Matrix{ComplexF64}

The Yukawa matrix that reproduces the oscillation parameters `params` (for
example the NuFIT best fit), Y = Λ^(-1/2) R m^(1/2) U† (Casas and Ibarra,
hep-ph/0103065; Toma and Vicente, arXiv:1312.2840 for the scotogenic model).
`values` holds the other parameters of the model; `R` is a complex
orthogonal matrix (see [`orthogonal_matrix`](@ref)).

```julia
mech = ScotogenicLoop()   # the mechanism, with the parameter names of the model file
values = Dict(:Mn => [1e3, 2e3, 3e3], :mEt2 => 1e6, :lambda3 => 0.1, :lambda4 => 0.1, :lambda5 => 1e-9, :v => 246.22)   # the other parameters, in GeV
values[:Yn] = casas_ibarra(mech, values, nufit_params("6.1"); lightest = 1e-3)   # mechanism, values, oscillation parameters to reproduce, lightest mass in eV
values[:Yn] = casas_ibarra(mech, values, nufit_params("6.1"); R = orthogonal_matrix(0.3, 0, 0))   # with a complex orthogonal matrix
neutrino_spectrum(mech, values).params   # reproduces NuFIT
```
"""
function casas_ibarra(
    mech::NeutrinoMassMechanism,
    values,
    p::OscParams;
    lightest::Real = 0.0,
    R = I,
    majorana_phases = (0.0, 0.0),
)
    Λ = loop_factors(mech, values) .* 1e9          # eV
    length(Λ) == 3 || throw(ArgumentError("Casas-Ibarra needs three heavy states"))
    m = light_masses(p; lightest)
    U = pmns_matrix(p; majorana_phases)
    return Diagonal(1 ./ sqrt.(complex.(Λ))) *
           (R * Matrix{ComplexF64}(I, 3, 3)) *
           Diagonal(sqrt.(m)) *
           U'
end

# ---------------------------------------------------------------------------
#  Lepton flavour violation
# ---------------------------------------------------------------------------

"""
    br_lfv(mechanism::ScotogenicLoop, values, from=:μ, to=:e) -> Float64

BR(ℓ_from → ℓ_to γ) in the scotogenic model at one loop (Toma and Vicente,
arXiv:1312.2840):

    BR = 3α/(64π G_F² m_η+⁴) |Σ_k Y_k,to* Y_k,from F₂(M_k²/m_η+²)|² BR(ℓ_from → ℓ_to ν ν̄).

The MEG II limit on BR(μ → eγ) is `GLoBES.MEG_LIMIT` = 1.5 × 10⁻¹³.

```julia
br_lfv(ScotogenicLoop(), values, :μ, :e)   # mechanism, values with the Yukawa matrix, from, to: BR(μ → eγ)
br_lfv(ScotogenicLoop(), values, :τ, :μ)   # BR(τ → μγ)
```
"""
function br_lfv(mech::ScotogenicLoop, values, from = :μ, to = :e)
    a, b = flavour(from), flavour(to)
    key = (FLAVOUR_NAMES[a], FLAVOUR_NAMES[b])
    haskey(LEPTONIC_BR, key) || throw(ArgumentError("no decay $from → $to γ"))
    Y = ComplexF64.(param(values, mech.yukawa))
    mp = scotogenic_masses(mech, values).mηp
    M = heavy_masses(values, mech.majorana)
    amp = sum(conj(Y[k, b]) * Y[k, a] * F2(M[k]^2 / mp^2) for k in eachindex(M))
    return 3ALPHA_EM / (64π * G_FERMI^2 * mp^4) * abs2(amp) * LEPTONIC_BR[key]
end

const FLAVOUR_NAMES = (:e, :μ, :τ)

# The loop function of ℓα → ℓβ γ; F₂(1) = 1/12
function F2(x)
    abs(x - 1) < 1e-3 && return 1 / 12 - (x - 1) / 20 + (x - 1)^2 / 30
    return (1 - 6x + 3x^2 + 2x^3 - 6x^2 * log(x)) / (6 * (1 - x)^4)
end
