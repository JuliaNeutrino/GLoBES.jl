# NuFIT global fits of three-flavour oscillations (www.nu-fit.org): the
# published best fits and ranges, and the Δχ² tables as priors.

# Each entry: (best fit, +1σ, -1σ, 3σ low, 3σ high); δ in degrees, Δm²21 in
# 1e-5 eV², Δm²3ℓ in 1e-3 eV² (Δm²31 for normal, Δm²32 for inverted ordering),
# and the Δχ² of the ordering. From the parameter tables of each release.
const NUFIT_TABLES = Dict(
    ("6.1", false, :normal) => (
        s12 = (0.3088, 0.0067, 0.0066, 0.2893, 0.3295),
        t12 = (33.76, 0.42, 0.41, 32.54, 35.03),
        s23 = (0.470, 0.017, 0.014, 0.432, 0.587),
        t23 = (43.27, 1.0, 0.82, 41.11, 50.02),
        s13 = (0.02249, 0.00057, 0.00057, 0.02070, 0.02420),
        t13 = (8.62, 0.11, 0.11, 8.27, 8.95),
        dcp = (207.0, 23.0, 20.0, 114.0, 405.0),
        dm21 = (7.537, 0.094, 0.10, 7.236, 7.823),
        dm3l = (2.521, 0.026, 0.018, 2.454, 2.592),
        chi2 = 0.0,
    ),
    ("6.1", false, :inverted) => (
        s12 = (0.3088, 0.0067, 0.0066, 0.2893, 0.3295),
        t12 = (33.76, 0.42, 0.41, 32.54, 35.03),
        s23 = (0.555, 0.013, 0.016, 0.437, 0.590),
        t23 = (48.15, 0.75, 0.92, 41.40, 50.21),
        s13 = (0.02261, 0.00056, 0.00056, 0.02091, 0.02433),
        t13 = (8.65, 0.11, 0.11, 8.31, 8.97),
        dcp = (283.0, 24.0, 28.0, 202.0, 347.0),
        dm21 = (7.537, 0.094, 0.10, 7.236, 7.822),
        dm3l = (-2.500, 0.024, 0.023, -2.569, -2.430),
        chi2 = 1.5,
    ),
    ("6.1", true, :normal) => (
        s12 = (0.3088, 0.0067, 0.0066, 0.2893, 0.3295),
        t12 = (33.76, 0.42, 0.41, 32.54, 35.03),
        s23 = (0.470, 0.017, 0.014, 0.435, 0.584),
        t23 = (43.29, 0.96, 0.79, 41.27, 49.86),
        s13 = (0.02248, 0.00055, 0.00059, 0.02064, 0.02418),
        t13 = (8.62, 0.11, 0.11, 8.26, 8.95),
        dcp = (212.0, 26.0, 36.0, 125.0, 365.0),
        dm21 = (7.537, 0.094, 0.10, 7.236, 7.823),
        dm3l = (2.511, 0.021, 0.020, 2.450, 2.576),
        chi2 = 0.0,
    ),
    ("6.1", true, :inverted) => (
        s12 = (0.3088, 0.0067, 0.0066, 0.2893, 0.3295),
        t12 = (33.76, 0.42, 0.41, 32.54, 35.03),
        s23 = (0.550, 0.013, 0.016, 0.439, 0.584),
        t23 = (47.90, 0.73, 0.92, 41.51, 49.83),
        s13 = (0.02262, 0.00057, 0.00056, 0.02093, 0.02441),
        t13 = (8.65, 0.11, 0.11, 8.32, 8.99),
        dcp = (274.0, 22.0, 25.0, 203.0, 335.0),
        dm21 = (7.537, 0.094, 0.10, 7.236, 7.822),
        dm3l = (-2.483, 0.020, 0.020, -2.547, -2.421),
        chi2 = 5.9,
    ),
    ("6.0", false, :normal) => (
        s12 = (0.307, 0.012, 0.011, 0.275, 0.345),
        t12 = (33.68, 0.73, 0.70, 31.63, 35.95),
        s23 = (0.561, 0.012, 0.015, 0.430, 0.596),
        t23 = (48.5, 0.7, 0.9, 41.0, 50.5),
        s13 = (0.02195, 0.00054, 0.00058, 0.02023, 0.02376),
        t13 = (8.52, 0.11, 0.11, 8.18, 8.87),
        dcp = (177.0, 19.0, 20.0, 96.0, 422.0),
        dm21 = (7.49, 0.19, 0.19, 6.92, 8.05),
        dm3l = (2.534, 0.025, 0.023, 2.463, 2.606),
        chi2 = 0.6,
    ),
    ("6.0", false, :inverted) => (
        s12 = (0.308, 0.012, 0.011, 0.275, 0.345),
        t12 = (33.68, 0.73, 0.70, 31.63, 35.95),
        s23 = (0.562, 0.012, 0.015, 0.437, 0.597),
        t23 = (48.6, 0.7, 0.9, 41.4, 50.6),
        s13 = (0.02224, 0.00056, 0.00057, 0.02053, 0.02397),
        t13 = (8.58, 0.11, 0.11, 8.24, 8.91),
        dcp = (285.0, 25.0, 28.0, 201.0, 348.0),
        dm21 = (7.49, 0.19, 0.19, 6.92, 8.05),
        dm3l = (-2.510, 0.024, 0.025, -2.584, -2.438),
        chi2 = 0.0,
    ),
    ("6.0", true, :normal) => (
        s12 = (0.308, 0.012, 0.011, 0.275, 0.345),
        t12 = (33.68, 0.73, 0.70, 31.63, 35.95),
        s23 = (0.470, 0.017, 0.013, 0.435, 0.585),
        t23 = (43.3, 1.0, 0.8, 41.3, 49.9),
        s13 = (0.02215, 0.00056, 0.00058, 0.02030, 0.02388),
        t13 = (8.56, 0.11, 0.11, 8.19, 8.89),
        dcp = (212.0, 26.0, 41.0, 124.0, 364.0),
        dm21 = (7.49, 0.19, 0.19, 6.92, 8.05),
        dm3l = (2.513, 0.021, 0.019, 2.451, 2.578),
        chi2 = 0.0,
    ),
    ("6.0", true, :inverted) => (
        s12 = (0.308, 0.012, 0.011, 0.275, 0.345),
        t12 = (33.68, 0.73, 0.70, 31.63, 35.95),
        s23 = (0.550, 0.012, 0.015, 0.440, 0.584),
        t23 = (47.9, 0.7, 0.9, 41.5, 49.8),
        s13 = (0.02231, 0.00056, 0.00056, 0.02060, 0.02409),
        t13 = (8.59, 0.11, 0.11, 8.25, 8.93),
        dcp = (274.0, 22.0, 25.0, 201.0, 335.0),
        dm21 = (7.49, 0.19, 0.19, 6.92, 8.05),
        dm3l = (-2.484, 0.020, 0.020, -2.547, -2.421),
        chi2 = 6.1,
    ),
    ("5.3", false, :normal) => (
        s12 = (0.307, 0.012, 0.011, 0.275, 0.344),
        t12 = (33.66, 0.73, 0.70, 31.60, 35.94),
        s23 = (0.572, 0.018, 0.023, 0.407, 0.620),
        t23 = (49.1, 1.0, 1.3, 39.6, 51.9),
        s13 = (0.02203, 0.00056, 0.00058, 0.02029, 0.02391),
        t13 = (8.54, 0.11, 0.11, 8.19, 8.89),
        dcp = (197.0, 41.0, 25.0, 108.0, 404.0),
        dm21 = (7.41, 0.21, 0.20, 6.81, 8.03),
        dm3l = (2.511, 0.027, 0.027, 2.428, 2.597),
        chi2 = 0.0,
    ),
    ("5.3", false, :inverted) => (
        s12 = (0.307, 0.012, 0.011, 0.275, 0.344),
        t12 = (33.67, 0.73, 0.71, 31.61, 35.94),
        s23 = (0.578, 0.016, 0.021, 0.412, 0.623),
        t23 = (49.5, 0.9, 1.2, 39.9, 52.1),
        s13 = (0.02219, 0.00059, 0.00057, 0.02047, 0.02396),
        t13 = (8.57, 0.11, 0.11, 8.23, 8.90),
        dcp = (286.0, 27.0, 32.0, 192.0, 360.0),
        dm21 = (7.41, 0.21, 0.20, 6.81, 8.03),
        dm3l = (-2.498, 0.032, 0.024, -2.581, -2.409),
        chi2 = 2.3,
    ),
    ("5.3", true, :normal) => (
        s12 = (0.307, 0.012, 0.011, 0.275, 0.344),
        t12 = (33.67, 0.73, 0.71, 31.61, 35.94),
        s23 = (0.454, 0.019, 0.016, 0.411, 0.606),
        t23 = (42.3, 1.1, 0.9, 39.9, 51.1),
        s13 = (0.02224, 0.00056, 0.00057, 0.02047, 0.02397),
        t13 = (8.58, 0.11, 0.11, 8.23, 8.91),
        dcp = (232.0, 39.0, 25.0, 139.0, 350.0),
        dm21 = (7.41, 0.21, 0.20, 6.81, 8.03),
        dm3l = (2.505, 0.024, 0.026, 2.426, 2.586),
        chi2 = 0.0,
    ),
    ("5.3", true, :inverted) => (
        s12 = (0.307, 0.012, 0.011, 0.275, 0.344),
        t12 = (33.67, 0.73, 0.71, 31.61, 35.94),
        s23 = (0.568, 0.016, 0.021, 0.412, 0.611),
        t23 = (48.9, 0.9, 1.2, 39.9, 51.4),
        s13 = (0.02222, 0.00069, 0.00057, 0.02049, 0.02420),
        t13 = (8.57, 0.13, 0.11, 8.23, 8.95),
        dcp = (273.0, 24.0, 26.0, 195.0, 342.0),
        dm21 = (7.41, 0.21, 0.20, 6.81, 8.03),
        dm3l = (-2.487, 0.027, 0.024, -2.566, -2.407),
        chi2 = 9.1,
    ),
)

const NUFIT_DATASETS = Dict(
    "6.1" =>
        ("IC23 without SK atmospheric data", "IC24 with SK atmospheric data", "November 2025"),
    "6.0" =>
        ("IC19 without SK atmospheric data", "IC24 with SK atmospheric data", "September 2024"),
    "5.3" => ("without SK atmospheric data", "with SK atmospheric data", "March 2024"),
)

const NUFIT_LATEST = "6.1"

"""
    NuFit

One column of a NuFIT parameter table: best fit, ±1σ and 3σ range of every
oscillation parameter for one release, ordering and data set. Fields are
named tuples `(best, plus, minus, lo, hi)`: `s12`, `s23`, `s13` (sin²θ),
`t12`, `t23`, `t13` (θ in degrees), `dcp` (δ in degrees), `dm21` (Δm²21 in
10⁻⁵ eV²), `dm3l` (Δm²3ℓ in 10⁻³ eV²), and `chi2`, the Δχ² of the ordering.
"""
struct NuFit
    version::String
    ordering::Symbol
    sk::Bool
    table::NamedTuple
end

"""
    nufit(version="6.1"; ordering=:best, sk=false) -> NuFit

The published NuFIT best fit and ranges (versions "5.3", "6.0", "6.1").
`ordering` is `:normal`, `:inverted` or `:best` (the better fit of the
two); `sk=true` selects the analysis with the Super-Kamiokande atmospheric
χ² table.

```julia
nf = nufit("6.1"; ordering = :normal)   # version, ordering (:normal, :inverted or :best); sk = true with SK atmospheric data
nf.s23   # (best = 0.47, plus = 0.017, minus = 0.014, lo = 0.432, hi = 0.587)
```
"""
function nufit(version::AbstractString = NUFIT_LATEST; ordering::Symbol = :best, sk::Bool = false)
    version = String(version)
    haskey(NUFIT_DATASETS, version) || throw(
        ArgumentError(
            "NuFIT $version is not built in; available: $(join(sort(collect(keys(NUFIT_DATASETS))), ", "))",
        ),
    )
    ordering =
        ordering === :NO || ordering === :NH ? :normal :
        ordering === :IO || ordering === :IH ? :inverted : ordering
    if ordering === :best
        ordering =
            NUFIT_TABLES[(version, sk, :normal)].chi2 <=
            NUFIT_TABLES[(version, sk, :inverted)].chi2 ? :normal : :inverted
    end
    ordering in (:normal, :inverted) ||
        throw(ArgumentError("ordering must be :normal, :inverted or :best"))
    raw = NUFIT_TABLES[(version, sk, ordering)]
    named(t) = (best = t[1], plus = t[2], minus = t[3], lo = t[4], hi = t[5])
    table = (; (k => (k === :chi2 ? v : named(v)) for (k, v) in pairs(raw))...)
    return NuFit(version, ordering, sk, table)
end

Base.getproperty(nf::NuFit, k::Symbol) =
    k in fieldnames(NuFit) ? getfield(nf, k) : getproperty(getfield(nf, :table), k)
Base.propertynames(nf::NuFit) = (fieldnames(NuFit)..., keys(getfield(nf, :table))...)

dataset(nf::NuFit) = NUFIT_DATASETS[nf.version][nf.sk ? 2 : 1]

function Base.show(io::IO, nf::NuFit)
    print(io, "NuFIT ", nf.version, " (", nf.ordering, " ordering, ", dataset(nf), ")")
end

function Base.show(io::IO, ::MIME"text/plain", nf::NuFit)
    println(
        io,
        "NuFIT ",
        nf.version,
        ", ",
        nf.ordering,
        " ordering (Δχ² = ",
        nf.chi2,
        "), ",
        dataset(nf),
    )
    rows = (
        ("sin²θ12", nf.s12, ""),
        ("sin²θ23", nf.s23, ""),
        ("sin²θ13", nf.s13, ""),
        ("δ", nf.dcp, "°"),
        ("Δm²21", nf.dm21, " × 10⁻⁵ eV²"),
        (nf.ordering === :normal ? "Δm²31" : "Δm²32", nf.dm3l, " × 10⁻³ eV²"),
    )
    for (i, (name, v, unit)) in enumerate(rows)
        @printf(
            io,
            "  %-8s %8.5g +%.3g -%.3g   3σ: %.5g → %.5g%s",
            name,
            v.best,
            v.plus,
            v.minus,
            v.lo,
            v.hi,
            unit
        )
        i < length(rows) && println(io)
    end
end

"""
    nufit_params(version="6.1"; ordering=:best, sk=false) -> OscParams
    nufit_params(nf::NuFit) -> OscParams

The NuFIT best fit as oscillation parameters (Δm²31 = Δm²32 + Δm²21 for
inverted ordering).

```julia
truth = nufit_params("6.1"; ordering = :inverted)   # version, ordering: the best fit as OscParams
truth = nufit_params(nufit("6.1"))   # from a NuFIT table
```
"""
nufit_params(version::AbstractString = NUFIT_LATEST; kwargs...) =
    nufit_params(nufit(version; kwargs...))

function nufit_params(nf::NuFit)
    Δm21 = nf.dm21.best * 1e-5
    Δm3l = nf.dm3l.best * 1e-3
    return OscParams(
        θ12 = asin(sqrt(nf.s12.best)),
        θ13 = asin(sqrt(nf.s13.best)),
        θ23 = asin(sqrt(nf.s23.best)),
        δ = deg2rad(nf.dcp.best),
        Δm21 = Δm21,
        Δm31 = nf.ordering === :normal ? Δm3l : Δm3l + Δm21,
    )
end

"""
    nufit_priors(nf::NuFit; params=(:θ12, :θ13, :Δm21)) -> NamedTuple

Gaussian 1σ errors (the mean of the + and − errors) for the `priors` of
[`marginalize`](@ref), in radians and eV². The Gaussian is good for θ12, θ13
and the mass splittings; for θ23 and δ use [`nufit_prior`](@ref) with the
Δχ² tables instead.

```julia
priors = nufit_priors(nufit("6.1"); params = (:θ12, :θ13, :Δm21))   # NuFIT table, parameters: 1σ errors for marginalize
```
"""
function nufit_priors(nf::NuFit; params = (:θ12, :θ13, :Δm21))
    σ(v) = (v.plus + v.minus) / 2
    all_ = (
        θ12 = deg2rad(σ(nf.t12)),
        θ13 = deg2rad(σ(nf.t13)),
        θ23 = deg2rad(σ(nf.t23)),
        δ = deg2rad(σ(nf.dcp)),
        Δm21 = σ(nf.dm21) * 1e-5,
        Δm31 = σ(nf.dm3l) * 1e-3,
    )
    return NamedTuple(canonical(k) => all_[canonical(k)] for k in params)
end

"""
    nufit_ranges(nf::NuFit; sigma=3) -> NamedTuple

The 1σ or 3σ ranges `(lo, hi)` of the parameters, in radians and eV², for
scan limits.

```julia
nufit_ranges(nufit("6.1"); sigma = 3)   # NuFIT table, 1 or 3: the range of each parameter
```
"""
function nufit_ranges(nf::NuFit; sigma::Integer = 3)
    sigma in (1, 3) || throw(ArgumentError("NuFIT publishes 1σ and 3σ ranges"))
    r(v) = sigma == 3 ? (v.lo, v.hi) : (v.best - v.minus, v.best + v.plus)
    deg(v) = deg2rad.(r(v))
    off = nf.ordering === :inverted ? nf.dm21.best * 1e-5 : 0.0
    return (
        θ12 = deg(nf.t12),
        θ13 = deg(nf.t13),
        θ23 = deg(nf.t23),
        δ = deg(nf.dcp),
        Δm21 = r(nf.dm21) .* 1e-5,
        Δm31 = r(nf.dm3l) .* 1e-3 .+ off,
    )
end

# ---------------------------------------------------------------------------
#  Δχ² tables
# ---------------------------------------------------------------------------

# Projection names in the NuFIT files and their parameters
const NUFIT_AXES =
    Dict(:θ12 => "T12", :θ13 => "T13", :θ23 => "T23", :δ => "DCP", :Δm21 => "DMS", :Δm31 => "DMA")

"""
    NuFitChi2

Δχ² projections of a NuFIT release, from [`nufit_chi2`](@ref). Index it
with the parameter names: `t[:δ]` is the one-dimensional projection
`(x, Δχ²)`, `t[:θ23, :δ]` the two-dimensional `(x, y, Δχ²)` for normal
ordering (`t[:θ23, :δ; ordering=:inverted]` for inverted). The axes are
sin²θ, δ in degrees, Δm²21 and Δm²3ℓ in eV². Call `t(params)` for the
Δχ² of a parameter point (see [`nufit_prior`](@ref)).
"""
struct NuFitChi2
    version::String
    sk::Bool
    tables::Dict{Tuple{Symbol, String}, Any}     # (ordering, "T23/DCP") => grid
    offset::Dict{Symbol, Float64}                # Δχ² of each ordering
end

Base.show(io::IO, t::NuFitChi2) = print(
    io,
    "NuFitChi2(NuFIT ",
    t.version,
    t.sk ? ", with SK" : "",
    ", ",
    length(t.tables) ÷ 2,
    " projections per ordering)",
)

struct Grid{N}
    axes::NTuple{N, Vector{Float64}}
    values::Array{Float64, N}
end

"""
    nufit_chi2(version="6.1"; sk=false, full=false, files=nothing) -> NuFitChi2

Load the NuFIT Δχ² tables. The one-dimensional projections of NuFIT 6.0 and
6.1 come with GLoBES.jl; `full=true` downloads the complete release files
(about 4 MB each) with the two- and three-dimensional projections from
www.nu-fit.org and keeps them in the Julia depot. `files` reads your own
copies instead: a pair `(normal, inverted)` of paths, `.txt` or `.txt.xz`.

```julia
t = nufit_chi2("6.1")   # version: the one-dimensional Δχ² projections that come with GLoBES.jl
t = nufit_chi2("6.1"; full = true)   # download the full tables, with two- and three-dimensional projections
x, χ² = t[:δ]   # a projection: sin²θ, δ in degrees, Δm² in eV²
```
"""
function nufit_chi2(
    version::AbstractString = NUFIT_LATEST;
    sk::Bool = false,
    full::Bool = false,
    files = nothing,
)
    version = String(version)
    tag = "v" * replace(version, "." => "") * (sk ? "-TByes" : "-TBoff")
    paths = if files !== nothing
        Dict(:normal => files[1], :inverted => files[2])
    elseif full
        Dict(o => nufit_download(version; sk, ordering = o) for o in (:normal, :inverted))
    else
        dir = joinpath(pkgdir(GLoBES), "data", "nufit")
        p = Dict(
            :normal => joinpath(dir, "$tag-NO-1d.txt"),
            :inverted => joinpath(dir, "$tag-IO-1d.txt"),
        )
        isfile(p[:normal]) || throw(
            ArgumentError(
                "the Δχ² tables of NuFIT $version are not built in; use full=true to download them",
            ),
        )
        p
    end
    tables = Dict{Tuple{Symbol, String}, Any}()
    offset = Dict{Symbol, Float64}()
    for (ordering, path) in paths
        for (name, grid) in read_nufit_file(path)
            tables[(ordering, name)] = grid
        end
        offset[ordering] = minimum(tables[(ordering, "T13")].values)
    end
    return NuFitChi2(version, sk, tables, offset)
end

nufit_cache() =
    joinpath(first(DEPOT_PATH), "scratchspaces", "367f8d2e-8d77-4eed-9bd4-73b32ac618fb", "nufit")

"""
    nufit_download(version="6.1"; sk=false, ordering=:normal) -> path

Download (once) and unpack a NuFIT release file; returns the path of the
text file.

```julia
nufit_download("6.1"; ordering = :inverted)   # version, ordering: the path of the downloaded table
```
"""
function nufit_download(
    version::AbstractString = NUFIT_LATEST;
    sk::Bool = false,
    ordering::Symbol = :normal,
)
    o = ordering === :normal ? "NO" : "IO"
    name = "v" * replace(version, "." => "") * ".release-TB" * (sk ? "yes" : "off") * "-$o.txt"
    dir = mkpath(nufit_cache())
    txt = joinpath(dir, name)
    isfile(txt) && return txt
    xz = txt * ".xz"
    url = "http://www.nu-fit.org/sites/default/files/$name.xz"
    @info "Downloading $url"
    Downloads.download(url, xz)
    run(pipeline(`$(p7zip_jll.p7zip()) x -y -o$dir $xz`; stdout = devnull))
    rm(xz)
    return txt
end

# Read the projections of a NuFIT release file (plain or .xz)
function read_nufit_file(path::AbstractString)
    if endswith(path, ".xz")
        dir = mktempdir()
        run(pipeline(`$(p7zip_jll.p7zip()) x -y -o$dir $path`; stdout = devnull))
        path = joinpath(dir, only(readdir(dir)))
    end
    out = Dict{String, Grid}()
    name, rows = nothing, Vector{Vector{Float64}}()
    function finish!()
        name === nothing || isempty(rows) || (out[name] = make_grid(rows))
        empty!(rows)
    end
    for line in eachline(path)
        if startswith(line, "#")
            m = match(r"^#\s*(\S+) projection", line)
            m === nothing && continue
            finish!()
            name = m.captures[1]
            rows = Vector{Vector{Float64}}()
        elseif name !== nothing
            nums = [tryparse(Float64, w) for w in split(line)]
            isempty(nums) || any(isnothing, nums) || push!(rows, Float64.(nums))
        end
    end
    finish!()
    return out
end

function make_grid(rows)
    n = length(rows[1]) - 1
    axes = ntuple(k -> sort!(unique([r[k] for r in rows])), n)
    values = fill(NaN, map(length, axes)...)
    index = ntuple(k -> Dict(x => i for (i, x) in enumerate(axes[k])), n)
    for r in rows
        values[(index[k][r[k]] for k = 1:n)...] = r[end]
    end
    return Grid{n}(axes, values)
end

# Multilinear interpolation; linear extrapolation outside a 1D table, the
# nearest edge outside 2D and 3D tables
function interpolate(g::Grid{N}, x::NTuple{N, Real}) where {N}
    cell = map(g.axes, x) do ax, xi
        i = clamp(searchsortedlast(ax, xi), 1, length(ax) - 1)
        t = (xi - ax[i]) / (ax[i + 1] - ax[i])
        N > 1 && (t = clamp(t, 0.0, 1.0))
        (i, t)
    end
    acc = 0.0
    for corner in Iterators.product(ntuple(_ -> (0, 1), N)...)
        w = prod(k -> corner[k] == 1 ? cell[k][2] : 1 - cell[k][2], 1:N)
        idx = ntuple(k -> cell[k][1] + corner[k], N)
        w == 0 || (acc += w * g.values[idx...])
    end
    return acc
end

# The table coordinate of a parameter
function nufit_coordinate(name::Symbol, p::OscParams)
    name === :θ12 && return sin(p.θ12)^2
    name === :θ13 && return sin(p.θ13)^2
    name === :θ23 && return sin(p.θ23)^2
    name === :Δm21 && return log10(p.Δm21)
    name === :Δm31 && return (p.Δm31 > 0 ? p.Δm31 : p.Δm31 - p.Δm21) * 1e3
    name === :δ && return rad2deg(p.δ)
    throw(ArgumentError("NuFIT has no projection on $name"))
end

function table_δ(g::Grid, k, x)
    lo, hi = first(g.axes[k]), last(g.axes[k])
    return lo + mod(x - lo, hi - lo)
end

function lookup_grid(t::NuFitChi2, names, ordering)
    key = (ordering, join((NUFIT_AXES[canonical(n)] for n in names), "/"))
    haskey(t.tables, key) && return t.tables[key], false
    rkey = (ordering, join((NUFIT_AXES[canonical(n)] for n in reverse(names)), "/"))
    haskey(t.tables, rkey) && return t.tables[rkey], true
    throw(
        ArgumentError(
            "no NuFIT projection on $(join(names, ", ")) (load full=true for 2D and 3D tables)",
        ),
    )
end

# The axis of a table in natural units
function natural_axis(name::Symbol, ax)
    name === :Δm21 && return exp10.(ax)
    name === :Δm31 && return ax .* 1e-3
    return copy(ax)
end

function Base.getindex(t::NuFitChi2, names::Symbol...; ordering::Symbol = :normal)
    names = canonical.(names)
    g, reversed = lookup_grid(t, names, ordering)
    if length(names) == 1
        return (natural_axis(names[1], g.axes[1]), copy(g.values))
    end
    axes = reversed ? reverse(g.axes) : g.axes
    vals = reversed ? permutedims(g.values, reverse(1:ndims(g.values))) : g.values
    return (map(natural_axis, names, axes)..., copy(vals))
end

"""
    (t::NuFitChi2)(params; use=(:θ12, :θ13, :θ23, :δ, :Δm21, :Δm31)) -> Float64

The NuFIT Δχ² of the oscillation parameters `params`, for the ordering given
by the sign of Δm²31. With the full tables, θ23, δ and Δm²3ℓ come from the
three-dimensional projection; otherwise, and for the other parameters, the
one-dimensional projections are added. Adding projections ignores the
correlations between parameters, which is a good approximation for θ12, θ13
and Δm²21. Only the parameters in `use` contribute.
"""
function (t::NuFitChi2)(p::OscParams; use = (:θ12, :θ13, :θ23, :δ, :Δm21, :Δm31))
    ordering = p.Δm31 > 0 ? :normal : :inverted
    off = t.offset[ordering]
    use = Set(canonical.(use))
    χ² = off
    key3 = (ordering, "T23/DMA/DCP")
    if haskey(t.tables, key3) && :θ23 in use && :δ in use && :Δm31 in use
        g = t.tables[key3]
        χ² +=
            interpolate(
                g,
                (
                    nufit_coordinate(:θ23, p),
                    nufit_coordinate(:Δm31, p),
                    table_δ(g, 3, rad2deg(p.δ)),
                ),
            ) - off
        setdiff!(use, (:θ23, :δ, :Δm31))
    end
    for name in use
        g = t.tables[(ordering, NUFIT_AXES[name])]
        x = name === :δ ? table_δ(g, 1, rad2deg(p.δ)) : nufit_coordinate(name, p)
        χ² += interpolate(g, (x,)) - off
    end
    return χ²
end

"""
    nufit_prior(t::NuFitChi2; use=(:θ12, :θ13, :Δm21)) -> Function

A prior for [`marginalize`](@ref)'s `prior` keyword: the NuFIT Δχ² of the
parameters in `use` (see `NuFitChi2`), relative to its minimum.

```julia
t = nufit_chi2("6.1")   # the NuFIT Δχ² tables
marginalize(test; fixed = (:δ,), prior = nufit_prior(t; use = (:θ12, :θ13, :Δm21, :θ23)))   # parameters whose NuFIT Δχ² is added
```
"""
function nufit_prior(t::NuFitChi2; use = (:θ12, :θ13, :Δm21))
    ref = min(t.offset[:normal], t.offset[:inverted])
    return p -> t(p; use) - ref
end
