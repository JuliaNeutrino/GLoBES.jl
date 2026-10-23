# Neutrino fluxes and cross sections (port of glb_fluxes.c)
#
# Units: energies in GeV, lengths in km, cross sections in 1e-38 cm^2,
# target mass in kt. Fluxes are normalized to a baseline of 1 km.

# C GLoBES uses π = 3.1415, an outdated Avogadro number and rounded lepton
# masses here; these are the CODATA 2018 values.
const FLUX_NORM = 6.02214076 * 1e-12    # 1e-38 / atomic mass unit
const GLB_ELECTRON_MASS = 0.51099895e-3
const MUON_MASS = 0.1056583755
const NFLUX = 500

# Parent muon kinematics for the built-in neutrino factory flux
mutable struct NuFactState
    beam_en::Float64
    rgamma::Float64
    beta::Float64
end
const nufact_state = NuFactState(20.0, 20.0 / MUON_MASS, -sqrt(1 - (MUON_MASS / 20.0)^2))

function glb_set_max_energy(en)
    nufact_state.beam_en = en
    nufact_state.rgamma = en / MUON_MASS
    nufact_state.beta = -sqrt(1 - 1 / nufact_state.rgamma / nufact_state.rgamma)
end

function fastflux_mu(en, baseline)
    rg = nufact_state.rgamma
    beta = nufact_state.beta
    baseline *= 1000
    x_lab = 2 * en / MUON_MASS
    x = x_lab / (rg * (1 - beta))
    xmax = rg * (1 - beta)
    x_lab <= xmax * (1 + 1e-12) ?
    rg * (1 - beta) / baseline / baseline / MUON_MASS * x * x / π * (3 - 2x) : 0.0
end

function fastflux_e(en, baseline)
    rg = nufact_state.rgamma
    beta = nufact_state.beta
    baseline *= 1000
    x_lab = 2 * en / MUON_MASS
    x = x_lab / (rg * (1 - beta))
    xmax = rg * (1 - beta)
    x_lab <= xmax * (1 + 1e-12) ?
    rg * (1 - beta) / baseline / baseline / MUON_MASS * 6 * x * x / π * (1 - x) : 0.0
end

function nufact_flux(data::Flux, en, polarity, l, anti)
    baseline = 1.0
    l == 1 && anti == 1 && polarity == 1 && return fastflux_e(en, baseline)
    l == 2 && anti == -1 && polarity == 1 && return fastflux_mu(en, baseline)
    l == 1 && anti == -1 && polarity == -1 && return fastflux_e(en, baseline)
    l == 2 && anti == 1 && polarity == -1 && return fastflux_mu(en, baseline)
    return 0.0
end

function bb_spectrum(gamma, end_point, enu)
    baseline = 1000.0
    ye = GLB_ELECTRON_MASS / (end_point + GLB_ELECTRON_MASS)
    # Phase space integral of beta decay (C GLoBES has log10 instead of log)
    gye =
        (1 / 60) * (sqrt(1 - ye^2) * (2 - 9ye^2 - 8ye^4) + 15ye^4 * log(ye / (1 - sqrt(1 - ye^2))))
    part = enu / (2 * gamma * (end_point + GLB_ELECTRON_MASS))
    part >= 1 && return 0.0
    root = (1 - part) * (1 - part) - ye * ye
    root <= 0 && return 0.0
    return (1 / π / gye / baseline^2) * (gamma * gamma / enu) * part^3 * (1 - part) * sqrt(root)
end

function bb_flux(data::Flux, en, polarity, l, anti)
    l == 1 && anti == 1 && polarity == 1 && return bb_spectrum(data.gamma, data.end_point, en)
    l == 1 && anti == -1 && polarity == -1 && return bb_spectrum(data.gamma, data.end_point, en)
    return 0.0
end

function builtin_flux_setup!(data::Flux, flx, pl)
    de = 0.0
    # The table covers the spectrum up to its endpoint (C GLoBES cuts it at
    # 500/501 of the endpoint and then reads it with the wrong step size)
    if data.builtin == 1 || data.builtin == 2
        de = data.parent_energy * (1 + sqrt(1 - (MUON_MASS / data.parent_energy)^2)) / 2 / NFLUX
    elseif data.builtin == 3 || data.builtin == 4
        de = 2 * data.gamma * (data.end_point + GLB_ELECTRON_MASS) / NFLUX
    end
    glb_set_max_energy(data.parent_energy)
    st = zeros(501, 7)
    for i = 0:500
        e = i * de
        st[i + 1, 1] = e
        st[i + 1, 2] = flx(data, e, pl, 1, +1)
        st[i + 1, 3] = flx(data, e, pl, 2, +1)
        st[i + 1, 4] = flx(data, e, pl, 3, +1)
        st[i + 1, 5] = flx(data, e, pl, 1, -1)
        st[i + 1, 6] = flx(data, e, pl, 2, -1)
        st[i + 1, 7] = flx(data, e, pl, 3, -1)
    end
    data.storage = st
end

# A flux or cross section table: rows of an energy and six columns. GLoBES
# 3.0 needs 501 (flux) or 1001 (cross section) rows; like GLoBES 3.2, any
# number of rows with increasing energies works.
function load_table(file_name)
    path = glb_find_file(file_name)
    path === nothing && (glb_error("Could not open $file_name"); return zeros(0, 7))
    rows = Vector{Float64}[]
    for line in eachline(path)
        s = lstrip(line)
        (isempty(s) || startswith(s, "#")) && continue
        vals = Float64[]
        for tok in split(s)
            v = tryparse(Float64, tok)
            v === nothing && break
            push!(vals, v)
        end
        length(vals) >= 7 && push!(rows, vals[1:7])
    end
    length(rows) >= 2 && issorted(first.(rows); lt = <=) ||
        (glb_error("Wrong format in file $file_name"); return zeros(0, 7))
    return permutedims(reduce(hcat, rows))
end

# Linear interpolation in column `col` of a table, zero outside its energies
function interpolate_table(st::Matrix{Float64}, x, col)
    n = size(st, 1)
    (n >= 2 && st[1, 1] <= x <= st[n, 1]) || return 0.0
    k = min(searchsortedlast(view(st, :, 1), x), n - 1)
    w = (x - st[k, 1]) / (st[k + 1, 1] - st[k, 1])
    return st[k, col] + w * (st[k + 1, col] - st[k, col])
end

function glb_init_fluxtable!(data::Flux)
    if data.builtin == 1
        builtin_flux_setup!(data, nufact_flux, +1)
    elseif data.builtin == 2
        builtin_flux_setup!(data, nufact_flux, -1)
    elseif data.builtin == 3
        builtin_flux_setup!(data, bb_flux, +1)
    elseif data.builtin == 4
        builtin_flux_setup!(data, bb_flux, -1)
    elseif data.builtin == 0 || data.builtin == GLB_OLD_NORM
        data.storage = load_table(data.file_name)
    end
end

function norm2(type)
    (type == 1 || type == 2 || type == 3 || type == 4) && return FLUX_NORM * 100
    type == GLB_OLD_NORM && return FLUX_NORM * 295 * 295 * 9.92033 * 1e6
    return 1.0
end

flux_column(l, anti) = anti == 1 ? l + 1 : l + 4

"""
    glb_flux_calc(en, baseline, l, anti, flux)

Flux of flavour `l` (1=e, 2=μ, 3=τ), `anti=+1/-1`, at energy `en` and the
given baseline, including running time, power and normalization.
"""
function glb_flux_calc(en, baseline, l, anti, data::Flux)
    res = interpolate_table(data.storage, en, flux_column(l, anti)) / baseline / baseline
    return res * norm2(data.builtin) * data.time * data.target_power * data.stored_muons * data.norm
end

function glb_init_xsectable!(data::XSec)
    data.builtin == 0 && (data.storage = load_table(data.file_name))
end

"""
    glb_xsec_calc(en, l, anti, xsec)

Total cross section (1e-38 cm²) for flavour `l` and `anti=±1` at energy `en`.
"""
function glb_xsec_calc(enl, l, anti, data::XSec)
    part = l + 3 * (1 - anti) ÷ 2 + 1
    res = interpolate_table(data.storage, log10(enl), part)
    return res * enl
end
