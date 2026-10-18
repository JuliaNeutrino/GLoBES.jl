# Matter density profiles (port of glb_profile.c)

include("stacey_data.jl")

const NSTACEY = 500

core_distance(l, bl) = sqrt(l^2 + GLB_EARTH_RADIUS^2 - bl * l)

function matter_profile(r)
    nr = floor(Int, NSTACEY * r / GLB_EARTH_RADIUS)
    return STACEY[nr + 1]
end

function earth_range(baseline)
    baseline < 0 && (glb_error("The baseline may not be negative"); return -1)
    baseline > 2 * GLB_EARTH_RADIUS &&
        (glb_error("The baseline may not exceed 2*GLB_EARTH_RADIUS"); return -1)
    return 0
end

"""
    glbStaceyProfile(baseline, layers) -> (lengths, densities)

Profile with `layers` equidistant steps along the baseline, densities taken
from the Stacey Earth model (as used for `\$profiletype = 2`).
"""
function glbStaceyProfile(baseline::Real, layers::Integer)
    earth_range(baseline) != 0 && return nothing
    delta_l = baseline / layers
    l = fill(Float64(delta_l), layers)
    d = [
        matter_profile(GLB_EARTH_RADIUS - core_distance((i + 0.5) * delta_l, baseline)) for
        i = 0:(layers - 1)
    ]
    return l, d
end

"""
    glbAverageDensityProfile(baseline) -> (lengths, densities)

Single-layer profile with the average Stacey density along the baseline
(as used for `\$profiletype = 1`).
"""
function glbAverageDensityProfile(baseline::Real)
    earth_range(baseline) != 0 && return nothing
    st = 1000
    l, d = glbStaceyProfile(baseline, st)
    s = 0.0
    for i = 1:st
        s += d[i]
    end
    return [Float64(baseline)], [s / st]
end

"""
    glbLoadProfileData(filename) -> (lengths, densities)

Read a two-column file (length in km, density in g/cm³).
"""
function glbLoadProfileData(filename::AbstractString)
    path = glb_find_file(filename)
    path === nothing && (glb_error("Failed to open matter profile file"); return nothing)
    l = Float64[]
    d = Float64[]
    for line in eachline(path)
        f = split(line)
        length(f) >= 2 || continue
        push!(l, parse(Float64, f[1]))
        push!(d, parse(Float64, f[2]))
    end
    return l, d
end
