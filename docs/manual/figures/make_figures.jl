# Compute the data of the figures in the GLoBES.jl manual. Every figure is
# drawn by pgfplots from a table written here to figures/data/.
#
# Run from anywhere with `julia --project=<path to GLoBES.jl> make_figures.jl`.

using GLoBES
using Printf

const OUT = joinpath(@__DIR__, "data")
const AEDL = joinpath(pkgdir(GLoBES), "test", "data")
mkpath(OUT)

# Write columns to a whitespace-separated table with a header line
function table(file, header, columns...)
    open(joinpath(OUT, file), "w") do io
        println(io, join(header, " "))
        for row in zip(columns...)
            println(io, join((@sprintf("%.8g", x) for x in row), " "))
        end
    end
    println("wrote ", file)
end

# Contour lines by marching squares: the segments where the bilinear
# interpolation of z (on the grid xs × ys) crosses `level`, each written as
# two lines followed by an empty line, which pgfplots draws as a separate line.
function write_contour(file, xs, ys, z, level)
    crossing(a, b, za, zb) = a + (level - za) / (zb - za) * (b - a)
    open(joinpath(OUT, file), "w") do io
        println(io, "x y")
        for i in 1:length(xs)-1, j in 1:length(ys)-1
            corners = ((xs[i], ys[j], z[i, j]), (xs[i+1], ys[j], z[i+1, j]),
                       (xs[i+1], ys[j+1], z[i+1, j+1]), (xs[i], ys[j+1], z[i, j+1]))
            points = Tuple{Float64,Float64}[]
            for k in 1:4
                (x1, y1, z1), (x2, y2, z2) = corners[k], corners[k % 4 + 1]
                if (z1 < level) != (z2 < level)
                    push!(points, (crossing(x1, x2, z1, z2), crossing(y1, y2, z1, z2)))
                end
            end
            # Four crossings (a saddle): pair them by the value at the centre
            if length(points) == 4 && (sum(c[3] for c in corners) / 4 < level) != (corners[1][3] < level)
                points = points[[1, 4, 2, 3]]
            end
            for k in 1:2:length(points)-1
                @printf(io, "%.6g %.6g\n%.6g %.6g\n\n", points[k]..., points[k+1]...)
            end
        end
    end
    println("wrote ", basename(file))
end

# The parameters used throughout the manual (normal ordering)
truth = OscParams(θ12=asin(sqrt(0.307)), θ13=asin(sqrt(0.022)), θ23=asin(sqrt(0.546)),
                  δ=1.5π, Δm21=7.53e-5, Δm31=2.5e-3)

# ---------------------------------------------------------------------------
#  Oscillation probabilities at DUNE-like L = 1300 km
# ---------------------------------------------------------------------------
E = range(0.3, 5.0, 471)
δs = (0.0, 0.5π, π, 1.5π)
for (anti, file) in ((false, "prob_nu.dat"), (true, "prob_nubar.dat"))
    P = [probability(OscParams(truth; δ), :μ, :e; E, L=1300, ρ=2.8, antineutrino=anti) for δ in δs]
    vac = probability(truth, :μ, :e; E, L=1300, antineutrino=anti)
    table(file, ("E", "d0", "d90", "d180", "d270", "vac270"), E, P..., vac)
end

# Disappearance and the other appearance channels
Pmm = probability(truth, :μ, :μ; E, L=1300, ρ=2.8)
Pmt = probability(truth, :μ, :τ; E, L=1300, ρ=2.8)
table("prob_mu.dat", ("E", "mumu", "mutau"), E, Pmm, Pmt)

# ---------------------------------------------------------------------------
#  PREM density profiles
# ---------------------------------------------------------------------------
for L in (3000, 7500, 11000)
    lengths, ρ = prem_profile(L, 400)
    x = cumsum(lengths) .- lengths          # start of every layer
    avg = average_profile(L)[2][1]
    table("prem_$L.dat", ("x", "rho", "avg"), [x; L], [ρ; ρ[end]], fill(avg, length(x) + 1))
end

# ---------------------------------------------------------------------------
#  NOvA event rates
# ---------------------------------------------------------------------------
cd(AEDL)
nova = Experiment("NOvA.glb")
simulate!(truth)
for (rule, file) in ((1, "nova_app.dat"), (2, "nova_dis.dat"))
    r = event_rates(nova; rule)
    alt = event_rates(nova, OscParams(truth; δ=0.5π); rule)
    ih = event_rates(nova, OscParams(truth; Δm31=-truth.Δm31 + truth.Δm21); rule)
    table(file, ("E", "w", "signal", "bg", "total", "total_d90", "total_io"),
          r.energy, r.width, r.signal, r.background, r.total, alt.total, ih.total)
end

# δ scan: systematics only, and minimised over θ23, Δm31 and ρ (normal and
# inverted ordering)
δgrid = range(0, 2π, 73)
sysonly = chi2_scan(truth, :δ => δgrid)
fitted = Float64[]; fitted_io = Float64[]
for δ in δgrid
    start = OscParams(truth; δ)
    push!(fitted, marginalize(start; free=(:θ23, :Δm31, :ρ), priors=(ρ=0.05,)).chi2)
    mirror = OscParams(start; Δm31=-truth.Δm31 + truth.Δm21)
    push!(fitted_io, marginalize(mirror; free=(:θ23, :Δm31, :ρ), priors=(ρ=0.05,)).chi2)
end
table("nova_delta.dat", ("delta", "sys", "fit", "fit_io"), δgrid ./ π, sysonly, fitted, fitted_io)

# Energy scale shift on the ν_μ disappearance spectrum
r = event_rates(nova; rule=2)
# The C GLoBES method (glbShiftEnergyScale in glb_sys.c), for comparison
function shift_c(b, rates, emin, emax)
    n = length(rates); t0 = n * emin / (emax - emin)
    out = zeros(n)
    for i in 0:n-1
        d = b * (i + t0 + 0.5) + i; k = floor(Int, d)
        out[i+1] = k < -1 || k > n - 1 ? 0.0 :
                   k == -1 ? (1 + b) * rates[k+2] * (d - k) :
                   k == n - 1 ? (1 + b) * (rates[k+1] * (d - k) + rates[k+1]) :
                   (1 + b) * ((rates[k+2] - rates[k+1]) * (d - k) + rates[k+1])
    end
    return out
end
emin, emax = energy_range(nova)
spec(b) = energy_scale_shift(nova, r.signal, b)
table("escale_spec.dat", ("E", "w", "nominal", "jl_p", "jl_m", "c_p", "c_m"),
      r.energy, r.width, r.signal, spec(0.05), spec(-0.05),
      shift_c(0.05, r.signal, emin, emax), shift_c(-0.05, r.signal, emin, emax))
bs = range(-0.01, 0.01, 201)
χ²(s) = sum(poisson_chi2(r.signal[i], s[i]) for i in eachindex(s))
table("escale_chi.dat", ("b", "jl", "c"), bs,
      [χ²(spec(b)) for b in bs], [χ²(shift_c(b, r.signal, emin, emax)) for b in bs])

# ---------------------------------------------------------------------------
#  T2K: energy smearing and a two-parameter scan with NOvA
# ---------------------------------------------------------------------------
both = Setup()
with_setup(both) do
    t2k = Experiment("T2K.glb")
    simulate!(truth)
    c = channel_names(t2k)
    ch = something(findfirst(==("nu_e_signal"), c), 1)
    sd = t2k.data.smear_data[t2k.data.channels[ch].energy]
    pre = channel_rates(t2k, ch; smeared=false)
    post = channel_rates(t2k, ch)
    table("smear_pre.dat", ("E", "w", "rate"), sd.simbincenter, sd.simbinsize, pre ./ sd.simbinsize)
    table("smear_post.dat", ("E", "w", "rate"), bin_centers(t2k), bin_widths(t2k), post ./ bin_widths(t2k))

    Experiment("NOvA.glb")
    simulate!(truth)
    s2 = range(0.06, 0.115, 56)
    δ2 = range(0, 2π, 73)
    grid = chi2_scan(truth, :sin²2θ13 => s2, :δ => δ2)
    for (level, tag) in ((2.30, "1s"), (6.18, "2s"), (11.83, "3s"))
        write_contour("t2k_nova_$tag.dat", collect(s2), collect(δ2 ./ π), grid, level)
    end
end

# ---------------------------------------------------------------------------
#  Poisson χ² of one bin
# ---------------------------------------------------------------------------
fitrate = range(0, 25, 501)
table("poisson.dat", ("fit", "n0", "n1", "n5", "n10", "gauss10"), fitrate,
      [[poisson_chi2(n, f) for f in fitrate] for n in (0, 1, 5, 10)]..., (fitrate .- 10) .^ 2 ./ 10)

# ---------------------------------------------------------------------------
#  Neutrino factory flux and cross sections
# ---------------------------------------------------------------------------
nf = Experiment("NFstandard.glb"; setup=Setup())
Ef = range(0.5, 50, 199)
table("nf_flux.dat", ("E", "nue", "numubar"), Ef,
      [flux(nf, x; flavour=:e, L=1) for x in Ef], [flux(nf, x; flavour=:μ, antineutrino=true, L=1) for x in Ef])
Ex = exp10.(range(-1, 2, 151))
table("xsec.dat", ("E", "numu", "numubar"), Ex,
      [cross_section(nf, x; flavour=:μ) / x for x in Ex],
      [cross_section(nf, x; flavour=:μ, antineutrino=true) / x for x in Ex])

