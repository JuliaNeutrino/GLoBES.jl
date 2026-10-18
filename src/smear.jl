# Energy resolution / migration matrices (port of glb_smear.c)

const GLB_TYPE_A = 1
const GLB_TYPE_C = 2
const SMEAR_TOL = 1.0e-7

# sigma(E) = a*E + b*sqrt(E) + c
standard_sigma(x, p) = p[1] * x + p[2] * sqrt(x) + p[3]

# Reactor (inverse beta decay) resolution
function inverse_beta_sigma(x, p)
    x <= 0.0018 && return p[1] / 1000.0
    return p[1] * sqrt((x - 0.0018 + 0.001) * 1000.0) / 1000.0
end

const SIGMA_FUNCTIONS =
    Dict{String, Function}("#standard" => standard_sigma, "#inverse_beta" => inverse_beta_sigma)

function init_bin_data!(sd::SmearData)
    sd.bincenter = zeros(sd.numofbins)
    med = sd.e_min
    for i = 1:sd.numofbins
        med += sd.binsize[i]
        sd.bincenter[i] = med - sd.binsize[i] * 0.5
    end
end

function init_sim_bin_data!(sd::SmearData)
    sd.simbincenter = zeros(sd.simbins)
    med = sd.e_sim_min
    for i = 1:sd.simbins
        med += sd.simbinsize[i]
        sd.simbincenter[i] = med - sd.simbinsize[i] * 0.5
    end
end

lower_bin_boundary(i, sd) = sd.bincenter[i] - sd.binsize[i] * 0.5
upper_bin_boundary(i, sd) = sd.bincenter[i] + sd.binsize[i] * 0.5
lower_sbin_boundary(i, sd) = sd.simbincenter[i] - sd.simbinsize[i] * 0.5
upper_sbin_boundary(i, sd) = sd.simbincenter[i] + sd.simbinsize[i] * 0.5

function default_smear!(sd::SmearData, exp::ExperimentData)
    s = 0
    sd.type == -1 && (exp_error(exp, "Smearing type has not been defined!"); s = -1)
    sd.numofbins <= 0 && (exp_error(exp, "Number of bins must be positive!"); s = -1)
    sd.e_min == -1 && (exp_error(exp, "No minimal energy defined!"); s = -1)
    sd.e_min <= 0 && (exp_error(exp, "Minimal energy must be positive!"); s = -1)
    sd.e_max == -1 && (exp_error(exp, "No maximal energy defined!"); s = -1)
    sd.e_max <= sd.e_min &&
        (glb_error("Maximal energy must be larger than the minimal energy!"); s = -1)
    if sd.e_sim_min != -1 || sd.e_sim_max != -1
        sd.e_sim_min <= 0 && (exp_error(exp, "Minimal sim energy must be positive!"); s = -1)
        sd.e_sim_max <= sd.e_sim_min && (
            exp_error(exp, "Maximal sim energy must be larger than the minimal sim energy!");
            s = -1
        )
    end
    sd.num_of_params == -1 && (sd.num_of_params = 3)
    sd.num_of_params == 0 && (sd.sigma = nothing)
    if sd.num_of_params > 0 && sd.sigma === nothing
        glb_error("No parameters for the energy resolution function specified!")
        s = -1
    end
    sd.sig_f === nothing && (sd.sig_f = standard_sigma)
    if sd.binsize === nothing && s == 0
        sd.binsize =
            exp.binsize === nothing ? fill((sd.e_max - sd.e_min) / sd.numofbins, sd.numofbins) :
            copy(exp.binsize)
    end
    s == 0 && init_bin_data!(sd)
    sd.corr_fac = exp.filter_value
    s != 0 && glb_fatal("glb_smear not properly set up!")
    return s
end

# Determine the sampling points for a computed smearing matrix
function setup_smear_matrix_A!(sd::SmearData, exp::ExperimentData)
    sd.e_sim_min == -1 && (sd.e_sim_min = exp.simtresh == -1 ? exp.emin : exp.simtresh)
    sd.e_sim_max == -1 && (sd.e_sim_max = exp.simbeam == -1 ? exp.emax : exp.simbeam)
    sd.simbins == -1 && (sd.simbins = exp.simbins == -1 ? exp.numofbins : exp.simbins)
    if sd.simbinsize === nothing && exp.simbinsize !== nothing
        sd.simbinsize = copy(exp.simbinsize[1:sd.simbins])
    end
    if sd.simbinsize === nothing
        sd.simbinsize = fill((sd.e_sim_max - sd.e_sim_min) / sd.simbins, sd.simbins)
    end
    init_sim_bin_data!(sd)
end

# Bin data for user-defined (@energy) migration matrices
function set_up_smear_data!(sd::SmearData, exp::ExperimentData)
    sd.e_min == -1 && (sd.e_min = exp.emin)
    sd.e_max == -1 && (sd.e_max = exp.emax)
    sd.numofbins == -1 && (sd.numofbins = exp.numofbins)
    if sd.binsize === nothing
        sd.binsize =
            exp.binsize === nothing ? fill((sd.e_max - sd.e_min) / sd.numofbins, sd.numofbins) :
            copy(exp.binsize)
    end
    init_bin_data!(sd)
    setup_smear_matrix_A!(sd, exp)
    sd.sig_f = nothing
    sd.num_of_params = 0
    sd.type = 0
end

function truncated_sigma(x, sd::SmearData)
    res = sd.sig_f(x, sd.sigma)^2 - sd.corr_fac^2
    res <= 0 && glb_fatal("Filter value too large")
    return sqrt(res)
end

function smear_element(i, j, sd::SmearData, truncated::Bool)
    ec = sd.simbincenter[j]
    sig = truncated ? truncated_sigma(ec, sd) : sd.sig_f(ec, sd.sigma)
    p1 = 0.5 * (1 + erf((upper_bin_boundary(i, sd) - ec) / (sqrt(2) * sig)))
    p2 = 0.5 * (1 + erf((lower_bin_boundary(i, sd) - ec) / (sqrt(2) * sig)))
    return p1 - p2
end

"""
    compute_smearing_matrix(sd, exp) -> (rows, low, up)

Gaussian smearing matrix in the compressed GLoBES format: for every
reconstructed bin `i`, `rows[i]` holds the nonzero entries for sampling
points `low[i]:up[i]` (1-based). Entries are rounded to 1e-7 as in GLoBES.
"""
function compute_smearing_matrix(sd::SmearData, exp::ExperimentData)
    if sd.type != GLB_TYPE_A && sd.type != GLB_TYPE_C
        glb_fatal("Smear type out of range")
    end
    setup_smear_matrix_A!(sd, exp)
    truncated = sd.type == GLB_TYPE_C
    rows = Vector{Vector{Float64}}(undef, sd.numofbins)
    low = zeros(Int, sd.numofbins)
    up = zeros(Int, sd.numofbins)
    for i = 1:sd.numofbins
        nonzero = 0
        zero = 0
        prev = false
        row = Float64[]
        for j = 1:sd.simbins
            erg =
                round((1 / SMEAR_TOL) * smear_element(i, j, sd, truncated), RoundNearestTiesAway) *
                SMEAR_TOL
            if erg != 0
                prev = true
                nonzero += 1
                push!(row, erg)
            elseif !prev
                zero += 1
            end
        end
        rows[i] = row
        low[i] = zero + 1
        up[i] = zero + nonzero
    end
    return rows, low, up
end
