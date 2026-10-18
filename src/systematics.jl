# Systematics and chi^2 functions (port of glb_sys.c)
#
# A chi^2 function has the signature
#     f(exp::Int, rule::Int, n_params::Int, x::AbstractVector, errors::AbstractVector, user_data) -> Float64
# with 1-based experiment and rule numbers.

const GLB_2011_COMPATIBLE = ("chiSpectrumTilt", "chiNoSysSpectrum", "chiSpectrumOnly",
                             "chiTotalRatesTilt", "chiNoSysTotalRates", "chiSpectrumCalib")

is_2011_compatible(name) = name !== nothing && name in GLB_2011_COMPATIBLE

const FIT_RATE_FLOOR = 1e-10   # events; below this the Poisson term is extrapolated

"""
    glb_likelihood(true_rate, fit_rate)

Poissonian χ² contribution `2(fit - true + true*log(true/fit))`.

C GLoBES returns 1e100 whenever `fit_rate <= 0`, even for empty bins
(true = fit = 0), which ruins the χ² of any rule with such a bin and makes
χ² jump for minimizers. Here an empty bin contributes 0, a negative fit rate
with no observed events contributes `2|fit|`, and with observed events the
function is continued below `FIT_RATE_FLOOR` by its second-order Taylor
expansion: a steep but finite, smooth and convex wall.
"""
@inline function glb_likelihood(true_rate, fit_rate)
    if true_rate <= 0.0
        return 2.0 * abs(fit_rate)
    elseif fit_rate >= FIT_RATE_FLOOR
        return 2.0 * (fit_rate - true_rate + true_rate * log(true_rate / fit_rate))
    end
    f0 = FIT_RATE_FLOOR
    d = fit_rate - f0
    value = f0 - true_rate + true_rate * log(true_rate / f0)
    slope = 1.0 - true_rate / f0
    curvature = true_rate / f0^2
    return 2.0 * (value + slope * d + 0.5 * curvature * d^2)
end

"""
    glb_prior(x, center, sigma)

Gaussian pull term `((x - center)/sigma)^2`.
"""
@inline glb_prior(x, center, sigma) = ((x - center) / sigma)^2

function find_chi_function(name)
    for s in state().sys_list
        s.name == name && return s
    end
    glb_error("glbFindChiFunctionByName: Unknown chi^2 function specified")
    return nothing
end

function convert_errordim(errordim)
    errordim == 0 && return "chiSpectrumTilt"
    errordim == 2 && return "chiNoSysSpectrum"
    errordim == 4 && return "chiTotalRatesTilt"
    errordim == 7 && return "chiSpectrumOnly"
    errordim == 8 && return "chiNoSysTotalRates"
    errordim == 9 && return "chiSpectrumCalib"
    glb_error("glbConvertErrorDim: Invalid errordim")
    return nothing
end

# ---------------------------------------------------------------------------
#  Built-in chi^2 functions
# ---------------------------------------------------------------------------

function glbChiSpectrumTilt(exp, rule, n_params, x, errors, user_data)
    e = state().experiments[exp]; r = e.rules[rule]
    true_rates = e.rates0[rule]; sig = e.rates1[rule]; bg = e.rates1BG[rule]
    bin_centers = e.smear_data[1].bincenter
    emin, emax = e.emin, e.emax
    ecenter = 0.5 * (emax + emin)
    bg_norm_center, bg_tilt_center = r.bg_centers[1], r.bg_centers[2]
    signal_norm = 1.0 + x[1]
    signal_tilt = x[2] / (emax - emin)
    bg_norm = bg_norm_center * (1.0 + x[3])
    bg_tilt = x[4] / (emax - emin)
    chi2 = 0.0
    for i in r.energy_window_bins[1]:r.energy_window_bins[2]
        fit_rate = signal_norm * sig[i] + signal_tilt * (bin_centers[i] - ecenter) * sig[i] +
                   bg_norm * bg[i] + bg_tilt * (bin_centers[i] - ecenter) * bg[i]
        chi2 += glb_likelihood(true_rates[i], fit_rate)
    end
    chi2 += glb_prior(x[1], 0.0, errors[1]) + glb_prior(x[2], 0.0, errors[2]) +
            glb_prior(bg_norm, bg_norm_center, errors[3]) + glb_prior(bg_tilt, bg_tilt_center, errors[4])
    return chi2
end

function glbChiNoSysSpectrum(exp, rule, n_params, x, errors, user_data)
    e = state().experiments[exp]; r = e.rules[rule]
    true_rates = e.rates0[rule]; sig = e.rates1[rule]; bg = e.rates1BG[rule]
    bgc = r.bg_centers[1]
    chi2 = 0.0
    for i in r.energy_window_bins[1]:r.energy_window_bins[2]
        chi2 += glb_likelihood(true_rates[i], sig[i] + bgc * bg[i])
    end
    return chi2
end

function glbChiSpectrumOnly(exp, rule, n_params, x, errors, user_data)
    e = state().experiments[exp]; r = e.rules[rule]
    true_rates = e.rates0[rule]; sig = e.rates1[rule]; bg = e.rates1BG[rule]
    signal_norm = 1.0 + x[1]
    bgc = r.bg_centers[1]
    chi2 = 0.0
    for i in r.energy_window_bins[1]:r.energy_window_bins[2]
        chi2 += glb_likelihood(true_rates[i], signal_norm * sig[i] + bgc * bg[i])
    end
    return chi2
end

function glbChiTotalRatesTilt(exp, rule, n_params, x, errors, user_data)
    e = state().experiments[exp]; r = e.rules[rule]
    true_rates = e.rates0[rule]; sig = e.rates1[rule]; bg = e.rates1BG[rule]
    bin_centers = e.smear_data[1].bincenter
    emin, emax = e.emin, e.emax
    ecenter = 0.5 * (emax + emin)
    bg_norm_center, bg_tilt_center = r.bg_centers[1], r.bg_centers[2]
    signal_norm = 1.0 + x[1]
    signal_tilt = x[2] / (emax - emin)
    bg_norm = bg_norm_center * (1.0 + x[3])
    bg_tilt = x[4] / (emax - emin)
    total_fit = 0.0; total_true = 0.0
    for i in r.energy_window_bins[1]:r.energy_window_bins[2]
        total_fit += signal_norm * sig[i] + signal_tilt * (bin_centers[i] - ecenter) * sig[i] +
                     bg_norm * bg[i] + bg_tilt * (bin_centers[i] - ecenter) * bg[i]
        total_true += true_rates[i]
    end
    chi2 = glb_likelihood(total_true, total_fit)
    chi2 += glb_prior(x[1], 0.0, errors[1]) + glb_prior(x[2], 0.0, errors[2]) +
            glb_prior(bg_norm, bg_norm_center, errors[3]) + glb_prior(bg_tilt, bg_tilt_center, errors[4])
    return chi2
end

function glbChiNoSysTotalRates(exp, rule, n_params, x, errors, user_data)
    e = state().experiments[exp]; r = e.rules[rule]
    true_rates = e.rates0[rule]; sig = e.rates1[rule]; bg = e.rates1BG[rule]
    bgc = r.bg_centers[1]
    total_fit = 0.0; total_true = 0.0
    for i in r.energy_window_bins[1]:r.energy_window_bins[2]
        total_fit += sig[i] + bgc * bg[i]
        total_true += true_rates[i]
    end
    return glb_likelihood(total_true, total_fit)
end

function glbChiSpectrumCalib(exp, rule, n_params, x, errors, user_data)
    e = state().experiments[exp]; r = e.rules[rule]
    n_bins = e.numofbins
    true_rates = e.rates0[rule]
    bg_norm_center, bg_tilt_center = r.bg_centers[1], r.bg_centers[2]
    signal_norm = 1.0 + x[1]
    bg_norm = bg_norm_center * (1.0 + x[3])
    bg_tilt = bg_tilt_center + x[4]
    centers, widths = e.smear_data[1].bincenter, e.smear_data[1].binsize
    sig = shift_energy_scale(x[2], view(e.rates1[rule], 1:n_bins), centers, widths)
    bg = shift_energy_scale(bg_tilt, view(e.rates1BG[rule], 1:n_bins), centers, widths)
    chi2 = 0.0
    for i in r.energy_window_bins[1]:r.energy_window_bins[2]
        chi2 += glb_likelihood(true_rates[i], signal_norm * sig[i] + bg_norm * bg[i])
    end
    chi2 += glb_prior(x[1], 0.0, errors[1]) + glb_prior(x[2], 0.0, errors[2]) +
            glb_prior(bg_norm, bg_norm_center, errors[3]) + glb_prior(bg_tilt, bg_tilt_center, errors[4])
    return chi2
end

glbChiZero(exp, rule, n_params, x, errors, user_data) = 0.0

"""
    shift_energy_scale(b, rates, centers, widths) -> rates_out

Apply an energy calibration error of relative size `b` to a binned spectrum
with the given bin centers and widths: bin `i` of the result holds the events
between the scaled bin edges, `F((1+b) hi[i]) - F((1+b) lo[i])`, where `F` is
the cumulative number of events. `F` is known exactly at the bin edges and is
interpolated in between with a monotone cubic (Fritsch-Carlson); as in C
GLoBES, there are no events outside the binned range.

This reproduces `rates` at `b = 0`, never gives negative rates, conserves the
number of events, handles unequal bins and is smooth in `b`, except in the
first and last bin if the spectrum does not vanish at the ends of the binned
range (shifting then moves events out of the range, but none in). C GLoBES
instead interpolates the rates linearly at the shifted bin centers, which
makes χ² kinked at `b = 0` in every bin, doubles the last bin when shifting
past the upper end, and assumes equal bins.
"""
function shift_energy_scale(b, rates, centers, widths)
    n = length(rates)
    edges = [centers[1] - widths[1] / 2; centers .+ widths ./ 2]
    cumulative = [0.0; cumsum(rates)]
    slopes = monotone_slopes(edges, cumulative)
    F(E) = E <= edges[1] ? 0.0 : E >= edges[end] ? cumulative[end] :
           hermite_eval(edges, cumulative, slopes, E)
    return [F((1 + b) * edges[i+1]) - F((1 + b) * edges[i]) for i in 1:n]
end

# Node slopes of the monotone piecewise cubic through (x, y) (Fritsch & Carlson
# 1980, with the weighted harmonic mean of Fritsch & Butland 1984, and the
# shape-preserving one-sided end slopes used by MATLAB's pchip).
function monotone_slopes(x, y)
    m = length(x)
    h = diff(x)
    s = diff(y) ./ h
    slopes = zeros(m)
    for k in 2:m-1
        if s[k-1] * s[k] > 0
            slopes[k] = 3(h[k-1] + h[k]) / ((2h[k] + h[k-1]) / s[k-1] + (h[k] + 2h[k-1]) / s[k])
        end
    end
    if m == 2
        slopes .= s[1]
    elseif m > 2
        slopes[1] = end_slope(h[1], h[2], s[1], s[2])
        slopes[m] = end_slope(h[m-1], h[m-2], s[m-1], s[m-2])
    end
    return slopes
end

function end_slope(h1, h2, s1, s2)
    d = ((2h1 + h2) * s1 - h1 * s2) / (h1 + h2)
    sign(d) != sign(s1) && return 0.0
    sign(s1) != sign(s2) && abs(d) > abs(3s1) && return 3s1
    return d
end

# Evaluate the cubic Hermite interpolant through (x, y) with the given node
# slopes at x[1] <= E <= x[end]
function hermite_eval(x, y, slopes, E)
    k = min(searchsortedlast(x, E), length(x) - 1)
    h = x[k+1] - x[k]
    t = (E - x[k]) / h
    t2 = t * t; t3 = t2 * t
    return (2t3 - 3t2 + 1) * y[k] + (t3 - 2t2 + t) * h * slopes[k] +
           (-2t3 + 3t2) * y[k+1] + (t3 - t2) * h * slopes[k+1]
end

"""
    glbShiftEnergyScale(b, rates_in, n_bins, emin, emax) -> rates_out

Apply an energy calibration error of relative size `b` to a spectrum of
`n_bins` equal bins between `emin` and `emax`. See [`shift_energy_scale`](@ref)
for the method and how it differs from C GLoBES.
"""
function glbShiftEnergyScale(b, rates_in, n_bins, emin, emax)
    w = (emax - emin) / n_bins
    centers = [emin + (i - 0.5) * w for i in 1:n_bins]
    return shift_energy_scale(b, view(rates_in, 1:n_bins), centers, fill(w, n_bins))
end

# ---------------------------------------------------------------------------
#  API
# ---------------------------------------------------------------------------

"""
    glbDefineChiFunction(chi_func, dim, name, user_data=nothing)

Register a χ² function under `name` for use in AEDL files (`@sys_on_function`)
or with `glbSetChiFunction`. `chi_func(exp, rule, n_params, x, errors, user_data)`
gets 1-based experiment/rule numbers and a vector `x` of `dim` nuisance parameters.
"""
function glbDefineChiFunction(chi_func::Function, dim::Integer, name::AbstractString, user_data=nothing)
    dim < 0 && (glb_error("glbDefineChiFunction: Dimension of parameter space must be >= 0"); return -1)
    pushfirst!(state().sys_list, ChiFunction(chi_func, dim, String(name), user_data))
    return 0
end

function set_chi_function_in_rule!(e::ExperimentData, rule, on_off, sys_id, errors)
    r = e.rules[rule]
    sys = find_chi_function(sys_id)
    sys === nothing && return -1
    if errors === nothing && sys.dim > 0
        if is_2011_compatible(sys_id)
            if r.signal_errors[1] > 0 && r.signal_errors[2] > 0 && r.bg_errors[1] > 0 && r.bg_errors[2] > 0
                errors = [r.signal_errors[1], r.signal_errors[2], r.bg_errors[1], r.bg_errors[2]]
            else
                glb_error("glbSetChiFunctionInRule: Incomplete systematics definition")
                return -1
            end
        else
            glb_error("glbSetChiFunctionInRule: Invalid error array")
            return -1
        end
    end
    if sys.dim > 0 && length(errors) < sys.dim
        glb_error("glbSetChiFunctionInRule: Error list too short")
        return -1
    end
    errs = sys.dim > 0 ? Float64.(errors[1:sys.dim]) : nothing
    starts = sys.dim > 0 ? zeros(sys.dim) : nothing
    if on_off == GLB_ON
        r.sys_on = sys; r.sys_on_string = sys.name
        r.sys_on_errors = errs; r.sys_on_startvals = starts
    elseif on_off == GLB_OFF
        r.sys_off = sys; r.sys_off_string = sys.name
        r.sys_off_errors = errs; r.sys_off_startvals = starts
    else
        glb_error("glbSetChiFunctionInRule: on_off must be GLB_ON or GLB_OFF")
        return -1
    end
    if is_2011_compatible(sys.name)
        sys.dim > 0 && (r.signal_errors[1] = errors[1]; r.signal_startvals[1] = 0.0)
        sys.dim > 1 && (r.signal_errors[2] = errors[2]; r.signal_startvals[2] = 0.0)
        sys.dim > 2 && (r.bg_errors[1] = errors[3]; r.bg_startvals[1] = 0.0)
        sys.dim > 3 && (r.bg_errors[2] = errors[4]; r.bg_startvals[2] = 0.0)
    end
    return 0
end

# Iterate over (experiment, rule) pairs, expanding GLB_ALL
function foreach_exp_rule(f, exp, rule)
    exps = exp == GLB_ALL ? (1:length(state().experiments)) : (exp:exp)
    for i in exps
        (1 <= i <= length(state().experiments)) || (glb_error("Invalid experiment number"); return -1)
        rules = rule == GLB_ALL ? (1:length(state().experiments[i].rules)) : (rule:rule)
        for j in rules
            (1 <= j <= length(state().experiments[i].rules)) || (glb_error("Invalid rule number"); return -1)
            f(i, j) != 0 && return -1
        end
    end
    return 0
end

"""
    glbSetChiFunction(exp, rule, on_off, name, errors)

Select the χ² function `name` for systematics state `on_off` (`GLB_ON`/`GLB_OFF`).
`exp` and `rule` may be `GLB_ALL`.
"""
glbSetChiFunction(exp, rule, on_off, sys_id, errors=nothing) =
    foreach_exp_rule((i, j) -> set_chi_function_in_rule!(state().experiments[i], j, on_off, sys_id, errors), exp, rule)

"""
    glbGetChiFunction(exp, rule, on_off) -> String
"""
function glbGetChiFunction(exp, rule, on_off)
    r = getrule(exp, rule); r === nothing && return nothing
    return on_off == GLB_ON ? r.sys_on_string : r.sys_off_string
end

glbGetChiFunctionPtr(name) = (s = find_chi_function(name); s === nothing ? nothing : s.func)

function glbGetChiFunctionPtrInExperiment(exp, rule, on_off)
    r = getrule(exp, rule); r === nothing && return nothing
    return (on_off == GLB_ON ? r.sys_on : r.sys_off).func
end

glbSetErrorDim(exp, rule, on_off, errordim) =
    (name = convert_errordim(errordim); name === nothing ? -1 : glbSetChiFunction(exp, rule, on_off, name, nothing))

function glbGetErrorDim(exp, rule, on_off)
    name = glbGetChiFunction(exp, rule, on_off)
    d = Dict("chiSpectrumTilt" => 0, "chiNoSysSpectrum" => 2, "chiTotalRatesTilt" => 4,
             "chiSpectrumOnly" => 7, "chiNoSysTotalRates" => 8, "chiSpectrumCalib" => 9)
    return get(d, name, -1)
end

"""
    glbGetSysDim(name) -> Int

Number of nuisance parameters of the registered χ² function `name`.
"""
glbGetSysDim(name) = (s = find_chi_function(name); s === nothing ? -1 : s.dim)

function glbGetSysDimInExperiment(exp, rule, on_off)
    r = getrule(exp, rule); r === nothing && return -1
    sys = on_off == GLB_ON ? r.sys_on : r.sys_off
    sys === nothing && return -1
    return sys.dim
end

"""
    glbSwitchSystematics(exp, rule, on_off)

Switch between the `GLB_ON` and `GLB_OFF` systematics definitions.
"""
function glbSwitchSystematics(exp, rule, on_off)
    (on_off == GLB_ON || on_off == GLB_OFF) || (glb_error("glbSwitchSystematics: Invalid value for on_off"); return -1)
    foreach_exp_rule((i, j) -> (state().experiments[i].rules[j].sys_on_off = on_off; 0), exp, rule)
end

glbGetSysOnOffState(exp, rule) = (r = getrule(exp, rule); r === nothing ? -1 : r.sys_on_off)

active_sys(r::Rule) = r.sys_on_off == GLB_ON ? r.sys_on : r.sys_off
active_errors(r::Rule) = r.sys_on_off == GLB_ON ? r.sys_on_errors : r.sys_off_errors
active_startvals(r::Rule) = r.sys_on_off == GLB_ON ? r.sys_on_startvals : r.sys_off_startvals

function sync_2011!(r::Rule, on_off, idx_offset, vals, field)
    name = on_off == GLB_ON ? r.sys_on_string : r.sys_off_string
    is_2011_compatible(name) || return
    sys = on_off == GLB_ON ? r.sys_on : r.sys_off
    target = field == :errors ? (on_off == GLB_ON ? r.sys_on_errors : r.sys_off_errors) :
                                (on_off == GLB_ON ? r.sys_on_startvals : r.sys_off_startvals)
    for k in 1:2
        sys.dim > idx_offset + k - 1 && (target[idx_offset+k] = vals[k])
    end
end

"""
    glbSetSignalErrors(exp, rule, norm, tilt)

Set signal normalization and tilt errors (for the GLoBES 2.0.11-style χ² functions).
"""
function glbSetSignalErrors(exp, rule, norm, tilt)
    (norm > 0 && tilt > 0) || (glb_error("glbSetSignalErrors: Errors must be > 0"); return -1)
    foreach_exp_rule(exp, rule) do i, j
        r = state().experiments[i].rules[j]
        r.signal_errors .= (norm, tilt)
        sync_2011!(r, GLB_ON, 0, r.signal_errors, :errors)
        sync_2011!(r, GLB_OFF, 0, r.signal_errors, :errors)
        0
    end
end
glbGetSignalErrors(exp, rule) = (r = getrule(exp, rule); r === nothing ? nothing : (r.signal_errors[1], r.signal_errors[2]))

"""
    glbSetBGErrors(exp, rule, norm, tilt)

Set background normalization and tilt errors.
"""
function glbSetBGErrors(exp, rule, norm, tilt)
    (norm > 0 && tilt > 0) || (glb_error("glbSetBGErrors: Errors must be > 0"); return -1)
    foreach_exp_rule(exp, rule) do i, j
        r = state().experiments[i].rules[j]
        r.bg_errors .= (norm, tilt)
        sync_2011!(r, GLB_ON, 2, r.bg_errors, :errors)
        sync_2011!(r, GLB_OFF, 2, r.bg_errors, :errors)
        0
    end
end
glbGetBGErrors(exp, rule) = (r = getrule(exp, rule); r === nothing ? nothing : (r.bg_errors[1], r.bg_errors[2]))

"""
    glbSetBGCenters(exp, rule, norm, tilt)

Central values of background normalization and tilt (`@backgroundcenter`).
"""
glbSetBGCenters(exp, rule, norm, tilt) =
    foreach_exp_rule((i, j) -> (state().experiments[i].rules[j].bg_centers .= (norm, tilt); 0), exp, rule)
glbGetBGCenters(exp, rule) = (r = getrule(exp, rule); r === nothing ? nothing : (r.bg_centers[1], r.bg_centers[2]))

function glbSetSignalStartingValues(exp, rule, norm, tilt)
    foreach_exp_rule(exp, rule) do i, j
        r = state().experiments[i].rules[j]
        r.signal_startvals .= (norm, tilt)
        sync_2011!(r, GLB_ON, 0, r.signal_startvals, :startvals)
        sync_2011!(r, GLB_OFF, 0, r.signal_startvals, :startvals)
        0
    end
end
glbGetSignalStartingValues(exp, rule) = (r = getrule(exp, rule); r === nothing ? nothing : Tuple(r.signal_startvals))

function glbSetBGStartingValues(exp, rule, norm, tilt)
    foreach_exp_rule(exp, rule) do i, j
        r = state().experiments[i].rules[j]
        r.bg_startvals .= (norm, tilt)
        sync_2011!(r, GLB_ON, 2, r.bg_startvals, :startvals)
        sync_2011!(r, GLB_OFF, 2, r.bg_startvals, :startvals)
        0
    end
end
glbGetBGStartingValues(exp, rule) = (r = getrule(exp, rule); r === nothing ? nothing : Tuple(r.bg_startvals))

function set_list_in_rule!(r::Rule, on_off, list, which)
    sys = on_off == GLB_ON ? r.sys_on : r.sys_off
    sys === nothing && return -1
    dim = sys.dim
    length(list) < dim && (glb_error("List too short"); return -1)
    dim == 0 && return 0
    target = which == :errors ? (on_off == GLB_ON ? r.sys_on_errors : r.sys_off_errors) :
                                (on_off == GLB_ON ? r.sys_on_startvals : r.sys_off_startvals)
    target[1:dim] .= list[1:dim]
    if is_2011_compatible(sys.name)
        a = which == :errors ? (r.signal_errors, r.bg_errors) : (r.signal_startvals, r.bg_startvals)
        dim > 0 && (a[1][1] = list[1]); dim > 1 && (a[1][2] = list[2])
        dim > 2 && (a[2][1] = list[3]); dim > 3 && (a[2][2] = list[4])
    end
    return 0
end

"""
    glbSetSysErrorsList(exp, rule, on_off, list)

Set the systematical errors passed to the χ² function of a rule.
"""
glbSetSysErrorsList(exp, rule, on_off, list) =
    foreach_exp_rule((i, j) -> set_list_in_rule!(state().experiments[i].rules[j], on_off, list, :errors), exp, rule)

glbGetSysErrorsListPtr(exp, rule, on_off) =
    (r = getrule(exp, rule); r === nothing ? nothing : (on_off == GLB_ON ? r.sys_on_errors : r.sys_off_errors))

"""
    glbSetSysStartingValuesList(exp, rule, on_off, list)

Starting values for the nuisance-parameter minimization of a rule.
"""
glbSetSysStartingValuesList(exp, rule, on_off, list) =
    foreach_exp_rule((i, j) -> set_list_in_rule!(state().experiments[i].rules[j], on_off, list, :startvals), exp, rule)

glbGetSysStartingValuesListPtr(exp, rule, on_off) =
    (r = getrule(exp, rule); r === nothing ? nothing : (on_off == GLB_ON ? r.sys_on_startvals : r.sys_off_startvals))

function register_builtin_chi_functions!()
    empty!(state().sys_list)
    glbDefineChiFunction(glbChiSpectrumTilt, 4, "chiSpectrumTilt")
    glbDefineChiFunction(glbChiNoSysSpectrum, 0, "chiNoSysSpectrum")
    glbDefineChiFunction(glbChiSpectrumOnly, 1, "chiSpectrumOnly")
    glbDefineChiFunction(glbChiTotalRatesTilt, 4, "chiTotalRatesTilt")
    glbDefineChiFunction(glbChiNoSysTotalRates, 0, "chiNoSysTotalRates")
    glbDefineChiFunction(glbChiSpectrumCalib, 4, "chiSpectrumCalib")
    glbDefineChiFunction(glbChiZero, 0, "chiZero")
end
