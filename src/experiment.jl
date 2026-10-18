# Experiment initialization (glbDefaultExp in glb_multiex.c) and the event
# rate engine (glb_rate_engine.c).

function exp_error(e::ExperimentData, msg)
    glb_error("In experiment $(basename(e.filename)): $msg")
end

function rule_error(e::ExperimentData, rule, msg)
    glb_error("In experiment $(basename(e.filename)), rule $rule: $msg")
end

function setup_density_profile!(e::ExperimentData)
    s = 0
    if e.density_profile_type == -1
        exp_error(e, "No profile type specified")
        s = -1
    elseif e.density_profile_type == 1
        e.psteps = 1
        if e.baseline > 0
            r = glbAverageDensityProfile(e.baseline)
            r === nothing ? (s = -1) : ((e.lengthtab, e.densitytab) = r)
        else
            exp_error(e, "Baseline must be a positive number")
            s = -1
        end
    elseif e.density_profile_type == 2
        if e.baseline > 0
            if e.psteps > 0
                r = glbStaceyProfile(e.baseline, e.psteps)
                r === nothing ? (s = -1) : ((e.lengthtab, e.densitytab) = r)
            else
                exp_error(e, "Densitysteps must be a positive number")
                s = -1
            end
        else
            exp_error(e, "Baseline must be a positive number")
            s = -1
        end
    elseif e.density_profile_type == 3
        ttl = 0.0
        if !isempty(e.lengthtab) && e.psteps > 0
            for i = 1:e.psteps
                ttl += e.lengthtab[i]
            end
        end
        ttl != e.baseline &&
            e.baseline != -1 &&
            glb_warning("Override baseline with sum of lengthtab")
        e.baseline = ttl
    end
    return s
end

function default_flux!(f::Flux)
    s = 0
    f.builtin == -1 && (f.builtin = 0)
    f.time == -1 && (f.time = 1)
    f.target_power == -1 && (f.target_power = 1)
    f.stored_muons == -1 && (f.stored_muons = 1)
    f.norm == -1 && (f.norm = 1)
    if f.builtin == 1 || f.builtin == 2
        f.parent_energy == -1 && (glb_error("No parent energy (parent_energy) defined!"); s = -1)
    end
    if f.builtin == 3 || f.builtin == 4
        f.gamma == -1 && (glb_error("No gamma defined!"); s = -1)
        f.end_point == -1 && (glb_error("No end point defined!"); s = -1)
    end
    if f.builtin == 0 || f.builtin == GLB_OLD_NORM
        f.file_name === nothing && (glb_error("No flux file specified!"); s = -1)
    end
    return s
end

function default_xsec!(x::XSec)
    x.builtin == -1 && (x.builtin = 0)
    if x.builtin == 0 && x.file_name === nothing
        glb_error("No cross section file specified!")
        return -1
    end
    return 0
end

"""
    default_exp!(e)

Check an experiment read from AEDL, fill in defaults, load fluxes and cross
sections, compute smearing matrices and allocate the rate buffers. Returns
0 on success, < 0 if defaults had to be used or errors occurred.
"""
function default_exp!(e::ExperimentData)
    status = 0
    def = 0
    status += setup_density_profile!(e)

    if e.version === nothing
        glb_error("Missing version in AEDL file")
        def = -1
        e.version = GLB_VERSION
    end

    isempty(e.xsecs) && (exp_error(e, "No X-section selected!"); status = -1)
    isempty(e.fluxes) && (exp_error(e, "No flux selected!"); status = -1)
    for f in e.fluxes
        default_flux!(f) == 0 ? glb_init_fluxtable!(f) : (status = -1)
    end
    for x in e.xsecs
        default_xsec!(x) == 0 ? glb_init_xsectable!(x) : (status = -1)
    end

    # Systematics
    for (i, r) in enumerate(e.rules)
        for on_off in (GLB_ON, GLB_OFF)
            name = on_off == GLB_ON ? r.sys_on_string : r.sys_off_string
            if name === nothing
                rule_error(e, i, "No chi^2 function specified")
                status = -1
                continue
            end
            tmp = on_off == GLB_ON ? r.sys_on_errors : r.sys_off_errors
            if on_off == GLB_ON
                r.sys_on_errors = nothing
            else
                r.sys_off_errors = nothing
            end
            dim = glbGetSysDim(name)
            if dim >= 0
                if tmp !== nothing
                    k = 0
                    while k < length(tmp) && tmp[k + 1] > 0.0
                        k += 1
                    end
                    if k != dim
                        rule_error(
                            e,
                            i,
                            "Invalid systematical error list @sys_$(on_off == GLB_ON ? "on" : "off")_errors",
                        )
                        status = -1
                    end
                end
            else
                rule_error(e, i, "Invalid systematics specification")
                status = -1
            end
            if set_chi_function_in_rule!(e, i, on_off, name, tmp) != 0
                rule_error(e, i, "Invalid systematics specification")
                status = -1
            end
        end
    end

    e.baseline == -1 && (exp_error(e, "No baseline specified!"); status = -1)
    e.emin == -1 && (exp_error(e, "No emin specified!"); status = -1)
    e.emax == -1 && (exp_error(e, "No emax specified!"); status = -1)
    e.emin > e.emax && (exp_error(e, "emin must be less than emax!"); status = -1)
    e.numofbins <= 0 && (exp_error(e, "numofbins is not set!"); status = -1)
    e.targetmass == -1 && (exp_error(e, "No targetmass specified!"); status = -1)
    isempty(e.channels) && (exp_error(e, "No channels specified!"); status = -1)
    isempty(e.rules) && (exp_error(e, "No rules specified!"); status = -1)
    for (i, r) in enumerate(e.rules)
        isempty(r.signal_channels) && (rule_error(e, i, "No signal specified!"); status = -1)
    end
    for (i, c) in enumerate(e.channels)
        if c.flux < 1 ||
           c.flux > length(e.fluxes) ||
           c.xsec < 1 ||
           c.xsec > length(e.xsecs) ||
           c.energy < 1 ||
           c.energy > length(e.smear_data)
            exp_error(e, "Channel $i is not properly defined!")
            status = -1
        end
    end

    e.filter_state == -1 && (e.filter_state = 1; def = -1)
    e.filter_value == -1 && (e.filter_value = 0; def = -1)

    isempty(e.smear_data) && (exp_error(e, "No smearing data specified!"); status = -1)
    status != 0 && glb_fatal("Incompletely defined experiment!")

    for i in eachindex(e.smear_data)
        sd = e.smear_data[i]
        if e.user_smear_given[i]
            set_up_smear_data!(sd, e)
            default_smear!(sd, e)
        else
            sd.numofbins = e.numofbins
            sd.e_min = e.emin
            sd.e_max = e.emax
            default_smear!(sd, e)
            e.smear[i], e.lowrange[i], e.uprange[i] = compute_smearing_matrix(sd, e)
            e.simtresh = sd.e_sim_min
            e.simbeam = sd.e_sim_max
            e.simbins = sd.simbins
        end
        if e.simtresh == -1 || e.simbeam == -1 || e.simbins == -1
            # User-defined matrices without explicit sampling settings
            e.simtresh = sd.e_sim_min
            e.simbeam = sd.e_sim_max
            e.simbins = sd.simbins
        end
        e.simbins < e.numofbins && (exp_error(e, "Less sampling points than bins"); status = -1)
    end

    # Buffers
    nb = e.numofbins
    ns = e.simbins
    nr = length(e.rules)
    nc = length(e.channels)
    e.rates0 = [zeros(nb) for _ = 1:nr]
    e.rates1 = [zeros(nb) for _ = 1:nr]
    e.rates1BG = [zeros(nb) for _ = 1:nr]
    e.SignalRates = [zeros(nb) for _ = 1:nr]
    e.BackgroundRates = [zeros(nb) for _ = 1:nr]
    e.energy_tab = zeros(nb)
    e.chrb_0 = [zeros(ns) for _ = 1:nc]
    e.chrb_1 = [zeros(ns) for _ = 1:nc]
    e.chra_0 = [zeros(nb) for _ = 1:nc]
    e.chra_1 = [zeros(nb) for _ = 1:nc]
    e.psteps <= 0 && (exp_error(e, "Too few density steps defined!"); status = -1)
    e.densitybuffer = copy(e.densitytab)
    length(e.lengthtab) != e.psteps && (exp_error(e, "lengthtab has wrong length"); status = -1)
    length(e.densitytab) != e.psteps && (exp_error(e, "densitytab has wrong length"); status = -1)

    for i = 1:nc
        fields = (
            (:user_pre_smearing_channel, ns, 1.0),
            (:user_pre_smearing_background, ns, 0.0),
            (:user_post_smearing_channel, nb, 1.0),
            (:user_post_smearing_background, nb, 0.0),
        )
        for (f, len, val) in fields
            v = getfield(e, f)
            if v[i] === nothing
                v[i] = fill(val, len)
                def = -1
            elseif length(v[i]) != len
                exp_error(e, "$(f) of channel $i has $(length(v[i])) instead of $len elements")
                status = -1
            end
        end
    end

    sd0 = e.smear_data[1]
    for r in e.rules
        if r.energy_window[1] == -1 || r.energy_window[1] < e.emin
            r.energy_window[1] = e.emin
            def = -1
        end
        if r.energy_window[2] == -1 || r.energy_window[2] > e.emax
            r.energy_window[2] = e.emax
            def = -1
        end
        if r.energy_window[1] > r.energy_window[2]
            r.energy_window[1], r.energy_window[2] = r.energy_window[2], r.energy_window[1]
        end
        # Bin range of the energy window (as in GLoBES' glb_window_function)
        k = 1
        while k <= nb && sd0.bincenter[k] <= r.energy_window[1]
            k += 1
        end
        r.energy_window_bins[1] = k
        while k <= nb && sd0.bincenter[k] < r.energy_window[2]
            k += 1
        end
        r.energy_window_bins[2] = k - 1
    end

    status != 0 && glb_fatal("Incompletely defined experiment!")
    def != 0 && glb_warning("Incompletely defined experiment! Defaults are going to be used!")
    return def + status
end

# ---------------------------------------------------------------------------
#  Rate engine
# ---------------------------------------------------------------------------

function calc_all_probs!(e::ExperimentData, en)
    st = state()
    fs = e.filter_state == GLB_ON ? e.filter_value : -1.0
    s1 = st.prob_func(st.P, +1, en, e.psteps, e.lengthtab, e.densitybuffer, fs, st.prob_user_data)
    s2 =
        st.prob_func(st.Pbar, -1, en, e.psteps, e.lengthtab, e.densitybuffer, fs, st.prob_user_data)
    (s1 != 0 || s2 != 0) && glb_error("Calculation of oscillation probabilities failed.")
end

function calc_channel(e::ExperimentData, c::Channel, en)
    st = state()
    fin = c.final > 9 ? c.final - 10 : c.final
    ini = c.initial > 9 ? c.initial - 10 : c.initial
    flux = e.fluxes[c.flux]
    xsec = e.xsecs[c.xsec]
    if is_nosc(c)
        fin != ini && return 0.0
        return glb_xsec_calc(en, fin, c.cp, xsec) *
               glb_flux_calc(en, e.baseline, ini, c.cp, flux) *
               e.targetmass
    end
    P = c.cp == 1 ? st.P : st.Pbar
    return glb_xsec_calc(en, fin, c.cp, xsec) *
           glb_flux_calc(en, e.baseline, ini, c.cp, flux) *
           P[ini, fin] *
           e.targetmass
end

function apply_smearing!(out, e::ExperimentData, c_idx, pre)
    c = e.channels[c_idx]
    s = c.energy
    sm = e.smear[s]
    low = e.lowrange[s]
    up = e.uprange[s]
    eff = e.user_post_smearing_channel[c_idx]
    bg = e.user_post_smearing_background[c_idx]
    for j = 1:e.numofbins
        acc = 0.0
        row = sm[j]
        lo = low[j]
        for k = lo:up[j]
            acc += row[k - lo + 1] * pre[k]
        end
        out[j] = acc * eff[j] + bg[j]
    end
end

"""
    set_rates!(e)

Compute the "true" event rates (as `glbSetRates` does) for the current
oscillation parameters.
"""
function set_rates!(e::ExperimentData)
    sd0 = e.smear_data[1]
    for j = 1:e.numofbins
        e.energy_tab[j] = sd0.bincenter[j]
    end
    for j = 1:e.simbins
        calc_all_probs!(e, sd0.simbincenter[j])
        for (i, c) in enumerate(e.channels)
            sd = e.smear_data[c.energy]
            e.chrb_0[i][j] =
                calc_channel(e, c, sd.simbincenter[j]) *
                e.user_pre_smearing_channel[i][j] *
                sd.simbinsize[j] + e.user_pre_smearing_background[i][j]
            is_nosc(c) && (e.chrb_1[i][j] = e.chrb_0[i][j])
        end
    end
    for (i, c) in enumerate(e.channels)
        apply_smearing!(e.chra_0[i], e, i, e.chrb_0[i])
        is_nosc(c) && (e.chra_1[i] .= e.chra_0[i])
    end
    for (i, r) in enumerate(e.rules)
        bg = e.BackgroundRates[i]
        sig = e.SignalRates[i]
        for j = 1:e.numofbins
            b = 0.0
            for k in eachindex(r.bg_channels)
                b += r.bg_coeff[k] * e.chra_0[r.bg_channels[k]][j]
            end
            bg[j] = b
            s = 0.0
            for k in eachindex(r.signal_channels)
                s += r.signal_coeff[k] * e.chra_0[r.signal_channels[k]][j]
            end
            sig[j] = s
        end
        fill!(e.rates0[i], 0.0)
        for j = r.energy_window_bins[1]:r.energy_window_bins[2]
            e.rates0[i][j] = sig[j] + bg[j] * r.bg_centers[1]
        end
    end
end

"""
    set_new_rates!(e)

Compute the fitted event rates for the current oscillation parameters
(channels flagged NOSC_ are not recomputed).
"""
function set_new_rates!(e::ExperimentData)
    sd0 = e.smear_data[1]
    for j = 1:e.simbins
        calc_all_probs!(e, sd0.simbincenter[j])
        for (i, c) in enumerate(e.channels)
            is_nosc(c) && continue
            sd = e.smear_data[c.energy]
            e.chrb_1[i][j] =
                calc_channel(e, c, sd.simbincenter[j]) *
                e.user_pre_smearing_channel[i][j] *
                sd.simbinsize[j] + e.user_pre_smearing_background[i][j]
        end
    end
    for (i, c) in enumerate(e.channels)
        is_nosc(c) && continue
        apply_smearing!(e.chra_1[i], e, i, e.chrb_1[i])
    end
    for (i, r) in enumerate(e.rules)
        rb = e.rates1BG[i]
        rs = e.rates1[i]
        for j = 1:e.numofbins
            b = 0.0
            for k in eachindex(r.bg_channels)
                b += r.bg_coeff[k] * e.chra_1[r.bg_channels[k]][j]
            end
            rb[j] = b
            s = 0.0
            for k in eachindex(r.signal_channels)
                s += r.signal_coeff[k] * e.chra_1[r.signal_channels[k]][j]
            end
            rs[j] = s
        end
    end
end

set_profile_scaling!(e::ExperimentData, scale) = (e.densitybuffer .= scale .* e.densitytab)
