# Minimizers (glb_min_sup.c) and chi^2 projections (glb_minimize.c)

const TOLSYS = 1.0e-8
const TOLOSC = 1.0e-5

# ---------------------------------------------------------------------------
#  Powell's method with Brent line minimization (Numerical Recipes), as
#  used by GLoBES. Vectors are 1-based like the NR originals.
# ---------------------------------------------------------------------------

const BRENT_TOL = 2.0e-4
const CGOLD = 0.3819660
const ZEPS = 1.0e-10
const GOLD = 1.618034
const GLIMIT = 100.0
const TINY = 1.0e-20

nrsign(a, b) = b > 0.0 ? abs(a) : -abs(a)

function bracket(f, ax, bx)
    fa = f(ax)
    fb = f(bx)
    if fb > fa
        ax, bx = bx, ax
        fa, fb = fb, fa
    end
    cx = bx + GOLD * (bx - ax)
    fc = f(cx)
    while fb > fc
        r = (bx - ax) * (fb - fc)
        q = (bx - cx) * (fb - fa)
        u = bx - ((bx - cx) * q - (bx - ax) * r) / (2.0 * nrsign(max(abs(q - r), TINY), q - r))
        ulim = bx + GLIMIT * (cx - bx)
        if (bx - u) * (u - cx) > 0.0
            fu = f(u)
            if fu < fc
                ax = bx
                bx = u
                fa = fb
                fb = fu
                return ax, bx, cx, fa, fb, fc
            elseif fu > fb
                cx = u
                fc = fu
                return ax, bx, cx, fa, fb, fc
            end
            u = cx + GOLD * (cx - bx)
            fu = f(u)
        elseif (cx - u) * (u - ulim) > 0.0
            fu = f(u)
            if fu < fc
                bx = cx
                cx = u
                u = cx + GOLD * (cx - bx)
                fb = fc
                fc = fu
                fu = f(u)
            end
        elseif (u - ulim) * (ulim - cx) >= 0.0
            u = ulim
            fu = f(u)
        else
            u = cx + GOLD * (cx - bx)
            fu = f(u)
        end
        ax = bx
        bx = cx
        cx = u
        fa = fb
        fb = fc
        fc = fu
    end
    return ax, bx, cx, fa, fb, fc
end

function brent_min(f, ax, bx, cx, tol)
    a = ax < cx ? ax : cx
    b = ax > cx ? ax : cx
    x = w = v = bx
    fw = fv = fx = f(x)
    e = 0.0
    d = 0.0
    for iter = 1:100
        xm = 0.5 * (a + b)
        tol1 = tol * abs(x) + ZEPS
        tol2 = 2.0 * tol1
        if abs(x - xm) <= (tol2 - 0.5 * (b - a))
            return fx, x
        end
        if abs(e) > tol1
            r = (x - w) * (fx - fv)
            q = (x - v) * (fx - fw)
            p = (x - v) * q - (x - w) * r
            q = 2.0 * (q - r)
            q > 0.0 && (p = -p)
            q = abs(q)
            etemp = e
            e = d
            if abs(p) >= abs(0.5 * q * etemp) || p <= q * (a - x) || p >= q * (b - x)
                e = x >= xm ? a - x : b - x
                d = CGOLD * e
            else
                d = p / q
                u = x + d
                (u - a < tol2 || b - u < tol2) && (d = nrsign(tol1, xm - x))
            end
        else
            e = x >= xm ? a - x : b - x
            d = CGOLD * e
        end
        u = abs(d) >= tol1 ? x + d : x + nrsign(tol1, d)
        fu = f(u)
        if fu <= fx
            u >= x ? (a = x) : (b = x)
            v = w
            w = x
            x = u
            fv = fw
            fw = fx
            fx = fu
        else
            u < x ? (a = u) : (b = u)
            if fu <= fw || w == x
                v = w
                w = u
                fv = fw
                fw = fu
            elseif fu <= fv || v == x || v == w
                v = u
                fv = fu
            end
        end
    end
    glb_warning("glb_brent_min")
    return fx, x
end

function line_minimization!(p, xi, func)
    n = length(p)
    pcom = copy(p)
    xicom = copy(xi)
    xt = similar(p)
    f1d(t) = (for j = 1:n
        xt[j] = pcom[j] + t * xicom[j]
    end; func(xt))
    ax, xx, bx, _, _, _ = bracket(f1d, 0.0, 1.0)
    fret, xmin = brent_min(f1d, ax, xx, bx, BRENT_TOL)
    for j = 1:n
        xi[j] *= xmin
        p[j] += xi[j]
    end
    return fret
end

"""
    glb_powell!(p, xi, ftol, func) -> (status, fret, iter)

Powell's direction-set minimization of `func` starting from `p` with initial
directions given by the columns of `xi`. `p` is overwritten by the minimum.
"""
function glb_powell!(p::Vector{Float64}, xi::Matrix{Float64}, ftol, func)
    n = length(p)
    fret = func(p)
    pt = copy(p)
    ptt = similar(p)
    xit = similar(p)
    iter = 1
    while true
        fp = fret
        ibig = 0
        del = 0.0
        for i = 1:n
            for j = 1:n
                xit[j] = xi[j, i]
            end
            fptt = fret
            fret = line_minimization!(p, xit, func)
            if abs(fptt - fret) > del
                del = abs(fptt - fret)
                ibig = i
            end
        end
        if 2.0 * abs(fp - fret) <= ftol * (abs(fp) + abs(fret))
            return 0, fret, iter
        end
        if iter == 5000
            glb_warning("Minimizer could not converge!")
            return -1, fret, iter
        end
        for j = 1:n
            ptt[j] = 2.0 * p[j] - pt[j]
            xit[j] = p[j] - pt[j]
            pt[j] = p[j]
        end
        fptt = func(ptt)
        if fptt < fp
            # Powell's criterion for replacing a direction (Numerical Recipes).
            # C GLoBES writes it with the macro SQR(a) (sqrarg=(a),sqrarg*sqrarg),
            # which is undefined behaviour when used twice in one expression; gcc
            # computes (2(fp-2fret+fptt)-del)(fp-fptt)² instead. We use the
            # intended formula, so fits can take a different path than in C.
            t = 2.0 * (fp - 2.0 * fret + fptt) * (fp - fret - del)^2 - del * (fp - fptt)^2
            if t < 0.0
                fret = line_minimization!(p, xit, func)
                for j = 1:n
                    xi[j, ibig] = xit[j]
                end
            end
        end
        iter += 1
    end
end

identity_matrix(n) = Matrix{Float64}([i == j ? 1.0 : 0.0 for i = 1:n, j = 1:n])

const MAX_POWELL_RESTARTS = 10

"""
    glb_powell_restarted!(p, ftol, func) -> (status, fret)

Powell's method, restarted from the result with a fresh set of directions
until a whole run improves `func` by no more than `ftol` (relative).

Powell's stopping rule only looks at the improvement during one sweep, so in
long curved valleys it can stop well above the minimum (C GLoBES does). A
restart costs one more sweep or so when the first run had converged.
"""
function glb_powell_restarted!(p::Vector{Float64}, ftol, func)
    n = length(p)
    status, fret, _ = glb_powell!(p, identity_matrix(n), ftol, func)
    for _ = 1:MAX_POWELL_RESTARTS
        status == 0 || break
        fprev = fret
        status, fret, _ = glb_powell!(p, identity_matrix(n), ftol, func)
        2.0 * abs(fprev - fret) <= ftol * (abs(fprev) + abs(fret)) && break
    end
    return status, fret
end

# ---------------------------------------------------------------------------
#  Hybrid minimizer (oscillation and systematics parameters together)
# ---------------------------------------------------------------------------

function hybrid_eval_1d(f, P, u, xi, new_rates)
    ux = P .+ u .* xi
    return f(ux, new_rates)
end

function hybrid_bracket(f, P, xi, new_rates)
    ax = 0.0
    bx = 1.0
    fa = hybrid_eval_1d(f, P, ax, xi, new_rates)
    fb = hybrid_eval_1d(f, P, bx, xi, new_rates)
    if fb > fa
        ax, bx = bx, ax
        fa, fb = fb, fa
    end
    cx = bx + GOLD * (bx - ax)
    fc = hybrid_eval_1d(f, P, cx, xi, new_rates)
    while fb > fc
        r = (bx - ax) * (fb - fc)
        q = (bx - cx) * (fb - fa)
        u = bx - ((bx - cx) * q - (bx - ax) * r) / (2.0 * nrsign(max(abs(q - r), TINY), q - r))
        ulim = bx + GLIMIT * (cx - bx)
        if (bx - u) * (u - cx) > 0.0
            fu = hybrid_eval_1d(f, P, u, xi, new_rates)
            if fu < fc
                ax = bx
                fa = fb
                bx = u
                fb = fu
                break
            elseif fu > fb
                cx = u
                fc = fu
                break
            else
                u = cx + GOLD * (cx - bx)
                fu = hybrid_eval_1d(f, P, u, xi, new_rates)
            end
        elseif (cx - u) * (u - ulim) > 0.0
            fu = hybrid_eval_1d(f, P, u, xi, new_rates)
            if fu < fc
                bx = cx
                cx = u
                u = cx + GOLD * (cx - bx)
                fb = fc
                fc = fu
                fu = hybrid_eval_1d(f, P, u, xi, new_rates)
            end
        elseif (u - ulim) * (ulim - cx) >= 0.0
            u = ulim
            fu = hybrid_eval_1d(f, P, u, xi, new_rates)
        else
            u = cx + GOLD * (cx - bx)
            fu = hybrid_eval_1d(f, P, u, xi, new_rates)
        end
        ax = bx
        bx = cx
        cx = u
        fa = fb
        fb = fc
        fc = fu
    end
    a = ax < cx ? ax : cx
    c = ax > cx ? ax : cx
    return a, bx, c
end

function hybrid_linmin!(P, xi, ftol, f, new_rates)
    a, x, b = hybrid_bracket(f, P, xi, new_rates)
    g(t) = hybrid_eval_1d(f, P, t, xi, new_rates)
    e = 0.0
    d = 0.0
    w = v = x
    fw = fv = fx = g(x)
    for iter = 1:100
        xm = 0.5 * (a + b)
        tol1 = ftol * abs(x) + ZEPS
        tol2 = 2.0 * tol1
        if abs(x - xm) <= tol2 - 0.5 * (b - a)
            P .+= x .* xi
            return fx, 0
        end
        if abs(e) > tol1
            r = (x - w) * (fx - fv)
            q = (x - v) * (fx - fw)
            p = (x - v) * q - (x - w) * r
            q = 2.0 * (q - r)
            q > 0.0 && (p = -p)
            q = abs(q)
            etemp = e
            e = d
            if abs(p) >= abs(0.5 * q * etemp) || p <= q * (a - x) || p >= q * (b - x)
                e = x >= xm ? a - x : b - x
                d = CGOLD * e
            else
                d = p / q
                u = x + d
                (u - a < tol2 || b - u < tol2) && (d = nrsign(tol1, xm - x))
            end
        else
            e = x >= xm ? a - x : b - x
            d = CGOLD * e
        end
        u = abs(d) >= tol1 ? x + d : x + nrsign(tol1, d)
        fu = g(u)
        if fu <= fx
            u >= x ? (a = x) : (b = x)
            v = w
            w = x
            x = u
            fv = fw
            fw = fx
            fx = fu
        else
            u < x ? (a = u) : (b = u)
            if fu <= fw || w == x
                v = w
                fv = fw
                w = u
                fw = fu
            elseif fu <= fv || v == x || v == w
                v = u
                fv = fu
            end
        end
    end
    glb_warning("glb_hybrid_linmin: No convergence in 1-dim Brent minimizer")
    P .+= x .* xi
    return fx, -1
end

function iterate_hybrid!(P, xi, n_min, n_max, ftol, fret, f, new_rates)
    n = length(P)
    fp = fret
    pt = copy(P)
    ibig = 0
    del = 0.0
    for i = n_min:n_max
        xit = xi[:, i]
        fptt = fret
        fret, _ = hybrid_linmin!(P, xit, BRENT_TOL, f, new_rates)
        if fptt - fret > del
            del = fptt - fret
            ibig = i
        end
    end
    if 2.0 * (fp - fret) <= ftol * (abs(fp) + abs(fret))
        return 0, fret
    end
    ptt = 2.0 .* P .- pt
    xit = P .- pt
    fptt = f(ptt, new_rates)
    if fptt < fp
        t = 2.0 * (fp - 2.0 * fret + fptt) * (fp - fret - del)^2 - del * (fp - fptt)^2
        if t < 0.0
            fret, _ = hybrid_linmin!(P, xit, BRENT_TOL, f, new_rates)
            xi[:, ibig] .= xi[:, n_max]
            xi[:, n_max] .= xit
        end
    end
    return 1, fret
end

"""
    glb_hybrid_minimizer!(P, n_osc, ftol, f) -> (status, fret, iter)

Powell minimizer that alternates between the first `n_osc` (oscillation)
directions, for which `f(x, true)` recomputes event rates, and the remaining
(systematics) directions, evaluated with `f(x, false)`.
"""
function glb_hybrid_minimizer!(P::Vector{Float64}, n_osc, ftol, f)
    n = length(P)
    xi = zeros(n, n)
    for i = 1:n
        xi[i, i] = abs(P[i]) > ftol ? 0.1 * abs(P[i]) : 1.0
    end
    fret = f(P, true)
    iter = 1
    while iter < n * 100
        converged = true
        s, fret = iterate_hybrid!(P, xi, 1, n_osc, ftol, fret, f, true)
        s > 0 && (converged = false)
        fret = f(P, true)
        s, fret = iterate_hybrid!(P, xi, n_osc + 1, n, ftol, fret, f, false)
        s > 0 && (converged = false)
        converged && return 0, fret, iter
        iter += 1
    end
    glb_warning("glb_hybrid_minimizer: No convergence in n-dim Powell minimizer")
    return -1, fret, iter
end

# ---------------------------------------------------------------------------
#  Projections
# ---------------------------------------------------------------------------

n_osc() = state().oscp
n_exps() = length(state().experiments)

function select_projection!(vec::AbstractVector{Int}, check::Bool)
    st = state()
    N = n_osc() + n_exps()
    resize!(st.para_tab, N)
    bad = false
    for i = 1:N
        if check && vec[i] != GLB_FREE && vec[i] != GLB_FIXED
            bad = true
            st.para_tab[i] = GLB_FREE
        else
            st.para_tab[i] = vec[i]
        end
    end
    bad && glb_error("SelectProjection: Projection partly undefined. Using default GLB_FREE")
    free = [i for i = 1:N if st.para_tab[i] == GLB_FREE]
    fixed = [i for i = 1:N if st.para_tab[i] == GLB_FIXED]
    st.index_tab = vcat(free, fixed)
    st.n_free = length(free)
    st.n_fix = length(fixed)
end

"""
    glbSetProjection(proj)

Select which parameters are minimized over (`GLB_FREE`) in `glbChiNP`.
"""
function glbSetProjection(proj::GlbProjection)
    select_projection!(vcat(proj.osc, proj.density), true)
    return 0
end

"""
    glbGetProjection(proj)

Copy the current projection into `proj`.
"""
function glbGetProjection(proj::GlbProjection)
    st = state()
    length(st.para_tab) < n_osc() + n_exps() && reset_projection!()
    proj.osc = st.para_tab[1:n_osc()]
    proj.density = st.para_tab[(n_osc() + 1):(n_osc() + n_exps())]
    return 0
end

function reset_projection!()
    st = state()
    st.para_tab = fill(GLB_UNDEFINED, n_osc() + n_exps())
    st.index_tab = Int[]
    st.n_free = 0
    st.n_fix = 0
end

"""
    glbSelectMinimizer(id)

`GLB_MIN_NESTED_POWELL` (default) or `GLB_MIN_POWELL` (hybrid).
"""
function glbSelectMinimizer(id)
    (id == GLB_MIN_NESTED_POWELL || id == GLB_MIN_POWELL) ||
        (glb_error("glbSelectMinimizer: Invalid minimization algorithm requested"); return -1)
    state().minimizer = id
    return 0
end

# ---------------------------------------------------------------------------
#  Nested minimization: systematics inside, oscillation parameters outside
# ---------------------------------------------------------------------------

function chi_rule(exp, rule)
    e = state().experiments[exp]
    r = e.rules[rule]
    sys = active_sys(r)
    dim = sys.dim
    errors = something(active_errors(r), Float64[])
    sp = dim > 0 ? copy(active_startvals(r)) : Float64[]
    f = x -> sys.func(exp, rule, dim, x, errors, sys.user_data)
    status, res, _ = glb_powell!(sp, identity_matrix(dim), TOLSYS, f)
    if status != 0
        glb_warning("Systematics minimization failed.")
        return -res
    end
    return res
end

chi_exp(exp) =
    sum((chi_rule(exp, j) for j in eachindex(state().experiments[exp].rules)); init = 0.0)

function set_osc_from_vector!(x)
    p = glbAllocParams()
    p.osc .= x[1:n_osc()]
    state().set_func(p, state().prob_user_data)
end

# chi^2 over all experiments for oscillation params x[1:oscp] and densities x[oscp+1:end]
function chi_all_exps(x)
    set_osc_from_vector!(x)
    for i = 1:n_exps()
        e = state().experiments[i]
        set_profile_scaling!(e, x[n_osc() + i])
        set_new_rates!(e)
    end
    return sum((chi_exp(i) for i = 1:n_exps()); init = 0.0)
end

function chi_single_exp(x, exp, rule = GLB_ALL)
    set_osc_from_vector!(x)
    e = state().experiments[exp]
    set_profile_scaling!(e, x[n_osc() + 1])
    set_new_rates!(e)
    rule == GLB_ALL && return chi_exp(exp)
    return chi_rule(exp, rule)
end

function density_or_default(in::GlbParams, i, msg)
    v = i <= length(in.density) ? in.density[i] : NaN
    if isnan(v)
        glb_warning("$msg: Density parameters incompletely defined. Using default 1.0")
        v = 1.0
    end
    return v
end

"""
    glbChiSys(params, exp, rule) -> Float64

χ² for fixed oscillation parameters, minimized over systematics only.
`exp` and `rule` may be `GLB_ALL`.
"""
function glbChiSys(in::GlbParams, exp::Integer = GLB_ALL, rule::Integer = GLB_ALL)
    st = state()
    x = zeros(n_osc() + max(n_exps(), 1))
    for i = 1:n_osc()
        x[i] = in.osc[i]
        isnan(x[i]) &&
            (glb_error("glbChiSys: Oscillation parameters incompletely defined"); return NaN)
    end
    if st.minimizer == GLB_MIN_POWELL
        for i = 1:n_exps()
            x[n_osc() + i] =
                (exp == GLB_ALL || i == exp) ? density_or_default(in, i, "glbChiSys") :
                (i <= length(in.density) ? in.density[i] : 1.0)
        end
        old = copy(st.para_tab)
        select_projection!(fill(GLB_FIXED, n_osc() + n_exps()), false)
        res = invoke_hybrid_minimizer!(exp, rule, x)
        select_projection!(old, false)
        return res
    end
    if exp == GLB_ALL
        if rule == GLB_ALL
            for i = 1:n_exps()
                x[n_osc() + i] = density_or_default(in, i, "glbChiSys")
            end
            return chi_all_exps(x)
        else
            res = 0.0
            for i = 1:n_exps()
                rule <= length(st.experiments[i].rules) ||
                    (glb_error("Invalid rule number"); return -1.0)
                x[n_osc() + 1] = density_or_default(in, i, "glbChiSys")
                res += chi_single_exp(x, i, rule)
            end
            return res
        end
    end
    (1 <= exp <= n_exps()) ||
        (glb_error("Failure in glbChiSys: invalid experiment number"); return -1.0)
    x[n_osc() + 1] = density_or_default(in, exp, "glbChiSys")
    if rule != GLB_ALL
        rule <= length(st.experiments[exp].rules) || (
            glb_error("Failure in glbChiSys: 3rd argument must not exceed numofrules");
            return -1.0
        )
    end
    return chi_single_exp(x, exp, rule)
end

function fill_prior_defaults!(proj::GlbProjection, msg)
    st = state()
    op = glbAllocParams()
    glbGetOscillationParameters(op)
    for i = 1:n_osc()
        if proj.osc[i] == GLB_FREE
            if isnan(st.input_errors.osc[i])
                glb_error("$msg: Input errors incompletely defined. Using default 0.0")
                st.input_errors.osc[i] = 0.0
            end
            if isnan(st.central_values.osc[i])
                glb_error(
                    "$msg: Central values incompletely defined. Using default (curr. osc. params)",
                )
                st.central_values.osc[i] = op.osc[i]
            end
        end
    end
    for i = 1:n_exps()
        if proj.density[i] == GLB_FREE
            if isnan(st.input_errors.density[i])
                glb_error("$msg: Density input errors incompletely defined. Using default 0.05")
                st.input_errors.density[i] = 0.05
            end
            if isnan(st.central_values.density[i])
                glb_warning("$msg: Density central values incompletely defined. Using default 1.0")
                st.central_values.density[i] = 1.0
            end
        end
    end
    st.input_errors_hook(st.input_errors, st.prior_user_data) != 0 && return false
    st.starting_values_hook(st.central_values, st.prior_user_data) != 0 && return false
    return true
end

function internal_chi_np(in::GlbParams, out)
    st = state()
    st.single_experiment = GLB_ALL
    select_projection!(st.para_tab, true)
    out !== nothing && glbCopyParams(in, out)
    N = n_osc() + n_exps()
    x = zeros(N)
    for i = 1:n_osc()
        x[i] = in.osc[i]
        isnan(x[i]) &&
            (glb_error("glbChiNP: Oscillation parameters incompletely defined"); return NaN)
    end
    for i = 1:n_exps()
        x[n_osc() + i] = density_or_default(in, i, "glbChiNP")
    end
    proj = glbAllocProjection()
    glbGetProjection(proj)
    fill_prior_defaults!(proj, "glbChiNP") || return NaN

    if st.minimizer == GLB_MIN_POWELL
        st.count = 0
        res = invoke_hybrid_minimizer!(GLB_ALL, GLB_ALL, x)
        if out !== nothing
            out.osc .= x[1:n_osc()]
            out.density .= x[(n_osc() + 1):end]
            out.iterations = st.count
        end
        return res
    end

    st.count = 0
    fix = copy(x)
    nf = st.n_free
    idx = copy(st.index_tab)
    sp = [x[idx[i]] for i = 1:nf]
    y = similar(x)
    prior_input = glbAllocParams()
    function md_chi_np(v)
        st.count += 1
        for i = 1:nf
            y[idx[i]] = v[i]
        end
        for i = (nf + 1):length(idx)
            y[idx[i]] = fix[idx[i]]
        end
        erg = chi_all_exps(y)
        prior_input.osc .= y[1:n_osc()]
        prior_input.density .= y[(n_osc() + 1):end]
        prior_input.iterations = st.count
        return erg + st.prior_func(prior_input, st.prior_user_data)
    end
    status, er1 = glb_powell_restarted!(sp, TOLOSC, md_chi_np)
    status != 0 && (st.count = -st.count)
    if out !== nothing
        for i = 1:nf
            if idx[i] <= n_osc()
                out.osc[idx[i]] = sp[i]
            else
                out.density[idx[i] - n_osc()] = sp[i]
            end
        end
        out.iterations = st.count
    end
    return er1
end

function internal_single_chi_np(in::GlbParams, out, exp)
    st = state()
    (1 <= exp <= n_exps()) || (glb_error("glbChiNP: invalid experiment number"); return -1.0)
    st.single_experiment = exp
    # Projection restricted to this experiment's density
    length(st.para_tab) < n_osc() + n_exps() && reset_projection!()
    s_para = vcat(st.para_tab[1:n_osc()], st.para_tab[n_osc() + exp])
    bad = false
    for i in eachindex(s_para)
        if s_para[i] != GLB_FREE && s_para[i] != GLB_FIXED
            s_para[i] = GLB_FREE
            bad = true
        end
    end
    bad && glb_error("single_SelectProjection: Projection partly undefined. Using default GLB_FREE")
    st.para_tab[1:n_osc()] .= s_para[1:n_osc()]
    st.para_tab[n_osc() + exp] = s_para[end]
    s_idx = vcat(
        [i for i in eachindex(s_para) if s_para[i] == GLB_FREE],
        [i for i in eachindex(s_para) if s_para[i] == GLB_FIXED],
    )
    s_nfree = count(==(GLB_FREE), s_para)

    fbuf = glbAllocProjection()
    glbGetProjection(fbuf)
    fnew = glbCopyProjection(fbuf, glbAllocProjection())
    glbSetDensityProjectionFlag(fnew, GLB_FIXED, GLB_ALL)
    glbSetDensityProjectionFlag(fnew, glbGetDensityProjectionFlag(fbuf, exp), exp)
    glbSetProjection(fnew)

    out !== nothing && glbCopyParams(in, out)
    x = zeros(n_osc() + max(n_exps(), 1))
    for i = 1:n_osc()
        x[i] = in.osc[i]
        isnan(x[i]) &&
            (glb_error("glbChiNP: Oscillation parameters incompletely defined"); return NaN)
    end
    if st.minimizer == GLB_MIN_POWELL
        for i = 1:n_exps()
            x[n_osc() + i] =
                i == exp ? density_or_default(in, i, "glbChiNP") :
                (i <= length(in.density) ? in.density[i] : NaN)
        end
    else
        x[n_osc() + 1] = density_or_default(in, exp, "glbChiNP")
    end
    fill_prior_defaults!(fnew, "glbChiNP") || return NaN

    local er1
    if st.minimizer == GLB_MIN_POWELL
        st.count = 0
        er1 = invoke_hybrid_minimizer!(exp, GLB_ALL, x)
        if out !== nothing
            out.osc .= x[1:n_osc()]
            out.density .= x[(n_osc() + 1):(n_osc() + n_exps())]
            out.iterations = st.count
        end
    else
        st.count = 0
        fix = x[1:(n_osc() + 1)]
        sp = [fix[s_idx[i]] for i = 1:s_nfree]
        y = similar(fix)
        prior_input = glbAllocParams()
        function chi_np(v)
            st.count += 1
            for i = 1:s_nfree
                y[s_idx[i]] = v[i]
            end
            for i = (s_nfree + 1):length(s_idx)
                y[s_idx[i]] = fix[s_idx[i]]
            end
            erg = chi_single_exp(y, exp)
            prior_input.osc .= y[1:n_osc()]
            fill!(prior_input.density, NaN)
            prior_input.density[exp] = y[n_osc() + 1]
            prior_input.iterations = st.count
            return erg + st.prior_func(prior_input, st.prior_user_data)
        end
        status, er1 = glb_powell_restarted!(sp, TOLOSC, chi_np)
        status != 0 && (st.count = -st.count)
        if out !== nothing
            for i = 1:s_nfree
                if s_idx[i] <= n_osc()
                    out.osc[s_idx[i]] = sp[i]
                else
                    out.density[exp] = sp[i]
                end
            end
            out.iterations = st.count
        end
    end
    st.single_experiment = GLB_ALL
    glbSetProjection(fbuf)
    return er1
end

"""
    glbChiNP(params, out, exp) -> Float64

χ² minimized over all `GLB_FREE` parameters of the current projection
(see `glbSetProjection`) and over systematics. The position of the minimum
is written to `out` (pass `nothing` to ignore it). `exp` may be `GLB_ALL`.
"""
function glbChiNP(in::GlbParams, out::Union{Nothing, GlbParams} = nothing, exp::Integer = GLB_ALL)
    exp == GLB_ALL ? internal_chi_np(in, out) : internal_single_chi_np(in, out, exp)
end

function chi_with_projection(in, out, exp, swit)
    st = state()
    length(st.para_tab) < n_osc() + n_exps() && reset_projection!()
    buff = copy(st.para_tab)
    select_projection!(swit, true)
    res = exp == GLB_ALL ? internal_chi_np(in, out) : internal_single_chi_np(in, out, exp)
    select_projection!(buff, false)
    return res
end

function chi_1p(in, out, exp, n)
    swit = fill(GLB_FREE, n_osc() + n_exps())
    swit[n] = GLB_FIXED
    return chi_with_projection(in, out, exp, swit)
end

"""
    glbChiTheta12(params, out, exp)

χ² with θ12 fixed and all other parameters (including matter density) free.
"""
glbChiTheta12(in, out = nothing, exp = GLB_ALL) = chi_1p(in, out, exp, GLB_THETA_12)
"""
    glbChiTheta13(params, out, exp)

χ² with θ13 fixed and all other parameters free.
"""
glbChiTheta13(in, out = nothing, exp = GLB_ALL) = chi_1p(in, out, exp, GLB_THETA_13)
"""
    glbChiTheta23(params, out, exp)

χ² with θ23 fixed and all other parameters free.
"""
glbChiTheta23(in, out = nothing, exp = GLB_ALL) = chi_1p(in, out, exp, GLB_THETA_23)
"""
    glbChiDelta(params, out, exp)

χ² with δCP fixed and all other parameters free.
"""
glbChiDelta(in, out = nothing, exp = GLB_ALL) = chi_1p(in, out, exp, GLB_DELTA_CP)
"""
    glbChiDm21(params, out, exp)

χ² with Δm²21 fixed and all other parameters free.
"""
glbChiDm21(in, out = nothing, exp = GLB_ALL) = chi_1p(in, out, exp, GLB_DM_21)
"""
    glbChiDm31(params, out, exp)

χ² with Δm²31 fixed and all other parameters free.
"""
glbChiDm31(in, out = nothing, exp = GLB_ALL) = chi_1p(in, out, exp, GLB_DM_31)

"""
    glbChiTheta13Delta(params, out, exp)

χ² with θ13 and δCP fixed and all other parameters free.
"""
function glbChiTheta13Delta(in, out = nothing, exp = GLB_ALL)
    swit = fill(GLB_FREE, n_osc() + n_exps())
    swit[GLB_THETA_13] = GLB_FIXED
    swit[GLB_DELTA_CP] = GLB_FIXED
    return chi_with_projection(in, out, exp, swit)
end

"""
    glbChiAll(params, out, exp)

χ² minimized over all oscillation parameters and densities.
"""
glbChiAll(in, out = nothing, exp = GLB_ALL) =
    chi_with_projection(in, out, exp, fill(GLB_FREE, n_osc() + n_exps()))

# ---------------------------------------------------------------------------
#  Hybrid minimizer driver
# ---------------------------------------------------------------------------

function hybrid_chi_callback(x, new_rates)
    st = state()
    if new_rates
        p = glbAllocParams()
        j = 1
        for i = 1:n_osc()
            p.osc[i] = st.para_tab[i] == GLB_FREE ? x[j] : st.fix_params[i]
            st.para_tab[i] == GLB_FREE && (j += 1)
        end
        for k = 1:n_exps()
            i = n_osc() + k
            p.density[k] = st.para_tab[i] == GLB_FREE ? x[j] : st.fix_params[i]
            st.para_tab[i] == GLB_FREE && (j += 1)
        end
        st.count += 1
        p.iterations = st.count
        st.set_func(p, st.prob_user_data)
        for i = 1:n_exps()
            (st.single_experiment != GLB_ALL && st.single_experiment != i) && continue
            e = st.experiments[i]
            set_profile_scaling!(e, p.density[i])
            set_new_rates!(e)
        end
        st.last_priors = st.prior_func(p, st.prior_user_data)
    end
    m = st.n_free + 1
    chi2 = 0.0
    for i = 1:n_exps()
        (st.single_experiment != GLB_ALL && st.single_experiment != i) && continue
        e = st.experiments[i]
        for j in eachindex(e.rules)
            (st.single_rule != GLB_ALL && st.single_rule != j) && continue
            r = e.rules[j]
            sys = active_sys(r)
            y = x[m:(m + sys.dim - 1)]
            m += sys.dim
            chi2 +=
                sys.func(i, j, sys.dim, y, something(active_errors(r), Float64[]), sys.user_data)
        end
    end
    return chi2 + st.last_priors
end

function invoke_hybrid_minimizer!(exp, rule, x)
    st = state()
    N = n_osc() + n_exps()
    st.fix_params = copy(x[1:N])
    P = Float64[]
    for i = 1:N
        st.para_tab[i] == GLB_FREE && push!(P, x[i])
    end
    for i = 1:n_exps()
        (exp != GLB_ALL && exp != i) && continue
        e = st.experiments[i]
        for j in eachindex(e.rules)
            (rule != GLB_ALL && rule != j) && continue
            sv = active_startvals(e.rules[j])
            sv !== nothing && append!(P, sv)
        end
    end
    st.single_experiment = exp
    st.single_rule = rule
    _, chi2, _ = glb_hybrid_minimizer!(P, st.n_free, TOLOSC, hybrid_chi_callback)
    k = 1
    for i = 1:N
        if st.para_tab[i] == GLB_FREE
            x[i] = P[k]
            k += 1
        end
    end
    st.single_experiment = GLB_ALL
    st.single_rule = GLB_ALL
    return chi2
end
