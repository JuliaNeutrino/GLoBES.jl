# Three-flavour oscillation probabilities in matter of piecewise constant
# density. Port of glb_probability.c, including the fast 3x3 Hermitian
# eigensolvers of J. Kopp, physics/0610206.

# Matter potential V = sqrt(2) G_F N_e: GLB_V_FACTOR = sqrt(2) G_F N_A (ħc)^3 in
# eV per g/cm^3 (CODATA 2018; C GLoBES uses 7.5e-14), times the electron
# fraction Y_e = 0.5
const GLB_V_FACTOR = sqrt(2) * 1.1663787e-23 * 6.02214076e23 * (1.973269804e-5)^3
const GLB_Ne_MANTLE = 0.5
const V_THRESHOLD = 0.001 * GLB_V_FACTOR * GLB_Ne_MANTLE
const M_SQRT3 = 1.73205080756887729352744634151

@inline sqr(x) = x * x
@inline sqrabs(z::Complex) = sqr(real(z)) + sqr(imag(z))

# ---------------------------------------------------------------------------
#  3 x 3   E I G E N S Y S T E M   F U N C T I O N S
# ---------------------------------------------------------------------------
# All matrices are 3x3 Matrix{ComplexF64}; only the diagonal and the upper
# triangle of A are accessed.

function zhetrd3!(A, Q, d, e)
    n = 3
    u1 = zero(ComplexF64)
    u2 = zero(ComplexF64)
    q1 = zero(ComplexF64)
    q2 = zero(ComplexF64)

    for i = 1:n, j = 1:n
        Q[i, j] = i == j ? 1.0 : 0.0
    end

    h = sqrabs(A[1, 2]) + sqrabs(A[1, 3])
    g = real(A[1, 2]) > 0 ? -sqrt(h) : sqrt(h)
    e[1] = g
    f = g * A[1, 2]
    u1 = conj(A[1, 2]) - g
    u2 = conj(A[1, 3])

    omega = h - f
    if real(omega) > 0.0
        omega = 0.5 * (1.0 + conj(omega) / omega) / real(omega)
        K = 0.0
        # i = 2
        f = conj(A[2, 2]) * u1 + A[2, 3] * u2
        q1 = omega * f
        K += real(conj(u1) * f)
        # i = 3
        f = conj(A[2, 3]) * u1 + A[3, 3] * u2
        q2 = omega * f
        K += real(conj(u2) * f)
        K *= 0.5 * sqrabs(omega)

        q1 = q1 - K * u1
        q2 = q2 - K * u2

        d[1] = real(A[1, 1])
        d[2] = real(A[2, 2]) - 2.0 * real(q1 * conj(u1))
        d[3] = real(A[3, 3]) - 2.0 * real(q2 * conj(u2))

        # Store inverse Householder transformation in Q
        f = omega * conj(u1)
        Q[2, 2] = Q[2, 2] - f * u1
        Q[3, 2] = Q[3, 2] - f * u2
        f = omega * conj(u2)
        Q[2, 3] = Q[2, 3] - f * u1
        Q[3, 3] = Q[3, 3] - f * u2

        f = A[2, 3] - q1 * conj(u2) - u1 * conj(q2)
    else
        for i = 1:n
            d[i] = real(A[i, i])
        end
        f = A[2, 3]
    end

    # Make (23) element real
    e[2] = abs(f)
    if e[2] != 0.0
        f = conj(f) / e[2]
        Q[2, 3] = Q[2, 3] * f
        Q[3, 3] = Q[3, 3] * f
    end
    return nothing
end

function zheevc3!(A, w)
    de = A[1, 2] * A[2, 3]
    dd = sqrabs(A[1, 2])
    ee = sqrabs(A[2, 3])
    ff = sqrabs(A[1, 3])
    a = real(A[1, 1])
    b = real(A[2, 2])
    c = real(A[3, 3])
    m = a + b + c
    c1 = (a * b + a * c + b * c) - (dd + ee + ff)
    c0 =
        c * dd + a * ee + b * ff - a * b * c -
        2.0 * (real(A[1, 3]) * real(de) + imag(A[1, 3]) * imag(de))

    p = sqr(m) - 3.0 * c1
    q = m * (p - (3.0 / 2.0) * c1) - (27.0 / 2.0) * c0
    sqrt_p = sqrt(abs(p))

    phi = 27.0 * (0.25 * sqr(c1) * (p - c1) + c0 * (q + 27.0 / 4.0 * c0))
    phi = (1.0 / 3.0) * atan(sqrt(abs(phi)), q)

    cc = sqrt_p * cos(phi)
    s = (1.0 / M_SQRT3) * sqrt_p * sin(phi)

    w[2] = (1.0 / 3.0) * (m - cc)
    w[3] = w[2] + s
    w[1] = w[2] + cc
    w[2] -= s
    return 0
end

function zheevq3!(A, Q, w)
    n = 3
    e = MVector3()
    zhetrd3!(A, Q, w, e)

    for l = 1:(n - 1)
        nIter = 0
        while true
            m = l
            while m <= n - 1
                g = abs(w[m]) + abs(w[m + 1])
                if abs(e[m]) + g == g
                    break
                end
                m += 1
            end
            if m == l
                break
            end

            nIter += 1
            nIter > 30 && return -1

            g = (w[l + 1] - w[l]) / (e[l] + e[l])
            r = sqrt(sqr(g) + 1.0)
            if g > 0
                g = w[m] - w[l] + e[l] / (g + r)
            else
                g = w[m] - w[l] + e[l] / (g - r)
            end

            s = 1.0
            c = 1.0
            p = 0.0
            for i = (m - 1):-1:l
                f = s * e[i]
                b = c * e[i]
                if abs(f) > abs(g)
                    c = g / f
                    r = sqrt(sqr(c) + 1.0)
                    e[i + 1] = f * r
                    s = 1.0 / r
                    c *= s
                else
                    s = f / g
                    r = sqrt(sqr(s) + 1.0)
                    e[i + 1] = g * r
                    c = 1.0 / r
                    s *= c
                end

                g = w[i + 1] - p
                r = (w[i] - g) * s + 2.0 * c * b
                p = s * r
                w[i + 1] = g + p
                g = c * r - b

                for k = 1:n
                    t = Q[k, i + 1]
                    Q[k, i + 1] = s * Q[k, i] + c * t
                    Q[k, i] = c * Q[k, i] - s * t
                end
            end
            w[l] -= p
            e[l] = g
            e[m] = 0.0
        end
    end
    return 0
end

# Small fixed-size scratch vector (avoids depending on StaticArrays)
MVector3() = zeros(Float64, 3)

function zheevh3!(A, Q, w)
    zheevc3!(A, w)

    t = abs(w[1])
    (u = abs(w[2])) > t && (t = u)
    (u = abs(w[3])) > t && (t = u)
    u = t < 1.0 ? t : sqr(t)
    err = 256.0 * eps(Float64) * sqr(u)

    Q[1, 2] = A[1, 2] * A[2, 3] - A[1, 3] * real(A[2, 2])
    Q[2, 2] = A[1, 3] * conj(A[1, 2]) - A[2, 3] * real(A[1, 1])
    Q[3, 2] = sqrabs(A[1, 2])

    # First eigenvector: v[0] = conj( (A - w[0]).e1 x (A - w[0]).e2 )
    Q[1, 1] = Q[1, 2] + A[1, 3] * w[1]
    Q[2, 1] = Q[2, 2] + A[2, 3] * w[1]
    Q[3, 1] = (real(A[1, 1]) - w[1]) * (real(A[2, 2]) - w[1]) - Q[3, 2]
    nrm = sqrabs(Q[1, 1]) + sqrabs(Q[2, 1]) + sqr(real(Q[3, 1]))

    if nrm <= err
        return zheevq3!(A, Q, w)
    else
        nrm = sqrt(1.0 / nrm)
        for j = 1:3
            Q[j, 1] = Q[j, 1] * nrm
        end
    end

    # Second eigenvector
    Q[1, 2] = Q[1, 2] + A[1, 3] * w[2]
    Q[2, 2] = Q[2, 2] + A[2, 3] * w[2]
    Q[3, 2] = (real(A[1, 1]) - w[2]) * (real(A[2, 2]) - w[2]) - real(Q[3, 2])
    nrm = sqrabs(Q[1, 2]) + sqrabs(Q[2, 2]) + sqr(real(Q[3, 2]))
    if nrm <= err
        return zheevq3!(A, Q, w)
    else
        nrm = sqrt(1.0 / nrm)
        for j = 1:3
            Q[j, 2] = Q[j, 2] * nrm
        end
    end

    # Third eigenvector: v[2] = conj(v[0] x v[1])
    Q[1, 3] = conj(Q[2, 1] * Q[3, 2] - Q[3, 1] * Q[2, 2])
    Q[2, 3] = conj(Q[3, 1] * Q[1, 2] - Q[1, 1] * Q[3, 2])
    Q[3, 3] = conj(Q[1, 1] * Q[2, 2] - Q[2, 1] * Q[1, 2])
    return 0
end

# ---------------------------------------------------------------------------
#  Standard three-flavour engine
# ---------------------------------------------------------------------------

mutable struct StandardEngine
    th12::Float64
    th13::Float64
    th23::Float64
    delta::Float64
    mq::Vector{Float64}
    U::Matrix{ComplexF64}           # vacuum mixing matrix
    H::Matrix{ComplexF64}           # Hamiltonian
    Q::Matrix{ComplexF64}           # eigenvectors of H
    lambda::Vector{Float64}         # eigenvalues of H
    S::Matrix{ComplexF64}           # S-matrix
    H0_template::Matrix{ComplexF64} # U . diag(m^2/2) . U^dagger
    S1::Matrix{ComplexF64}
    T0::Matrix{ComplexF64}
end

function StandardEngine()
    z() = zeros(ComplexF64, 3, 3)
    StandardEngine(0.0, 0.0, 0.0, 0.0, zeros(3), z(), z(), z(), zeros(3), z(), z(), z(), z())
end

function set_oscillation_parameters!(eng::StandardEngine, osc::AbstractVector{<:Real})
    th12, th13, th23, delta = osc[1], osc[2], osc[3], osc[4]
    eng.th12 = th12
    eng.th13 = th13
    eng.th23 = th23
    eng.delta = delta
    eng.mq[1] = 0.0
    eng.mq[2] = osc[5]
    eng.mq[3] = osc[6]

    U = eng.U
    s12, c12 = sincos(th12)
    s13, c13 = sincos(th13)
    s23, c23 = sincos(th23)
    ep = cis(delta)
    em = cis(-delta)
    U[1, 1] = c12 * c13
    U[1, 2] = s12 * c13
    U[1, 3] = s13 * em
    U[2, 1] = -s12 * c23 - c12 * s23 * s13 * ep
    U[2, 2] = c12 * c23 - s12 * s23 * s13 * ep
    U[2, 3] = s23 * c13
    U[3, 1] = s12 * s23 - c12 * c23 * s13 * ep
    U[3, 2] = -c12 * s23 - s12 * c23 * s13 * ep
    U[3, 3] = c23 * c13

    D = zeros(ComplexF64, 3, 3)
    for i = 1:3
        D[i, i] = 0.5 * eng.mq[i]
    end
    eng.H0_template .= U * D * U'
    fill!(eng.H, 0)
    return 0
end

function get_oscillation_parameters(eng::StandardEngine)
    return [eng.th12, eng.th13, eng.th23, eng.delta, eng.mq[2] - eng.mq[1], eng.mq[3] - eng.mq[1]]
end

function hamiltonian_cd!(eng::StandardEngine, E, V, cp_sign)
    inv_E = 1.0 / E
    H = eng.H
    H0 = eng.H0_template
    if cp_sign > 0
        for i = 1:3, j = 1:3
            H[i, j] = H0[i, j] * inv_E
        end
    else
        for i = 1:3, j = 1:3
            H[i, j] = conj(H0[i, j] * inv_E)   # delta_CP -> -delta_CP
        end
    end
    H[1, 1] = H[1, 1] + cp_sign * V
    return 0
end

# Puts eigenvalues into eng.lambda and eigenvectors into eng.Q
function diagonalize!(eng::StandardEngine, E, V, cp_sign)
    if V < V_THRESHOLD
        inv_E = 0.5 / E
        for i = 1:3
            eng.lambda[i] = eng.mq[i] * inv_E
        end
        if cp_sign > 0
            eng.Q .= eng.U
        else
            for i = 1:3, j = 1:3
                eng.Q[i, j] = conj(eng.U[i, j])
            end
        end
        return 0
    else
        hamiltonian_cd!(eng, E, V, cp_sign)
        return zheevh3!(eng.H, eng.Q, eng.lambda)
    end
end

function S_matrix_cd!(eng::StandardEngine, E, L, V, cp_sign)
    status = diagonalize!(eng, E, V, cp_sign)
    status != 0 && return status
    Q = eng.Q
    T0 = eng.T0
    S = eng.S
    # S = Q . diag(exp(-i L lambda)) . Q^dagger
    ph1 = cis(-L * eng.lambda[1])
    ph2 = cis(-L * eng.lambda[2])
    ph3 = cis(-L * eng.lambda[3])
    for i = 1:3, j = 1:3
        T0[i, j] = (i == 1 ? ph1 : i == 2 ? ph2 : ph3) * conj(Q[j, i])
    end
    for i = 1:3, j = 1:3
        acc = zero(ComplexF64)
        for k = 1:3
            acc += Q[i, k] * T0[k, j]
        end
        S[i, j] = acc
    end
    return 0
end

function filtered_probability_matrix_cd!(P, eng::StandardEngine, E, L, V, sigma, cp_sign)
    status = diagonalize!(eng, E, V, cp_sign)
    status != 0 && return status
    Q = eng.Q
    T0 = eng.T0
    lam = eng.lambda
    fill!(T0, 0)
    for i = 1:3, j = (i + 1):3
        phase = -L * (lam[i] - lam[j])
        filter_factor = exp(-0.5 * sqr(phase * sigma) / sqr(1.0e-9 * E))
        T0[i, j] = filter_factor * cis(phase)
    end
    for k = 1:3, l = 1:3
        acc = 0.0
        for i = 1:3
            for j = (i + 1):3
                acc += 2 * real(Q[k, j] * conj(Q[l, j]) * conj(Q[k, i]) * Q[l, i] * T0[i, j])
            end
            acc += sqrabs(Q[k, i]) * sqrabs(Q[l, i])
        end
        P[k, l] = acc
    end
    return 0
end

km_to_ev(x) = x / GLB_EV_TO_KM_FACTOR

"""
    probability_matrix!(P, eng, cp_sign, E, psteps, length, density, filter_sigma)

Fill the 3x3 matrix `P` with `P[α,β] = P(ν_α → ν_β)` for energy `E` (GeV)
through `psteps` layers of the given lengths (km) and densities (g/cm³).
`filter_sigma > 0` activates the low-pass filter (only for one layer).
"""
function probability_matrix!(P, eng::StandardEngine, cp_sign, E, psteps, len, dens, filter_sigma)
    E *= 1.0e9
    if filter_sigma > 0.0
        psteps == 1 || return -1
        return filtered_probability_matrix_cd!(
            P,
            eng,
            E,
            km_to_ev(len[1]),
            dens[1] * GLB_V_FACTOR * GLB_Ne_MANTLE,
            filter_sigma,
            cp_sign,
        )
    end
    if psteps > 1
        S1 = eng.S1
        fill!(S1, 0)
        S1[1, 1] = S1[2, 2] = S1[3, 3] = 1
        for i = 1:psteps
            status = S_matrix_cd!(
                eng,
                E,
                km_to_ev(len[i]),
                dens[i] * GLB_V_FACTOR * GLB_Ne_MANTLE,
                cp_sign,
            )
            status != 0 && return status
            mul!(eng.T0, eng.S, S1)
            S1 .= eng.T0
        end
        eng.S .= S1
    else
        status =
            S_matrix_cd!(eng, E, km_to_ev(len[1]), dens[1] * GLB_V_FACTOR * GLB_Ne_MANTLE, cp_sign)
        status != 0 && return status
    end
    S = eng.S
    for i = 1:3, j = 1:3
        P[j, i] = sqrabs(S[i, j])
    end
    return 0
end
