"""
zak_kernel.py -- per-k dressed-band Zak phase of the cavity-coupled SSH chain.

Self-contained (numpy only). All chain / cavity parameters are passed explicitly,
so every routine is reentrant and safe under multiprocessing / job arrays.

Physical convention: quasiparticle frequencies are 2 x (positive eigenvalues of
tau3 H_C)  (bare transmon at 1, matching PHOTONICChain._nambu_prefactor = 2).

Reference: notebooks/unitary/zak_phase/K_SPACE_CALCULATIONS.md
"""
import numpy as np

_TWO_PI = 2.0 * np.pi


# --------------------------------------------------------------------------
#  H_C / H_J sector blocks  (generic {k,-k}, 8x8, honest 1/2 Psi^dag H Psi)
#  basis (b_kA, b_kB, b_-kA, b_-kB) then daggers
# --------------------------------------------------------------------------
def sector_HC_HJ(k, gc1, gc2, gJ1, gJ2):
    gJ = gJ1 + gJ2
    dC = gc1 + gc2 * np.exp(-1j * k)
    dJ = gJ1 + gJ2 * np.exp(-1j * k)
    dCc, dJc = np.conj(dC), np.conj(dJ)

    A = np.zeros((4, 4), dtype=complex)
    A[0, 0] = A[1, 1] = A[2, 2] = A[3, 3] = 0.5
    A[0, 1] = 0.5 * dC;  A[1, 0] = 0.5 * dCc
    A[2, 3] = 0.5 * dCc; A[3, 2] = 0.5 * dC
    B = np.zeros((4, 4), dtype=complex)
    B[0, 3] = -0.5 * dC;  B[1, 2] = -0.5 * dCc
    B[2, 1] = -0.5 * dCc; B[3, 0] = -0.5 * dC

    AJ = np.zeros((4, 4), dtype=complex)
    AJ[0, 0] = AJ[1, 1] = AJ[2, 2] = AJ[3, 3] = 0.5 * gJ
    AJ[0, 1] = -0.5 * dJ;  AJ[1, 0] = -0.5 * dJc
    AJ[2, 3] = -0.5 * dJc; AJ[3, 2] = -0.5 * dJ
    BJ = np.zeros((4, 4), dtype=complex)
    BJ[0, 2] = BJ[1, 3] = BJ[2, 0] = BJ[3, 1] = -0.5 * gJ
    BJ[0, 3] = -0.5 * dJ;  BJ[1, 2] = -0.5 * dJc
    BJ[2, 1] = -0.5 * dJc; BJ[3, 0] = -0.5 * dJ

    HC = np.zeros((8, 8), dtype=complex)
    HC[:4, :4] = A;  HC[:4, 4:] = B;  HC[4:, :4] = B.conj();  HC[4:, 4:] = A.conj()
    HJ = np.zeros((8, 8), dtype=complex)
    HJ[:4, :4] = AJ; HJ[:4, 4:] = BJ; HJ[4:, :4] = BJ.conj(); HJ[4:, 4:] = AJ.conj()
    return HC, HJ


def diagonalize_C_sector(HC, k, gc1, gc2):
    """(T, T^dag HC T): T paraunitary, result numerically diagonal. Needs g_c1+g_c2 < 1/2."""
    d = gc1 + gc2 * np.exp(-1j * k)
    ad, phi = abs(d), np.angle(d)
    za, zb = 0.25 * np.log(1 + 2 * ad), 0.25 * np.log(1 - 2 * ad)
    ep, em = np.exp(-1j * phi / 2), np.exp(1j * phi / 2)
    r2 = 1.0 / np.sqrt(2.0)

    # phased sublattice rotation, rows (alpha,beta) at +k ; (alpha,beta) at -k has phi -> -phi
    Mrot = np.zeros((4, 4), dtype=complex)
    Mrot[0, 0] = ep * r2;  Mrot[0, 1] = em * r2
    Mrot[1, 0] = ep * r2;  Mrot[1, 1] = -em * r2
    Mrot[2, 2] = em * r2;  Mrot[2, 3] = ep * r2
    Mrot[3, 2] = em * r2;  Mrot[3, 3] = -ep * r2
    Minv = np.linalg.inv(Mrot)
    Trot = np.zeros((8, 8), dtype=complex)
    Trot[:4, :4] = Minv;  Trot[4:, 4:] = Minv.conj()

    # two-mode squeeze on (alpha_k, alpha_-k) = (0,2) and (beta_k, beta_-k) = (1,3)
    ca, sa, cb, sb = np.cosh(za), np.sinh(za), np.cosh(zb), np.sinh(zb)
    Tsq = np.eye(8, dtype=complex)
    Tsq[0, 0] = Tsq[2, 2] = ca;  Tsq[1, 1] = Tsq[3, 3] = cb
    Tsq[4, 4] = Tsq[6, 6] = ca;  Tsq[5, 5] = Tsq[7, 7] = cb
    Tsq[0, 6] = Tsq[2, 4] = sa;  Tsq[1, 7] = Tsq[3, 5] = sb
    Tsq[6, 0] = Tsq[4, 2] = sa;  Tsq[7, 1] = Tsq[5, 3] = sb

    T = Trot @ Tsq
    return T, T.conj().T @ HC @ T


# --------------------------------------------------------------------------
#  Bosonic-Fock data for the 3-QP virtual sector  (cached per Nq)
# --------------------------------------------------------------------------
_FOCK = {}


def _get_fock(Nq):
    if Nq in _FOCK:
        return _FOCK[Nq]
    m = 4
    adg = np.diag(np.sqrt(np.arange(1, Nq)), 1).T
    ident = np.eye(Nq)

    def op(idx):
        o = np.array([[1.0]])
        for j in range(m):
            o = np.kron(o, adg if j == idx else ident)
        return o

    gd = [op(i) for i in range(m)]
    states = np.array(list(np.ndindex(*([Nq] * m))))
    nexc = states.sum(1)
    Xi = np.where(nexc == 1)[0]
    Xi = Xi[np.argsort([int(np.argmax(states[x])) for x in Xi])]   # alpha_k, beta_k, alpha_-k, beta_-k
    Ri = np.where(nexc == 3)[0]
    # sub-tensor  G[i,j, r, X] = (gd[i] @ gd[j])[Ri, Xi]     shape (4,4, |Ri|, 4)
    G = np.empty((m, m, len(Ri), len(Xi)), dtype=float)
    for i in range(m):
        for j in range(m):
            G[i, j] = (gd[i] @ gd[j])[np.ix_(Ri, Xi)]
    _FOCK[Nq] = (states.astype(float), G, Xi, Ri)
    return _FOCK[Nq]


# --------------------------------------------------------------------------
#  Schrieffer-Wolff correction  (verbatim schrieffer_wolff_correction)
#  -> H_SW on (cavity Fock d_ph) x (alpha_k, beta_k, alpha_-k, beta_-k)
# --------------------------------------------------------------------------
def H_SW_composite(HJ_d, wq, omega_r, gamma, d_ph=4, Nq=4, eps=1e-6):
    m = d_X = 4
    BJ = 2.0 * HJ_d[:m, m:]                                   # pairing block ~B_J
    states, G, Xi, Ri = _get_fock(Nq)
    Vr = np.tensordot(BJ, G, axes=([0, 1], [0, 1]))           # <r|B_J|X>, (|Ri|, 4)
    E = states @ wq
    E_X, E_r = E[Xi], E[Ri]

    photons = np.arange(1, d_ph)
    s_off = gamma * np.sqrt(photons)                          # s[n, n+1]
    dEXr = E_X[:, None] - E_r[None, :]                        # (4, |Ri|)

    H = np.zeros((d_ph * d_X, d_ph * d_X), dtype=complex)
    for n1 in range(d_ph):
        for n2 in range(d_ph):
            for nm in range(d_ph):
                sf = 0.0
                if abs(n1 - nm) == 1 and abs(nm - n2) == 1:
                    sf = s_off[min(n1, nm)] * s_off[min(nm, n2)]
                if sf == 0.0:
                    continue
                d1 = (n1 - nm) * omega_r + dEXr
                d2 = (nm - n2) * omega_r - dEXr.T
                ok1, ok2 = np.abs(d1) > eps, np.abs(d2) > eps
                i1 = np.where(ok1, 1.0 / np.where(ok1, d1, 1.0), 0.0)
                i2 = np.where(ok2, 1.0 / np.where(ok2, d2, 1.0), 0.0)
                blk = sf * ((Vr.T * i1) @ Vr.conj() - Vr.T @ (Vr.conj() * i2))
                H[n1 * d_X:(n1 + 1) * d_X, n2 * d_X:(n2 + 1) * d_X] += blk
    H = 0.5 * (H + H.conj().T)
    return 0.5 * H


# --------------------------------------------------------------------------
#  Effective per-k Hamiltonian on {alpha_k, beta_k} x cavity Fock
# --------------------------------------------------------------------------
def ssh_bloch_U(k, gc1, gc2):
    """[u_alpha | u_beta]: eigenvectors of h(k), eigenvalues +|delta_C_k|, -|delta_C_k|."""
    d = gc1 + gc2 * np.exp(-1j * k)
    _, v = np.linalg.eigh(np.array([[0.0 + 0j, d], [np.conj(d), 0.0]]))
    return v[:, [1, 0]]


def H_eff(k, gc1, gc2, gJ1, gJ2, omega_r, gamma, d_ph, with_sw):
    """(2 d_ph)x(2 d_ph) Hermitian; ordering index = n*2 + s, s: 0=alpha, 1=beta."""
    HC, HJ = sector_HC_HJ(k, gc1, gc2, gJ1, gJ2)
    T, HC_d = diagonalize_C_sector(HC, k, gc1, gc2)
    HJ_d = T.conj().T @ HJ @ T
    wq = 2.0 * np.diag(HC_d[:4, :4]).real
    w = wq[:2]
    M = HJ_d[:2, :2] + HJ_d[:2, :2].conj().T                   # A_J + D_J^T, +k band block

    a = np.diag(np.sqrt(np.arange(1, d_ph)), 1)
    diag = np.array([omega_r * n + w[s] for n in range(d_ph) for s in (0, 1)], dtype=complex)
    H = np.diag(diag) - gamma * np.kron(a + a.conj().T, M)
    if with_sw:
        Hsw = H_SW_composite(HJ_d, wq, omega_r, gamma, d_ph=d_ph)
        keep = [n * 4 + s for n in range(d_ph) for s in (0, 1)]
        H = H + Hsw[np.ix_(keep, keep)]
    return H


# --------------------------------------------------------------------------
#  Wilson-loop Zak phase of the dressed bands, at one parameter point
# --------------------------------------------------------------------------
def zak_point(gc1, gc2, gJ1, gJ2, omega_r, gamma, d_ph, with_sw,
              Nk=96, s=(-0.25, 0.25)):
    """Returns (zak_b0, zak_b1, zak_comb12, gap01, gap12).
    Composite state lifted to (sublattice x cavity): Psi_n = U(k) C_n(k),
    flattened sublattice-outer; the BZ-boundary factor acts on the sublattice half."""
    ks = np.linspace(0.0, _TWO_PI, Nk, endpoint=False)
    nb = 3
    Psi = np.empty((Nk, 2 * d_ph, nb), dtype=complex)
    E = np.empty((Nk, nb))
    for j, k in enumerate(ks):
        H = H_eff(k, gc1, gc2, gJ1, gJ2, omega_r, gamma, d_ph, with_sw)
        w, V = np.linalg.eigh(H)
        U = ssh_bloch_U(k, gc1, gc2)
        E[j] = w[:nb]
        for col in range(nb):
            C = V[:, col].reshape(d_ph, 2).T
            Psi[j, :, col] = (U @ C).reshape(-1)

    gph = np.repeat(np.exp(-1j * _TWO_PI * np.asarray(s)), d_ph)

    def wilson1(col):
        pr = 1.0 + 0j
        for j in range(Nk - 1):
            pr *= np.vdot(Psi[j, :, col], Psi[j + 1, :, col])
        pr *= np.vdot(Psi[Nk - 1, :, col], gph * Psi[0, :, col])
        return float(-np.angle(pr))

    W = np.eye(2, dtype=complex)
    for j in range(Nk):
        B = (gph[:, None] * Psi[0, :, 1:3]) if j == Nk - 1 else Psi[(j + 1) % Nk, :, 1:3]
        W = W @ (Psi[j, :, 1:3].conj().T @ B)
    zak_comb = float(-np.angle(np.linalg.det(W)))

    return (wilson1(0), wilson1(1), zak_comb,
            float((E[:, 1] - E[:, 0]).min()), float((E[:, 2] - E[:, 1]).min()))
