#!/usr/bin/env python3
"""R7 (docs/15_R7.md): min traces of indecomposables in totally real cubic fields (and orders).

For an order O of a totally real field and alpha in O+ put
    m(alpha) := min{Tr(delta alpha) : delta in O^vee, delta >> 0} = min{Tr gamma : gamma in alpha O^vee, gamma >> 0},
    t(O)     := max of m(alpha) over the indecomposables of O        (t = 1 for every real quadratic O_K).
Exact method: gamma >> 0 with Tr gamma <= t has Tr(gamma^2) <= t^2, so enumerate the integral lattice
(alpha O^vee, Tr(xy)) (scaled to an integer Gram matrix) by qfminim, test total positivity exactly
(a real-rooted char poly has only positive roots iff its coefficients alternate), take the least trace,
and raise t until there is a hit.
Indecomposables (mod totally positive units) and facets come from the sails package
(github.com/pav10/sails, used read-only; SAILS_PATH, default ../sails).
Usage: python3 scripts/r7_cubic.py fields Dmin Dmax         (all cubic fields with Dmin <= disc <= Dmax, PARI nflist)
       python3 scripts/r7_cubic.py ennola amin amax         (x^3 + (a-1)x^2 - ax - 1: O_K and Z[rho])
       python3 scripts/r7_cubic.py simplest amin amax       (x^3 - ax^2 - (a+3)x - 1)
       python3 scripts/r7_cubic.py crosscheck 0 0           (qfminim min traces vs an LLL-box enumeration)
       python3 scripts/r7_cubic.py thomas amin amax [sails_amax] [box_amax]
                                   (x^3 - (a+b)x^2 + abx - 1: least b >= 2a-1 with Z[rho] = O_K, m(alpha) = a)"""
import os, sys, json, time
from fractions import Fraction as F
from math import lcm
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.environ.get("SAILS_PATH", os.path.join(HERE, "..", "..", "sails")))
import cypari2
pari = cypari2.Pari()
pari.allocatemem(2 * 10 ** 9, silent=True)
from sails.compute import compute_field


def _fr(g):
    return F(int(pari.numerator(g)), int(pari.denominator(g)))


class Order:
    """The order with Z-basis zk (power basis polynomials) of Q[x]/(pol); its codifferent O^vee."""

    def __init__(self, pol, zk):
        self.pol, self.n = pol, len(zk)
        self.zk = [pari.Mod(b, pol) for b in zk]
        T = pari.matrix(self.n, self.n, [pari.trace(a * b) for a in self.zk for b in self.zk])
        Ti = T ** -1                                              # dual basis w*_j = sum_i Ti[i, j] zk_i
        self.dual = [sum(Ti[i, j] * self.zk[i] for i in range(self.n)) for j in range(self.n)]

    def elt(self, x):
        return sum(int(c) * b for c, b in zip(x, self.zk))

    def min_trace(self, alpha, t0=3, tmax=200):
        """(m(alpha), a minimizing delta in O^vee) for alpha >> 0."""
        g = [alpha * w for w in self.dual]
        G = [[_fr(pari.trace(a * b)) for b in g] for a in g]
        tr = [_fr(pari.trace(a)) for a in g]
        den = lcm(*[x.denominator for r in G for x in r])
        Gi = pari.matrix(self.n, self.n, [int(x * den) for r in G for x in r])
        U = pari.qflllgram(Gi)                       # LLL first: qfminim needs a reduced Gram matrix
        Gr = U.mattranspose() * Gi * U
        t = t0
        while t <= tmax:
            V = U * _qfminim(Gr, den, t * t)
            best = None
            for k in range(len(V)):
                v = [int(V[i, k]) for i in range(self.n)]
                s = sum(c * x for c, x in zip(v, tr))
                if s < 0:
                    v, s = [-c for c in v], -s
                if s <= 0 or s > t or (best is not None and s >= best[0]):
                    continue
                gam = sum(c * b for c, b in zip(v, g))
                cp = pari.charpoly(gam)
                if all((-1) ** (self.n - i) * cp[i] > 0 for i in range(self.n)):   # roots all > 0
                    best = (s, sum(c * w for c, w in zip(v, self.dual)))
            if best is not None:
                assert best[0] == int(best[0])
                return int(best[0]), best[1]
            t *= 2
        raise RuntimeError("min trace above tmax")


def _qfminim(Gr, den, bound):
    """All v with v^T Gr v <= den*bound (one of each +-v).  Exact integer Fincke-Pohst when PARI can do it,
    otherwise PARI's floating-point version at 200 digits on Gr/den with a slack (callers re-check exactly)."""
    try:
        return pari.qfminim(Gr, den * bound)[2]
    except cypari2.PariError:
        old = pari.default("realprecision")
        pari.set_real_precision(200)
        try:
            R = Gr * pari("1.") / den
            return pari.qfminim(R, bound + pari("1e-20"), flag=2)[2]
        finally:
            pari.set_real_precision(int(old))


def field_data(coeffs, order_zk=None, label=None):
    """sails (totally positive sign class) + min traces.  order_zk: basis of a suborder in
    O_K-coordinates (sails convention), None for O_K."""
    t0 = time.time()
    res = compute_field(coeffs, label=label, tp_only=True, order=order_zk)
    K = res.K
    O = Order(K.pol, K.zk)
    orth = res.orthants[0]
    ms = []
    for N, x in orth.indecomposables:
        m, _ = O.min_trace(O.elt(x))
        ms.append((int(N), m, list(map(int, x))))
    rec = res.record()["orthants"][0]
    d1 = [len(Fc["lattice_points"]) for Fc in rec["facets"] if Fc["d"] == 1]
    hist = {}
    for _, m, _ in ms:
        hist[m] = hist.get(m, 0) + 1
    worst = max(ms, key=lambda r: (r[1], r[0]))
    return dict(label=label, coeffs=list(coeffs), disc=int(K.disc), field_disc=int(K.field_disc),
                index=int(K.order_index), h=int(K.class_number), iota=len(ms), t=worst[1],
                hist={str(k): v for k, v in sorted(hist.items())}, worst=worst,
                facets=len(rec["facets"]), max_d=max(Fc["d"] for Fc in rec["facets"]),
                facets_d1=len(d1), mu1=max(d1, default=0), sec=round(time.time() - t0, 2))


def ennola(a):
    return [-1, -a, a - 1, 1]                       # x^3 + (a-1)x^2 - a x - 1, ascending


def simplest(a):
    return [-1, -(a + 3), -a, 1]                    # x^3 - a x^2 - (a+3) x - 1


def thomas(a, b):
    return [-1, a * b, -(a + b), 1]                 # x(x-a)(x-b) - 1 = x^3 - (a+b)x^2 + ab x - 1


def box_indecomposable(O, alpha):
    """Independent test: alpha >> 0 is decomposable iff some beta in O has 0 << beta << alpha; such beta has
    Tr(beta^2) < Tr(alpha^2), so enumerate those (qfminim on the trace form of O) and test the signs
    (numerically, exact char poly check near the boundary).  Returns (indecomposable?, #vectors scanned)."""
    import numpy as np
    n = O.n
    T = pari.matrix(n, n, [int(pari.trace(u * v)) for u in O.zk for v in O.zk])
    B = int(pari.trace(alpha * alpha)) - 1
    U = pari.qflllgram(T)
    V = _qfminim(U.mattranspose() * T * U, 1, B)
    if len(V):
        V = U * V
    if len(V) == 0:
        return True, 0
    roots = [float(pari.real(r)) for r in pari.polroots(O.pol)]   # totally real
    E = np.array([[float(pari.subst(pari.lift(b), "x", r)) for b in O.zk] for r in roots])
    al = np.array([float(pari.subst(pari.lift(alpha), "x", r)) for r in roots])
    M = np.array([[int(V[i, k]) for i in range(n)] for k in range(len(V))], dtype=float)
    C = M @ E.T                                      # conjugates of the found beta (one of +-beta)
    for sgn in (1, -1):
        X = sgn * C
        near = np.all(X > -1e-6, axis=1) & np.all(X < al + 1e-6, axis=1)
        for k in np.nonzero(near)[0]:
            beta = sgn * sum(int(M[k, i]) * O.zk[i] for i in range(n))
            if all(_tp(z) for z in (beta, alpha - beta)):
                return False, len(V)
    return True, len(V)


def _tp(z):
    cp = pari.charpoly(z)
    n = int(pari.poldegree(cp))
    return all((-1) ** (n - i) * cp[i] > 0 for i in range(n))


def thomas_cert(a, bmax_extra=200, sails_max=0, box_max=0):
    """Theorem 15.6 input: the least b >= 2a-1 with Z[rho] = O_K, and the checks at that b:
    m(alpha) = a for alpha = 1 + (a-1)(b - rho) (exact), optional box test and sails t(K)."""
    t0 = time.time()
    for b in range(2 * a - 1, 2 * a + bmax_extra):
        c = thomas(a, b)
        pol = pari.Pol(list(reversed(c)))
        if not pari.polisirreducible(pol):
            continue
        nf = pari.nfinit(pol)
        if abs(int(nf[3])) == 1:
            break
    else:
        return dict(a=a, b=None)
    O = Order(pol, list(pari.nfbasis(pol)))          # an integral basis of O_K (= Z[rho] here)
    rho = pari.Mod(pari("x"), pol)
    alpha = 1 + (a - 1) * (b - rho)
    assert _tp(alpha)
    m, delta = O.min_trace(alpha, t0=a + 1)
    r = dict(a=a, b=b, b_skipped=b - (2 * a - 1), disc=int(nf[2]), m_alpha=m, delta_fprime=str(pari.lift(delta * pari.deriv(pol))))
    if a <= box_max:
        r["box_indecomposable"], r["box_vectors"] = box_indecomposable(O, alpha)
    if a <= sails_max:
        fd = field_data(c, label=f"thomas{a},{b}")
        r.update(t_sails=fd["t"], iota=fd["iota"], max_d=fd["max_d"], hist=fd["hist"])
    r["sec"] = round(time.time() - t0, 2)
    return r


def min_trace_box(O, alpha, m):
    """Independent check of m(alpha) <= m (no qfminim): gamma = delta*alpha ranges over the lattice
    alpha O^vee; take an LLL basis g_j (embedding space) and its dual basis g*_j; gamma >> 0 with Tr gamma <= m
    has |gamma| <= m, so its coordinates satisfy |c_j| <= m |g*_j|.  Enumerate that box (numpy) and return the
    least Tr gamma over gamma >> 0 found (float)."""
    import numpy as np
    roots = [float(pari.real(r)) for r in pari.polroots(O.pol)]
    ev = lambda z: [float(pari.subst(pari.lift(z), "x", r)) for r in roots]
    G = np.array([ev(alpha * w) for w in O.dual])
    U = pari.qflllgram(pari.matrix(O.n, O.n, [float(t) for t in (G @ G.T).ravel()]))
    U = np.array([[int(U[i, j]) for j in range(O.n)] for i in range(O.n)])
    Gl = U.T @ G
    Gd = np.linalg.inv(Gl).T
    B = [int(m * np.linalg.norm(Gd[j])) + 1 for j in range(O.n)]
    assert O.n == 3
    C1, C2 = np.meshgrid(np.arange(-B[1], B[1] + 1), np.arange(-B[2], B[2] + 1), indexing="ij")
    C1, C2 = C1.ravel(), C2.ravel()
    best = float("inf")
    for c0 in range(-B[0], B[0] + 1):
        g = c0 * Gl[0][None, :] + C1[:, None] * Gl[1][None, :] + C2[:, None] * Gl[2][None, :]
        ok = np.all(g > 1e-9, axis=1)
        if ok.any():
            best = min(best, float(g[ok].sum(axis=1).min()))
    return best


def crosscheck(coeffs, nmax=4):
    """m(alpha) (qfminim) vs min_trace_box for the nmax indecomposables of largest norm."""
    res = compute_field(coeffs, tp_only=True)
    O = Order(res.K.pol, res.K.zk)
    out = []
    for N, x in sorted(res.orthants[0].indecomposables, key=lambda t: -t[0])[:nmax]:
        al = O.elt(x)
        m, _ = O.min_trace(al)
        b = min_trace_box(O, al, m)
        assert abs(b - m) < 1e-4, (coeffs, x, m, b)
        out.append((int(N), m))
    return out


def zrho_basis(coeffs):
    nf = pari.nfinit(pari.Pol(list(reversed(coeffs))))
    x = pari("x")
    idx = abs(int(nf[3]))                            # nf.index = [O_K : Z[rho]]
    return [[int(c) for c in pari.nfalgtobasis(nf, x ** j)] for j in range(3)], idx


if __name__ == "__main__":
    mode = sys.argv[1]
    lo, hi = int(sys.argv[2]), int(sys.argv[3])
    if mode == "fields":
        pols = list(pari(f'nflist("C3",[{lo},{hi}])')) + list(pari(f'nflist("S3",[{lo},{hi}],0)'))
        pols.sort(key=lambda p: abs(int(pari.nfdisc(p))))
        for p in pols:
            coeffs = [int(c) for c in pari.Vecrev(p)]
            print(json.dumps(field_data(coeffs, label=str(pari.nfdisc(p)))), flush=True)
    elif mode == "crosscheck":                     # python3 scripts/r7_cubic.py crosscheck 0 0
        for c in ([-2, -4, 2, 1], [1, -8, -3, 1], [1, -12, -6, 1], [1, -2, -1, 1], thomas(4, 7), thomas(5, 9)):
            print(c, crosscheck(c), flush=True)
    elif mode == "thomas":
        smax = int(sys.argv[4]) if len(sys.argv) > 4 else 0
        bmax = int(sys.argv[5]) if len(sys.argv) > 5 else 0
        for a in range(lo, hi + 1):
            print(json.dumps(thomas_cert(a, sails_max=smax, box_max=bmax)), flush=True)
    elif mode in ("ennola", "simplest"):
        fam = ennola if mode == "ennola" else simplest
        for a in range(lo, hi + 1):
            c = fam(a)
            r = field_data(c, label=f"{mode}{a}")
            r["a"] = a
            zb, idx = zrho_basis(c)
            r["zrho_index"] = idx
            if idx > 1:
                rz = field_data(c, order_zk=zb, label=f"{mode}{a}-Zrho")
                r["zrho"] = {k: rz[k] for k in ("disc", "iota", "t", "hist", "worst", "max_d")}
            print(json.dumps(r), flush=True)
