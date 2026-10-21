#!/usr/bin/env python3
"""R3 certificate without class groups of E (docs/11_R3.md, Corollary 11.5): for a prime p = 3 (mod 4)
and a form class f of discriminant -8p in the principal genus with nu(f) = nu*(p) (scripts/r3_tools.py),
take g with g^2 = f, an ideal b = [a, -b/2 + sqrt(-2p)] of k2 = Q(sqrt -2p) in the class g, and the
ideal lattice  M = (b O_E, N_{E/K}/a),  E = K(sqrt -2) = Q(sqrt p, sqrt -2):  Q(x + y sqrt -2) =
(x^2 + 2y^2)/a on K^2.  Theorems 10.2(c) and 11.3 give s_sq(M) = n(j(g)) = nu(f).  Vectors are found
only on the edges carrying a represented class (Fincke-Pohst for Tr(delta Q(v)) = 1); each value gamma
comes with its edge functional delta (totally positive, in the codifferent, Tr(delta gamma) = 1), which
certifies indecomposability by itself (gamma = b + c with b, c >> 0 would give Tr >= 2).  The file is
re-checked by tests/test_r3.py::test_stored_form_certificates.  Writes data/certificates/D{p}_form.json.
Usage: python3 scripts/certify_form.py p [a b c]   (default: the best form found by r3_tools)"""
import os, sys, json, time
if hasattr(sys, "set_int_max_str_digits"):
    sys.set_int_max_str_digits(0)          # certificates hold units with > 4300 digits
from fractions import Fraction as F
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "src"))
sys.path.insert(0, HERE)
from rqf import QF
from lattice import hnf
from rqf import same_square_class
from lattice import Bform, enum_short_reduced, enum_short, qnorm


def short1(T):
    """All +-pairs v with v^T T v = 1 (T integer, positive definite).  The reduction is only a search
    aid -- every vector found is re-checked exactly -- so PARI's qflllgram is used when available
    (exact pure-Python LLL is slow for the 10^3-digit entries of trace forms of deep edges)."""
    try:
        import cypari2
    except ImportError:
        return [v for _, v in enum_short_reduced(T, 1)]
    P = cypari2.Pari()
    n = len(T)
    M = P.matrix(n, n, [x for r in T for x in r])
    Um = P.qflllgram(M)                                        # columns = new basis
    U = [[int(Um[i, j]) for i in range(n)] for j in range(n)]  # rows = new basis
    Tr = [[sum(U[a][i] * T[i][k] * U[b][k] for i in range(n) for k in range(n)) for b in range(n)]
          for a in range(n)]
    out = []
    for nv, w in enum_short(Tr, 1):
        v = tuple(sum(w[i] * U[i][j] for i in range(n)) for j in range(n))
        assert qnorm(T, v) == nv
        out.append(v)
    return out
from r3_tools import (Sieve, field, norm_classes, qf_compose, qf_pow, h8, _is_prime,
                      convergents_big, cf_sqrt)


def forms_of_disc(disc):
    from math import isqrt, gcd
    out = []
    for a in range(1, isqrt(-disc // 3) + 1):
        for b in range(-a + 1, a + 1):
            if (b * b - disc) % (4 * a):
                continue
            c = (b * b - disc) // (4 * a)
            if c >= a and not (c == a and b < 0) and gcd(gcd(a, abs(b)), c) == 1:
                out.append((a, b, c))
    return out


def square_root_class(f, disc, h):
    """g in Cl(disc) with g^2 = f (None if f is not a square)."""
    if (h // 2) % 2 == 1:                   # Cl^2 has odd order h/2 and contains f
        g = qf_pow(f, (h // 2 + 1) // 2, disc)
        assert qf_compose(g, g) == f
        return g
    for g in forms_of_disc(disc):
        if qf_compose(g, g) == f:
            return g
    return None


def ideal_lattice(p, g):
    """Z-basis (K^2 vectors) and Gram of (b O_E, N/a), b = [a, -b/2 + sqrt(-2p)] for g = (a, b, c)."""
    a, b, _ = g
    D = p
    # E-elements as pairs (x, y) in K^2 meaning x + y sqrt(-2); sqrt(-2p) = sqrt p * sqrt -2
    gens_b = [(QF(D, a), QF(D, 0)), (QF(D, F(-b, 2)), QF(D, 0, 1))]
    OE = [(QF(D, 1), QF(D, 0)), (QF(D, 0, 1), QF(D, 0)), (QF(D, 0), QF(D, 1)),
          (QF(D, 0), QF(D, F(1, 2), F(1, 2)))]
    mul = lambda u, v: (u[0] * v[0] - 2 * u[1] * v[1], u[0] * v[1] + u[1] * v[0])
    coords = lambda z: [2 * z[0].a, 2 * z[0].b, 2 * z[1].a, 2 * z[1].b]       # in (1, sqrt p, sqrt-2, sqrt-2p)/2
    rows = [[int(t) for t in coords(mul(x, y))] for x in gens_b for y in OE]
    H = hnf(rows)
    assert len(H) == 4
    basis = [(QF(D, F(r[0], 2), F(r[1], 2)), QF(D, F(r[2], 2), F(r[3], 2))) for r in H]
    G = [[QF(D, F(1, a)), QF(D, 0)], [QF(D, 0), QF(D, F(2, a))]]
    return G, basis


def edge_delta(D, conv, j):
    """delta with Tr(delta alpha_{j-1}) = 1, Tr(delta alpha_j) = 0: the functional of the edge
    {alpha_{j-1} + t alpha_j} (lattice.edge_deltas, same convention)."""
    p, q = conv
    A = [[2 * p[j], 2 * q[j] * D], [2 * p[j + 1], 2 * q[j + 1] * D]]
    det = A[0][0] * A[1][1] - A[0][1] * A[1][0]
    return QF(D, F(A[1][1], det), F(-A[1][0], det))


def class_forms(p, rec, sv):
    """Per square class (field(..., need_reps=True)): {form: member (j, t, twice) whose norm it represents}."""
    disc, cache, out = -8 * p, {}, []
    for (xy, tw0, N0, inE), mem in zip(rec["reps"], rec["members"]):
        d = {}
        if inE:
            for (j, t, tw, n) in mem:
                n <<= 2 * tw
                for g in norm_classes(n, sv.factor(n), disc, cache):
                    d.setdefault(g, (j, t, tw))
        out.append(d)
    return out


def main(p, f=None):
    assert p % 4 == 3 and _is_prime(p)
    t0 = time.time()
    sv = Sieve(4 * p)
    rec = field(p, sv, need_reps=True)
    cfs = class_forms(p, rec, sv)
    nu = {}
    for d in cfs:
        for g in d:
            nu[g] = nu.get(g, 0) + 1
    if f is None:
        f = max(nu, key=lambda g: (nu[g], g))
    disc, h = -8 * p, h8(p)
    g = square_root_class(f, disc, h)
    if g is None:
        raise SystemExit(f"form {f} is not a square in Cl({disc}); pick another best form")
    print(f"  [{time.time() - t0:.1f}s] field, nu, g done", file=sys.stderr, flush=True)
    G, basis = ideal_lattice(p, g)
    for u in basis:
        for v in basis:
            assert Bform(G, u, v).is_integral()
    cf = cf_sqrt(p)
    conv = convergents_big(p, cf)
    s = len(cf[0]) - 1
    eps = QF(p, conv[0][s], conv[1][s])                       # alpha_{s-1} = fundamental unit
    assert eps.norm() == 1
    # mu with mu^2 = 2 eps (H), from the convergent of norm +-2 (Q_{k+1} = 2); theta = (mu/2) sqrt-2
    k2 = next(k for k in range(s) if cf[2][k + 1] == 2)
    mu = QF(p, conv[0][k2 + 1], conv[1][k2 + 1])
    while mu * mu != 2 * eps:                                  # mu^2 / 2 = eps^j, j odd
        mu = mu * eps.conj()
        assert (mu * mu).emb1_gt(0)
    theta = lambda v: (-mu * v[1], mu * v[0] / 2)              # Q(theta v) = eps Q(v), theta M = M
    Bm = [[F(c) for c in (b[0].a, b[0].b, b[1].a, b[1].b)] for b in basis]

    def coords(v):
        """Integer coordinates of v (K^2) in the Z-basis."""
        t = [F(c) for c in (v[0].a, v[0].b, v[1].a, v[1].b)]
        A = [[Bm[i][j] for i in range(4)] + [t[j]] for j in range(4)]     # A c = t
        for c in range(4):
            r = next(r for r in range(c, 4) if A[r][c] != 0)
            A[c], A[r] = A[r], A[c]
            A[c] = [x / A[c][c] for x in A[c]]
            for r in range(4):
                if r != c and A[r][c] != 0:
                    A[r] = [x - A[r][c] * y for x, y in zip(A[r], A[c])]
        out = [A[r][4] for r in range(4)]
        assert all(x.denominator == 1 for x in out)
        return [int(x) for x in out]

    targets = {}                                               # j -> [(class id, gamma, tw)]
    for cid, d in enumerate(cfs):
        if f in d:
            j, t, tw = d[f]
            x, y = conv[0][j] + t * conv[0][j + 1], conv[1][j] + t * conv[1][j + 1]
            targets.setdefault(j, []).append((cid, QF(p, x, y), tw))
    found = {}
    for j, lst in targets.items():
        dl, kk = edge_delta(p, conv, j), 0                    # balance by eps (odd powers allowed)
        while True:
            if (dl * eps).trace() < dl.trace():
                dl, kk = dl * eps, kk + 1
            elif (dl * eps.conj()).trace() < dl.trace():
                dl, kk = dl * eps.conj(), kk - 1
            else:
                break
        # values on this edge: beta * eps^{-kk}, of class [beta eps^kk]
        T = [[(dl * Bform(G, basis[a], basis[b])).trace() for b in range(4)] for a in range(4)]
        assert all(x.denominator == 1 for r in T for x in r)
        T = [[int(x) for x in r] for r in T]
        vecs = []
        for v in short1(T):
            x = sum((basis[i][0] * v[i] for i in range(4)), QF(p, 0))
            y = sum((basis[i][1] * v[i] for i in range(4)), QF(p, 0))
            vecs.append((v, (x, y), Bform(G, (x, y), (x, y))))
        for cid, beta, tw in lst:
            for v, xy, val in vecs:
                if same_square_class(val, beta * (eps if kk % 2 else QF(p, 1))):
                    if (kk - tw) % 2 == 0:                     # val is already in [beta eps^tw]
                        found[cid] = (v, val, dl)
                    else:                                      # theta v: value eps*val, witness dl/eps
                        w = theta(xy)
                        found[cid] = (coords(w), eps * val, dl * eps.conj())
                    break
    print(f"  [{time.time() - t0:.1f}s] edges done", file=sys.stderr, flush=True)
    k = len(found)
    vals = [found[c][1] for c in sorted(found)]
    for c in sorted(found):                                     # independent re-checks
        v, val, dl = found[c]
        x = sum((basis[i][0] * v[i] for i in range(4)), QF(p, 0))
        y = sum((basis[i][1] * v[i] for i in range(4)), QF(p, 0))
        assert Bform(G, (x, y), (x, y)) == val
        assert val.is_tot_pos() and dl.is_tot_pos() and (dl * val).trace() == 1
        assert dl.trace().denominator == 1 and (dl * QF(p, 0, 1)).trace().denominator == 1
    for i in range(k):
        for j in range(i):
            assert not same_square_class(vals[i], vals[j])
    q = lambda x: [str(x.a), str(x.b)]
    cert = dict(D=p, family="ideal lattice (b O_E, N/a), b in the class g of Q(sqrt -2p), g^2 = f",
                form_f=list(f), form_g=list(g), nu_f=nu[f], nu_star=max(nu.values()), h8=h, s=s,
                gram=[[q(G[i][j]) for j in range(2)] for i in range(2)], glue=[],
                zbasis=[[q(v[0]), q(v[1])] for v in basis], s_sq_lower_bound=k,
                vectors=[list(found[c][0]) for c in sorted(found)], values=[q(x) for x in vals],
                deltas=[q(found[c][2]) for c in sorted(found)])
    path = os.path.join(HERE, "..", "data", "certificates", f"D{p}_form.json")
    json.dump(cert, open(path, "w"), indent=1)
    print(f"p={p}: s={s} nu(f)={nu[f]} nu*={max(nu.values())}; CERTIFIED s_sq >= {k} "
          f"({len(targets)} edges); {time.time() - t0:.1f}s", flush=True)
    return k, nu[f]

if __name__ == "__main__":
    p = int(sys.argv[1])
    f = tuple(map(int, sys.argv[2:5])) if len(sys.argv) > 4 else None
    main(p, f)
