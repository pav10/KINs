"""
Binary O_K-lattices over real quadratic K = Q(sqrt D): realization of abstract Gram
matrices, the rank-4 Z-trace lattice, exact short-vector enumeration, and the
represented indecomposable classes / orbits.

All arithmetic exact (Fraction / int) except the Cholesky step of Fincke-Pohst,
which is only used to bound the search box (every candidate's norm is then
recomputed exactly in integers, and the bound is padded), so results are exact
provided the float bounds are not off by more than the padding (they are not for
the sizes used here; see tests).
"""
from fractions import Fraction as F
from functools import lru_cache
from math import gcd, floor, ceil, sqrt
from rqf import QF, tp_unit, same_square_class
from indec import indecomposables, normalize_unit
from gram import rank_over_K, valid_off_diagonal, omega


# ----------------------------------------------------------------------------
# Square-class representatives with small trace
# ----------------------------------------------------------------------------
def balance(alpha, eta):
    """Minimal-trace element of alpha * eta^(2k), k in Z (same square class,
    and same O_K-line: alpha*eta^2 = Q(eta v) if alpha = Q(v))."""
    e2 = eta * eta
    e2c = e2.conj()
    cur = alpha
    while True:
        for nxt in (cur * e2, cur * e2c):
            if nxt.trace() < cur.trace():
                cur = nxt
                break
        else:
            return cur


def balanced_classes(D):
    """One small-trace rep per square class of indecomposables (kappa_sq(K) many)."""
    eta = tp_unit(D)
    pool = []
    for a in indecomposables(D):
        pool += [a, a * eta]
    classes = []
    for x in pool:
        if not any(same_square_class(x, r) for r in classes):
            classes.append(x)
    return [balance(c, eta) for c in classes]


# ----------------------------------------------------------------------------
# Max rank<=2 totally PSD Gram with prescribed diagonal set (subset search)
# ----------------------------------------------------------------------------
def find_winner(D, reps, sizecap=64):
    """
    Largest subset S of `reps` and off-diagonal entries c_ij in O_K such that the
    Gram matrix (diag = S) is totally PSD of rank <= 2 over K.
    Exact acceptance: rank_K <= 2 (sufficient given tot. pos. diagonal and
    tot. nonneg. 2x2 minors). 3x3 principal-minor pruning is a sound necessary test.
    NOTE: a LOWER bound for S(K,2): reps are fixed per class; classes reached only
    via a rep alpha*xi^2 with xi a NON-unit are not explored (unit-square changes
    of rep are harmless: they rescale a basis vector by a unit).
    """
    cache = {}

    def vc(a, b):
        key = (a.a, a.b, b.a, b.b)
        if key not in cache:
            cache[key] = valid_off_diagonal(D, a, b)
        return cache[key]

    best = [None]
    two = QF(D, 2, 0)

    def search(M, start):
        k = len(M)
        if k >= 4 and rank_over_K(M) > 2:
            return
        if best[0] is None or k > len(best[0]):
            best[0] = [row[:] for row in M]
        if k >= sizecap:
            return
        for i in range(start, len(reps)):
            vnew = reps[i]
            if k == 0:
                search([[vnew]], i + 1)
                continue
            lists = [vc(M[j][j], vnew) for j in range(k)]

            def bt(idx, cur):
                if idx == k:
                    nM = [row[:] + [cur[j]] for j, row in enumerate(M)]
                    nM.append(list(cur) + [vnew])
                    search(nM, i + 1)
                    return
                for c in lists[idx]:
                    ok = True
                    for j in range(idx):
                        det = (M[j][j] * M[idx][idx] * vnew + two * M[j][idx] * cur[j] * c
                               - M[j][j] * c * c - M[idx][idx] * cur[j] * cur[j]
                               - vnew * M[j][idx] * M[j][idx])
                        if not det.is_zero():
                            ok = False
                            break
                    if ok:
                        cur.append(c)
                        bt(idx + 1, cur)
                        cur.pop()
            bt(0, [])

    search([], 0)
    return best[0]


# ----------------------------------------------------------------------------
# Realization of an abstract Gram matrix as vectors in K^2
# ----------------------------------------------------------------------------
def realize_vectors(M):
    """Pick an invertible 2x2 principal block (p,q) of minimal trace; coordinates of
    every v_i w.r.t. the dual-of-pivot frame, so that B(v_i,v_j) = M[i][j] with
    B given by Gp. The v_i need NOT be O_K-integral in this frame: the honest
    lattice is L = sum_i O_K v_i (see realized_module), NOT O_K v_p + O_K v_q."""
    s = len(M)
    cands = []
    for p in range(s):
        for q in range(p + 1, s):
            det = M[p][p] * M[q][q] - M[p][q] * M[q][p]
            if not det.is_zero():
                cands.append((M[p][p].trace() + M[q][q].trace(), p, q))
    cands.sort(key=lambda t: t[0])
    _, p, q = cands[0]
    Gp = [[M[p][p], M[p][q]], [M[q][p], M[q][q]]]
    di = (Gp[0][0] * Gp[1][1] - Gp[0][1] * Gp[1][0]).inv()
    vs = []
    for i in range(s):
        r0, r1 = M[i][p], M[i][q]
        vs.append(((Gp[1][1] * r0 - Gp[0][1] * r1) * di,
                   (Gp[0][0] * r1 - Gp[1][0] * r0) * di))
    return vs, Gp


def Bform(G, P, R):
    (x1, y1), (x2, y2) = P, R
    return G[0][0] * x1 * x2 + G[0][1] * (x1 * y2 + x2 * y1) + G[1][1] * y1 * y2


# ----------------------------------------------------------------------------
# Exact integer lattice tools
# ----------------------------------------------------------------------------
def _xgcd(a, b):
    x0, y0, x1, y1 = 1, 0, 0, 1
    while b:
        qq = a // b
        a, b = b, a - qq * b
        x0, x1 = x1, x0 - qq * x1
        y0, y1 = y1, y0 - qq * y1
    return a, x0, y0


def hnf(rows):
    """Row Hermite normal form (echelon, positive pivots, reduced above) of the
    Z-row-span of integer `rows`. Returns list of nonzero rows (a Z-basis)."""
    A = [list(r) for r in rows if any(r)]
    if not A:
        return []
    ncol = len(A[0])
    H = []
    r = 0
    for c in range(ncol):
        # gather rows (from r on) with nonzero in column c; gcd-combine into row r
        idx = [i for i in range(r, len(A)) if A[i][c] != 0]
        if not idx:
            continue
        i0 = idx[0]
        A[r], A[i0] = A[i0], A[r]
        for i in range(r + 1, len(A)):
            if A[i][c] == 0:
                continue
            g, x, y = _xgcd(A[r][c], A[i][c])
            a_, b_ = A[r][c] // g, A[i][c] // g
            new_r = [x * A[r][j] + y * A[i][j] for j in range(ncol)]
            new_i = [-b_ * A[r][j] + a_ * A[i][j] for j in range(ncol)]
            A[r], A[i] = new_r, new_i
        if A[r][c] < 0:
            A[r] = [-t for t in A[r]]
        for i in range(r):                      # reduce above
            qq = A[i][c] // A[r][c]
            if qq:
                A[i] = [A[i][j] - qq * A[r][j] for j in range(ncol)]
        r += 1
        if r == len(A):
            break
    return [row for row in A[:r] if any(row)]


def in_lattice(H, v):
    """Is integer vector v in the Z-span of HNF basis H?"""
    v = list(v)
    for row in H:
        c = next(j for j, t in enumerate(row) if t != 0)
        if v[c] % row[c]:
            return False
        qq = v[c] // row[c]
        if qq:
            v = [v[j] - qq * row[j] for j in range(len(v))]
    return not any(v)


def det_hnf(H):
    if not H or len(H) < len(H[0]):
        return 0
    d = 1
    for row in H:
        d *= next(t for t in row if t != 0)
    return abs(d)


# ----------------------------------------------------------------------------
# The realized module L = sum O_K v_i and its trace lattice
# ----------------------------------------------------------------------------
def _to_q4(v):
    x, y = v
    return [x.a, x.b, y.a, y.b]


def realized_module(D, M):
    """Return (Gp, vs, basis) with basis a Z-basis (as K^2 vectors) of the
    O_K-module L = sum_i O_K v_i  (generators v_i and omega v_i over Z; O_K = Z + Z omega).
    (Before the R1 session this used sqrt(D) v_i, i.e. Z[sqrt D]-span: an index-4 non-O_K-lattice
    when D = 1 mod 4 -- see docs/03_corrections.md C17.)"""
    vs, Gp = realize_vectors(M)
    sd = omega(D)
    gens = []
    for v in vs:
        gens.append(_to_q4(v))
        gens.append(_to_q4((sd * v[0], sd * v[1])))
    den = 1
    for g in gens:
        for c in g:
            den = den * F(c).denominator // gcd(den, F(c).denominator)
    H = hnf([[int(F(c) * den) for c in g] for g in gens])
    assert len(H) == 4, "L not of Z-rank 4"
    basis = [(QF(D, F(h[0], den), F(h[1], den)), QF(D, F(h[2], den), F(h[3], den))) for h in H]
    return Gp, vs, basis


def trace_gram(Gp, basis):
    T = [[Bform(Gp, basis[a], basis[b]).trace() for b in range(4)] for a in range(4)]
    for row in T:
        for t in row:
            assert t.denominator == 1, "trace form not integral on L"
    return [[int(t) for t in row] for row in T]


def vec_value(D, Gp, basis, vv):
    x = QF(D, 0, 0)
    y = QF(D, 0, 0)
    for a in range(4):
        x = x + basis[a][0] * vv[a]
        y = y + basis[a][1] * vv[a]
    return Bform(Gp, (x, y), (x, y))


# ----------------------------------------------------------------------------
# Fincke-Pohst enumeration (exact acceptance)
# ----------------------------------------------------------------------------
def _cholesky_upper(T):
    n = len(T)
    R = [[0.0] * n for _ in range(n)]
    for i in range(n):
        s = T[i][i] - sum(R[k][i] ** 2 for k in range(i))
        R[i][i] = sqrt(max(s, 1e-300))
        for j in range(i + 1, n):
            R[i][j] = (T[i][j] - sum(R[k][i] * R[k][j] for k in range(i))) / R[i][i]
    return R


def qnorm(T, v):
    n = len(T)
    return sum(T[i][j] * v[i] * v[j] for i in range(n) for j in range(n))


def enum_short(T, cap, pad=1e-6):
    """All nonzero integer v with v^T T v <= cap, one per +-pair (first nonzero > 0).
    Returns list of (norm, v) sorted by norm. Norms exact."""
    n = len(T)
    R = _cholesky_upper(T)
    d = [R[i][i] ** 2 for i in range(n)]
    q = [[(R[i][j] / R[i][i] if j > i else 0.0) for j in range(n)] for i in range(n)]
    v = [0] * n
    out = []

    def rec(i, rem):
        if i < 0:
            if any(v):
                fnz = next(x for x in v if x)
                if fnz > 0:
                    nv = qnorm(T, v)
                    if 0 < nv <= cap:
                        out.append((nv, tuple(v)))
            return
        c = sum(q[i][l] * v[l] for l in range(i + 1, n))
        bound = sqrt(max(rem, 0.0) / d[i]) + pad
        for vi in range(ceil(-bound - c), floor(bound - c) + 1):
            v[i] = vi
            rec(i - 1, rem - d[i] * (vi + c) ** 2)
        v[i] = 0

    rec(n - 1, cap + pad * cap + pad)
    out.sort()
    return out


# ----------------------------------------------------------------------------
# Represented indecomposables, orbits, square classes, linear irreducibility
# ----------------------------------------------------------------------------
def indecomposable_test(D):
    """Return predicate is_indec(alpha) for K = Q(sqrt D) (exact, via orbit reps)."""
    eta = tp_unit(D)
    keys = set()
    for b in indecomposables(D):
        nb = normalize_unit(b, eta)
        keys.add((nb.a, nb.b))
    normcap = D if D % 4 != 1 else F(D, 4)      # Dress-Scharlau: N(alpha) <= Delta/4

    def is_indec(al):
        if not al.is_tot_pos() or al.norm() > normcap:
            return False
        na = normalize_unit(al, eta)
        return (na.a, na.b) in keys
    return is_indec, eta


def analyze_lattice(D, Gp, basis, cap):
    """
    Enumerate the trace lattice up to trace-norm `cap`. Returns dict with
      orbits   : set of tp-unit-orbit keys of represented indecomposables
      classes  : list of square-class reps represented (first hit, min trace)
      lin_irred_classes : classes having a Z-linearly-irreducible representative
                          (v not in Z-span of trace-shorter vectors)
      nstar    : trace-norm at which shorter vectors first Z-span L (det 1), or None
    Exhaustive for vectors of trace <= cap. A class whose minimal representing
    trace exceeds cap is missed: report cap with results.
    """
    T = trace_gram(Gp, basis)
    is_indec, eta = indecomposable_test(D)
    sv = enum_short(T, cap)
    by = {}
    for nv, vv in sv:
        by.setdefault(nv, []).append(vv)
    span_rows = []
    H = []
    orbits, classes, lin = set(), [], []
    nstar = None
    for nv in sorted(by):
        for vv in by[nv]:
            val = vec_value(D, Gp, basis, vv)
            if not is_indec(val):
                continue
            na = normalize_unit(val, eta)
            orbits.add((na.a, na.b))
            if not any(same_square_class(val, c) for c, _ in classes):
                classes.append((val, nv))
            if not in_lattice(H, vv):          # H = span of strictly shorter vectors
                if not any(same_square_class(val, c) for c, _ in lin):
                    lin.append((val, nv))
        span_rows += [list(x) for x in by[nv]]
        H = hnf(span_rows + H)
        span_rows = []
        if nstar is None and len(H) == 4 and det_hnf(H) == 1:
            nstar = nv
    return dict(T=T, orbits=orbits, classes=classes, lin_irred_classes=lin,
                nstar=nstar, cap=cap, n_short=len(sv))


# ----------------------------------------------------------------------------
# Exact, cap-free s_sq(M) via the per-face lemma (ledger L1.4/L1.5)
# ----------------------------------------------------------------------------
@lru_cache(maxsize=None)
def edge_deltas(D):
    """Codifferent functionals of the sail edges, one per edge modulo (O^x)^2.
    Edge E_i = {alpha_i + t*alpha_{i+1} : 0<=t<=u_{i+1}} (i even, indec.py conventions);
    delta_E is the unique delta with Tr(delta*alpha_i)=1, Tr(delta*alpha_{i+1})=0.
    Asserts delta_E is totally positive and lies in the codifferent. Cached (do not mutate)."""
    from indec import cf_list, convergents
    from rqf import fundamental_unit
    eps, _ = fundamental_unit(D)
    u = eps * eps                                  # generator of (O^x)^2 (tot. pos.)
    cf, dl, s = cf_list(D, 2 * 2 * 60 + 14)
    nterms = 4 * s + 12
    cvs = convergents(cf, nterms + 1)
    alpha = [QF(D, 1, 0)] + [QF(D, p, 0) + QF(D, q, 0) * dl for (p, q) in cvs]
    w = omega(D)
    reps, keys = [], set()
    for i in range(0, 4 * s + 8, 2):
        a0, a1 = alpha[i], alpha[i + 1]
        # unknown delta = x + y sqrt D ; Tr(delta*a) = 2(x*a.a + y*a.b*D)
        A = [[2 * a0.a, 2 * a0.b * D], [2 * a1.a, 2 * a1.b * D]]
        det = A[0][0] * A[1][1] - A[0][1] * A[1][0]
        x = (1 * A[1][1] - 0 * A[0][1]) / det
        y = (A[0][0] * 0 - A[1][0] * 1) / det
        dlt = QF(D, x, y)
        assert dlt.is_tot_pos(), (D, i, dlt)
        assert (dlt.trace()).denominator == 1 and ((dlt * w).trace()).denominator == 1, (D, i, dlt)
        nd = normalize_unit(dlt, u)
        if (nd.a, nd.b) not in keys:
            keys.add((nd.a, nd.b))
            reps.append(nd)
    return reps


def exact_sq_classes(D, Gp, basis):
    """All square classes of indecomposables represented by the lattice with Z-basis `basis`
    (as K^2 vectors) and form Gp. Exact and complete: by L1.4, an indecomposable on edge E is
    Q(v) with T_{delta_E}(v) = 1, and conversely T_delta(v)=1 forces Q(v) indecomposable.
    Returns (classes, per_edge) with per_edge[k] = #values on the k-th edge rep."""
    classes, per_edge = [], []
    for dlt in edge_deltas(D):
        Td = [[(dlt * Bform(Gp, basis[a], basis[b])).trace() for b in range(4)] for a in range(4)]
        assert all(t.denominator == 1 for r in Td for t in r)
        Td = [[int(t) for t in r] for r in Td]
        vals = []
        for nv, vv in enum_short_reduced(Td, 1):
            val = vec_value(D, Gp, basis, vv)
            if val not in vals:
                vals.append(val)
            if not any(same_square_class(val, c) for c in classes):
                classes.append(val)
        per_edge.append(len(vals))
    return classes, per_edge


def lll_gram(G, delta=F(3, 4)):
    """Exact LLL on a positive-definite integer Gram matrix G (all arithmetic in Fraction).
    Returns (Gred, U) with Gred = U G U^T, U unimodular (rows = new basis in old coords)."""
    n = len(G)
    U = [[int(i == j) for j in range(n)] for i in range(n)]

    def gram(U):
        return [[F(sum(U[i][a] * G[a][b] * U[j][b] for a in range(n) for b in range(n)))
                 for j in range(n)] for i in range(n)]

    def gso(Gc):
        mu = [[F(0)] * n for _ in range(n)]
        Bn = [F(0)] * n
        for i in range(n):
            for j in range(i):
                mu[i][j] = (Gc[i][j] - sum(mu[j][k] * mu[i][k] * Bn[k] for k in range(j))) / Bn[j]
            Bn[i] = Gc[i][i] - sum(mu[i][k] ** 2 * Bn[k] for k in range(i))
        return mu, Bn

    k = 1
    while k < n:
        for j in range(k - 1, -1, -1):
            mu, Bn = gso(gram(U))
            q = round(mu[k][j])
            if q:
                U[k] = [U[k][t] - q * U[j][t] for t in range(n)]
        mu, Bn = gso(gram(U))
        if Bn[k] >= (delta - mu[k][k - 1] ** 2) * Bn[k - 1]:
            k += 1
        else:
            U[k], U[k - 1] = U[k - 1], U[k]
            k = max(k - 1, 1)
    Gred = [[int(x) for x in row] for row in gram(U)]
    return Gred, U


def enum_short_reduced(T, cap):
    """enum_short after exact LLL; returns vectors in the ORIGINAL coordinates."""
    Tr, U = lll_gram(T)
    out = []
    n = len(T)
    for nv, w in enum_short(Tr, cap):
        v = tuple(sum(w[i] * U[i][j] for i in range(n)) for j in range(n))
        assert qnorm(T, v) == nv
        out.append((nv, v))
    return out
