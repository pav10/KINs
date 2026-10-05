"""
Roadmap R1: provably complete computation of
    S(K,2) = max s_sq(M)   over ALL classically integral (B(M,M) in O), totally positive
                           definite binary O_K-lattices M (free or not),  K = Q(sqrt D).
Completeness argument: docs/09_R1_complete.md (ledger L7.1).

Pipeline
  1. unit_square_reps(D): indecomposables modulo (O^x)^2 (rescaling v by a unit eps changes
     Q(v) by eps^2 and leaves the line O v unchanged).
  2. frames(D): (i, j, b): alpha = R[i], beta = R[j] in distinct square classes, b in O with
     alpha*beta - b^2 >> 0. Reduced modulo b -> -b (same lattice), swap (isometric lattice) and
     Galois conjugation (conjugate lattice, same s_sq).
  3. Frame lattice L = O e1 + O e2, Gram G = [[alpha, b], [b, beta]].  L^# = G^{-1} O^2, so
     A := L^#/L = O^2 / G O^2 (|A| = N(alpha*beta - b^2)) with the K/O-valued pairing
     b_A(y, y') = y^T G^{-1} y' mod O.  Integral overlattices M of L  <->  b_A-isotropic
     O-submodules of A.  Only the MAXIMAL ones are needed (s_sq is monotone under inclusion).
  4. s_sq of each maximal overlattice exactly and cap-free (lattice.exact_sq_classes, L6.1).

All decisions are exact (int / Fraction).  Floats only bound search boxes, with exact re-checks
(gram.valid_off_diagonal for b; principal-ideal search in is_free).
"""
from fractions import Fraction as F
from itertools import product
from math import sqrt, gcd, floor, ceil, isqrt
from rqf import QF, fundamental_unit, same_square_class
from indec import indecomposables
from gram import valid_off_diagonal, omega   # valid_off_diagonal: cross-check only
from lattice import hnf, exact_sq_classes


# ----------------------------------------------------------------------------
# O_K in integer coordinates: (c0, c1) <-> c0 + c1*omega,  omega^2 = t*omega + n
# ----------------------------------------------------------------------------
class OArith:
    def __init__(self, D):
        self.D = D
        if D % 4 == 1:
            self.t, self.n = 1, (D - 1) // 4
        else:
            self.t, self.n = 0, D

    def mul(self, x, y):
        a, b = x
        c, d = y
        bd = b * d
        return (a * c + bd * self.n, a * d + b * c + bd * self.t)

    def add(self, x, y): return (x[0] + y[0], x[1] + y[1])
    def sub(self, x, y): return (x[0] - y[0], x[1] - y[1])
    def neg(self, x):    return (-x[0], -x[1])
    def conj(self, x):   return (x[0] + x[1] * self.t, -x[1])      # omega' = t - omega

    def norm(self, x):
        p = self.mul(x, self.conj(x))
        assert p[1] == 0
        return p[0]

    def from_qf(self, q):
        if self.D % 4 == 1:                     # a + b sqrt D = (a - b) + 2b omega
            c0, c1 = q.a - q.b, 2 * q.b
        else:
            c0, c1 = q.a, q.b
        assert c0.denominator == 1 and c1.denominator == 1, ("not in O_K", q)
        return (int(c0), int(c1))

    def to_qf(self, x):
        if self.D % 4 == 1:
            return QF(self.D, x[0] + F(x[1], 2), F(x[1], 2))
        return QF(self.D, x[0], x[1])


# ----------------------------------------------------------------------------
# Step 1: indecomposables modulo (O^x)^2, with exponent tracking
# ----------------------------------------------------------------------------
def normalize_track(x, u):
    """(rep, k) with rep = x * u^k in the fundamental domain ratio in [1, u1^2) for the
    totally positive unit u > 1 (same domain as indec.normalize_unit)."""
    uc = u.conj()
    k = 0
    for _ in range(100000):
        if x.b < 0:
            x, k = x * u, k + 1
        elif (x * uc).b >= 0:
            x, k = x * uc, k - 1
        else:
            return x, k
    raise RuntimeError("normalization did not converge")


def unit_power(eps, Neps, k):
    """eps^k for k in Z (eps^-1 = N(eps) * eps')."""
    base = eps if k >= 0 else eps.conj() * Neps
    r = QF(eps.D, 1, 0)
    for _ in range(abs(k)):
        r = r * base
    return r


def emb_stable(x):
    """Float embeddings (sigma1, sigma2) of x != 0, the smaller one recovered as N(x)/larger
    (avoids the cancellation in a - b sqrt D for lopsided x).  Floats: box bounds only."""
    e1, e2 = x.emb()
    n = float(x.norm())
    if abs(e1) >= abs(e2):
        return e1, n / e1
    return n / e2, e2


def balance_track(x, u):
    """(r, j): r = x * u^j of minimal trace (u totally positive unit > 1)."""
    uc = u.conj()
    j = 0
    while True:
        if (x * u).trace() < x.trace():
            x, j = x * u, j + 1
        elif (x * uc).trace() < x.trace():
            x, j = x * uc, j - 1
        else:
            return x, j


class FieldData:
    """Everything about K = Q(sqrt D) the R1 search needs."""

    def __init__(self, D):
        self.D = D
        self.ar = OArith(D)
        self.eps, self.Neps = fundamental_unit(D)
        self.u = self.eps * self.eps                       # generator of (O^x)^2
        eta = self.eps if self.Neps > 0 else self.u        # generator of O^{x,+}
        pool = []
        for a in indecomposables(D):                       # reps modulo O^{x,+}
            pool.append(a)
            if self.Neps > 0:                              # [O^{x,+} : (O^x)^2] = 2
                pool.append(a * eta)
        reps, seen = [], set()
        for a in pool:
            c, _ = normalize_track(a, self.u)              # canonical form of the orbit
            if (c.a, c.b) not in seen:
                seen.add((c.a, c.b))
                r, j = balance_track(c, self.u)            # r = c * u^j, minimal trace
                reps.append((r, c, j))
        reps.sort(key=lambda t: (t[0].trace(), t[0].b))
        self.reps = [r for r, _, _ in reps]                # indecomposables mod (O^x)^2 (balanced)
        self.index = {(c.a, c.b): (i, j) for i, (_, c, j) in enumerate(reps)}
        # square class id of each rep (mod (K^x)^2)
        self.class_of = []
        class_reps = []
        for r in self.reps:
            for c, rc in enumerate(class_reps):
                if same_square_class(r, rc):
                    self.class_of.append(c)
                    break
            else:
                class_reps.append(r)
                self.class_of.append(len(class_reps) - 1)
        self.class_reps = class_reps                       # kappa_sq(K) many

    def normalize(self, x):
        """(i, k): x * u^k = reps[i]  (x indecomposable)."""
        c, k = normalize_track(x, self.u)
        i, j = self.index[(c.a, c.b)]
        return i, k + j


# ----------------------------------------------------------------------------
# Step 2: frames modulo symmetry
# ----------------------------------------------------------------------------
def _canon(fd, i, j, b):
    """Canonical form of the frame (R[i], R[j], b) under swap and b -> -b."""
    if i > j:
        i, j = j, i
    bc = fd.ar.from_qf(b)
    nb = fd.ar.neg(bc)
    return (i, j) + max(bc, nb)


def _galois(fd, i, j, b):
    """Canonical key of the Galois-conjugate frame, renormalized to the reps."""
    ia, ka = fd.normalize(fd.reps[i].conj())
    ib, kb = fd.normalize(fd.reps[j].conj())
    # Q(eps^ka v') = R[ia], Q(eps^kb w') = R[ib], B = eps^(ka+kb) b'
    b2 = b.conj() * unit_power(fd.eps, fd.Neps, ka + kb)
    return _canon(fd, ia, ib, b2)


def admissible_b(fd, al, be):
    """All b in O with al*be - b^2 totally positive (exact).  The box |sigma_k(b)| < B_k,
    B_k = sqrt(sigma_k(al*be)), is first balanced by b = eps^-k * c (floats only choose k and
    bound the search box, padded by 2; every candidate is re-checked exactly)."""
    from math import log
    ar = fd.ar
    p = ar.from_qf(al * be)
    e1 = max(abs(x) for x in fd.eps.emb())
    P1, P2 = emb_stable(al * be)
    B1, B2 = sqrt(P1), sqrt(P2)
    k = round(-log(B1 / B2) / (2 * log(e1)))
    ek = unit_power(fd.eps, fd.Neps, k)                # c = b * eps^k
    ekc = ar.from_qf(unit_power(fd.eps, fd.Neps, -k))  # b = c * eps^-k
    s1, s2 = (abs(x) for x in emb_stable(ek))
    C1, C2 = B1 * s1, B2 * s2
    w1, w2 = omega(fd.D).emb()
    dw = abs(w1 - w2)
    ymax = (C1 + C2) / dw
    out = []
    for y in range(-int(ymax) - 2, int(ymax) + 3):
        lo = max(-C1 - y * w1, -C2 - y * w2)
        hi = min(C1 - y * w1, C2 - y * w2)
        if lo > hi + 4:
            continue
        for x in range(floor(lo) - 2, ceil(hi) + 3):
            b = ar.mul((x, y), ekc)
            d = ar.sub(p, ar.mul(b, b))
            if 2 * d[0] + d[1] * ar.t > 0 and ar.norm(d) > 0:
                out.append(ar.to_qf(b))
    return out


def frames(fd, galois=True):
    """Yield (key, alpha, beta, b) for one frame per symmetry orbit."""
    R = fd.reps
    for i in range(len(R)):
        for j in range(i + 1, len(R)):
            if fd.class_of[i] == fd.class_of[j]:
                continue
            al, be = R[i], R[j]
            seen = set()
            for b in admissible_b(fd, al, be):
                key = _canon(fd, i, j, b)
                if key in seen:
                    continue
                seen.add(key)
                if galois and _galois(fd, i, j, b) < key:
                    continue
                yield key, al, be, b


# ----------------------------------------------------------------------------
# Step 3: discriminant module A = O^2 / G O^2 and its maximal isotropic submodules
# ----------------------------------------------------------------------------
class DiscModule:
    """A = L^#/L for L = O^2 with Gram G = [[al, b], [b, be]], realized as O^2/G O^2 via
    y <-> G^{-1} y.  Elements: canonical integer 4-tuples (y1c0, y1c1, y2c0, y2c1)."""

    def __init__(self, ar, al, b, be):
        self.ar = ar
        m = ar.mul
        w = (0, 1)
        self.al, self.b, self.be = al, b, be
        self.det = ar.sub(m(al, be), m(b, b))
        self.N = ar.norm(self.det)
        assert self.N > 0
        cols = [(al, b), (m(al, w), m(b, w)), (b, be), (m(b, w), m(be, w))]
        H = hnf([[c[0][0], c[0][1], c[1][0], c[1][1]] for c in cols])
        assert len(H) == 4 and all(H[i][i] > 0 for i in range(4))
        dprod = H[0][0] * H[1][1] * H[2][2] * H[3][3]
        assert dprod == self.N, (dprod, self.N)
        self.H = H
        self.zero = (0, 0, 0, 0)
        # pairing: b_A(y, y') = y'^T adj(G) y / det; test via det' * adj(G) y, coords mod N
        self.detc = ar.conj(self.det)

    def reduce(self, y):
        y = list(y)
        H = self.H
        for i in range(4):
            q = y[i] // H[i][i]
            if q:
                h = H[i]
                for j in range(i, 4):
                    y[j] -= q * h[j]
        return tuple(y)

    def elements(self):
        H = self.H
        return product(range(H[0][0]), range(H[1][1]), range(H[2][2]), range(H[3][3]))

    def add(self, x, y):
        return self.reduce((x[0] + y[0], x[1] + y[1], x[2] + y[2], x[3] + y[3]))

    def omega(self, y):
        m = self.ar.mul
        a, c = m((y[0], y[1]), (0, 1)), m((y[2], y[3]), (0, 1))
        return self.reduce((a[0], a[1], c[0], c[1]))

    def dual_vec(self, y):
        """w(y) = det' * adj(G) y (two O-elements, coords mod N)."""
        m, ar = self.ar.mul, self.ar
        y1, y2 = (y[0], y[1]), (y[2], y[3])
        z1 = ar.sub(m(self.be, y1), m(self.b, y2))
        z2 = ar.sub(m(self.al, y2), m(self.b, y1))
        z1, z2 = m(self.detc, z1), m(self.detc, z2)
        N = self.N
        return (z1[0] % N, z1[1] % N, z2[0] % N, z2[1] % N)

    def pair_zero(self, x, wy):
        """b_A(x, y) == 0 in K/O, given wy = dual_vec(y)."""
        m = self.ar.mul
        p = self.ar.add(m((x[0], x[1]), (wy[0], wy[1])), m((x[2], x[3]), (wy[2], wy[3])))
        return p[0] % self.N == 0 and p[1] % self.N == 0

    def cyclic(self, x):
        """O x = Z x + Z omega x inside A."""
        g = [x, self.omega(x)]
        out = {self.zero}
        frontier = [self.zero]
        while frontier:
            nxt = []
            for e in frontier:
                for h in g:
                    f = self.add(e, h)
                    if f not in out:
                        out.add(f)
                        nxt.append(f)
            frontier = nxt
        return out

    def primary_parts(self):
        """{p: elements of the p-primary part A_p}.  A = (+)_p A_p, and the A_p are mutually
        b_A-orthogonal (p^a b, q^c b in O with gcd(p^a, q^c) = 1 => b in O)."""
        N, parts = self.N, {}
        n, p = N, 2
        primes = []
        while p * p <= n:
            if n % p == 0:
                primes.append(p)
                while n % p == 0:
                    n //= p
            p += 1
        if n > 1:
            primes.append(n)
        for p in primes:
            pv = 1
            while N % (pv * p) == 0:
                pv *= p
            m = N // pv
            gens = []
            for e in ((1, 0, 0, 0), (0, 1, 0, 0), (0, 0, 1, 0), (0, 0, 0, 1)):
                g = self.reduce(tuple(m * c for c in e))
                if g != self.zero:
                    gens.append(g)
            Ap = {self.zero}
            frontier = [self.zero]
            while frontier:
                nxt = []
                for x in frontier:
                    for g in gens:
                        y = self.add(x, g)
                        if y not in Ap:
                            Ap.add(y)
                            nxt.append(y)
                frontier = nxt
            assert len(Ap) == pv, (p, len(Ap), pv)
            parts[p] = Ap
        return parts

    def _isotropic_part(self, elems):
        """ALL isotropic O-submodules inside the O-submodule `elems` (a set), as a list of
        (S, gens, is_maximal), S a frozenset.  DFS from 0; complete: if S < S' are isotropic,
        any x in S' minus S is isotropic and orthogonal to S with S + Ox <= S', and every such
        x is tried at node S -- so every isotropic submodule (in particular every maximal one)
        is reached."""
        iso = []
        for y in elems:
            wy = self.dual_vec(y)
            if y != self.zero and self.pair_zero(y, wy):
                iso.append(y)
        start = frozenset([self.zero])
        stack = [(start, [])]
        seen = {start}
        result = []
        while stack:
            S, gens = stack.pop()
            gw = [self.dual_vec(g) for g in gens]
            cands = [x for x in iso if x not in S and all(self.pair_zero(x, w) for w in gw)]
            result.append((S, gens, not cands))
            for x in cands:
                C = self.cyclic(x)
                S2 = frozenset(self.add(s, c) for s in S for c in C)
                if S2 not in seen:
                    seen.add(S2)
                    stack.append((S2, gens + [x]))
        return result

    def isotropic_submodules(self):
        """All b_A-isotropic O-submodules of A as (parts, gens, is_maximal); parts = tuple of
        the p-primary components (frozensets).  A submodule is the sum of its p-parts and the
        A_p are mutually orthogonal, so the isotropic submodules are exactly the products."""
        out = [((), [], True)]
        for p, Ap in sorted(self.primary_parts().items()):
            loc = self._isotropic_part(Ap)
            out = [(P + (S,), g + h, m1 and m2) for P, g, m1 in out for S, h, m2 in loc]
        return out

    def maximal_isotropic(self):
        """All maximal b_A-isotropic O-submodules of A, as lists of O-generators."""
        return [g for _, g, m in self.isotropic_submodules() if m]


def _contains(P, Q):
    """parts tuple P contains parts tuple Q (componentwise)."""
    return all(q <= p for p, q in zip(P, Q))


# ----------------------------------------------------------------------------
# Step 4: realize overlattices, evaluate exactly
# ----------------------------------------------------------------------------
def _q4(v):
    return [v[0].a, v[0].b, v[1].a, v[1].b]


def overlattice_basis(fd, G, gens):
    """Z-basis (K^2 vectors) of M = O^2 + sum O G^{-1} y_g, y_g in gens."""
    D = fd.D
    one, z, w = QF(D, 1), QF(D, 0), omega(D)
    det = G[0][0] * G[1][1] - G[0][1] * G[1][0]
    di = det.inv()
    vecs = [(one, z), (w, z), (z, one), (z, w)]
    for y in gens:
        y1, y2 = fd.ar.to_qf((y[0], y[1])), fd.ar.to_qf((y[2], y[3]))
        x = ((G[1][1] * y1 - G[0][1] * y2) * di, (G[0][0] * y2 - G[1][0] * y1) * di)
        vecs += [x, (w * x[0], w * x[1])]
    rows = [_q4(v) for v in vecs]
    den = 1
    for r in rows:
        for c in r:
            d = F(c).denominator
            den = den * d // gcd(den, d)
    Hm = hnf([[int(F(c) * den) for c in r] for r in rows])
    assert len(Hm) == 4
    return [(QF(D, F(h[0], den), F(h[1], den)), QF(D, F(h[2], den), F(h[3], den))) for h in Hm]


def steinitz_ideal(fd, G, basis):
    """Z-basis of the ideal generated by det(x, y), x, y in M (Steinitz class of M)."""
    D = fd.D
    dets = []
    for i in range(4):
        for j in range(i + 1, 4):
            dets.append(basis[i][0] * basis[j][1] - basis[i][1] * basis[j][0])
    rows = [[d.a, d.b] for d in dets]
    den = 1
    for r in rows:
        for c in r:
            dd = F(c).denominator
            den = den * dd // gcd(den, dd)
    Hm = hnf([[int(F(c) * den) for c in r] for r in rows])
    assert len(Hm) == 2
    return [QF(D, F(h[0], den), F(h[1], den)) for h in Hm]


def _cf_cycle(theta):
    """Complete quotients (as exact K-elements (a, b) = a + b sqrt D) of the periodic part of the
    continued fraction of the real quadratic irrational theta = a + b sqrt D (b != 0)."""
    D = theta.D
    A, B = theta.a, theta.b
    den = A.denominator * B.denominator // gcd(A.denominator, B.denominator)
    P, m, Q = int(A * den), int(B * den), den          # theta = (P + m sqrt D) / Q
    if m < 0:
        P, m, Q = -P, -m, -Q
    d = m * m * D                                      # theta = (P + sqrt d) / Q
    if (d - P * P) % Q:
        P, d, m, Q = P * abs(Q), d * Q * Q, m * abs(Q), Q * abs(Q)
    r = isqrt(d)
    seen, states = {}, []
    while (P, Q) not in seen:
        seen[(P, Q)] = len(states)
        states.append((P, Q))
        a = (P + r) // Q if Q > 0 else -((P + r) // (-Q) + 1)      # floor((P + sqrt d)/Q)
        P = a * Q - P
        Q = (d - P * P) // Q
    return {(F(Pi, Qi), F(m, Qi)) for Pi, Qi in states[seen[(P, Q)]:]}


_omega_cycle = {}


def is_principal(fd, I):
    """Is the fractional O_K-ideal with Z-basis I = [g1, g2] principal?  I = g1 (Z + Z theta),
    theta = g2/g1; principal  <=>  Z + Z theta = lambda O  <=>  theta ~ omega under GL_2(Z)
    <=>  (Serret) their continued fractions share a complete quotient in the period.  Exact."""
    D = fd.D
    if D not in _omega_cycle:
        _omega_cycle[D] = _cf_cycle(omega(D))
    return not _omega_cycle[D].isdisjoint(_cf_cycle(I[1] / I[0]))


def is_principal_box(fd, I):
    """Reference (slow for large eps): search x in I with |N(x)| = N(I) in the box
    |sigma_i(x)| <= sqrt(N(I) eps1) (a generator can be moved there by a power of eps)."""
    D = fd.D
    g1, g2 = I
    w = omega(D)
    NI = abs(g1.a * g2.b - g1.b * g2.a) / w.b           # covolume ratio vs O
    e1 = max(abs(x) for x in fd.eps.emb())
    R = sqrt(float(NI) * e1) * (1 + 1e-9) + 1e-9
    (a11, a21), (a12, a22) = emb_stable(g1), emb_stable(g2)   # sigma_k(m g1 + n g2)
    det = a11 * a22 - a12 * a21
    mb = R * (abs(a22) + abs(a12)) / abs(det) + 2
    nb = R * (abs(a21) + abs(a11)) / abs(det) + 2
    for mm in range(-int(mb) - 1, int(mb) + 2):
        for nn in range(-int(nb) - 1, int(nb) + 2):
            x = g1 * mm + g2 * nn
            if not x.is_zero() and abs(x.norm()) == NI:
                return True
    return False


def evaluate(fd, G, gens):
    basis = overlattice_basis(fd, G, gens)
    cl, pe = exact_sq_classes(fd.D, G, basis)
    return len(cl), cl, pe, basis


# ----------------------------------------------------------------------------
# Driver
# ----------------------------------------------------------------------------
def complete_S(D, galois=True, free_too=True, verbose=False, record_all=False):
    """Exact S(K,2) for K = Q(sqrt D) over ALL integral binary lattices (max over the maximal
    integral overlattices of all frames) and, if free_too, exact S_free(K,2) over FREE lattices
    (= binary forms): a free M is contained in a lattice maximal among the FREE integral
    overlattices of its frame, so those are evaluated (a maximal overlattice need not be free).
    Returns a dict."""
    import time
    t0 = time.time()
    fd = FieldData(D)
    ar = fd.ar
    best, best_free = 1, 1                     # <1> + <1> represents the class of 1
    best_rec, best_free_rec = None, None
    nfr = nlat = nfree = 0
    hist = {}
    allrec = []
    for key, al, be, b in frames(fd, galois=galois):
        nfr += 1
        A = DiscModule(ar, ar.from_qf(al), ar.from_qf(b), ar.from_qf(be))
        G = [[al, b], [b, be]]
        subs = A.isotropic_submodules()
        cache = {}

        def ev(P, gens, maximal):
            if P not in cache:
                n, cl, pe, basis = evaluate(fd, G, gens)
                rec = dict(frame=[str(al), str(be), str(b)], gens=[list(g) for g in gens],
                           s_sq=n, per_edge=pe, classes=[str(c) for c in cl], detN=A.N,
                           maximal=maximal)
                cache[P] = (n, rec)
            return cache[P]

        for P, gens, m in subs:
            if not m:
                continue
            nlat += 1
            n, rec = ev(P, gens, True)
            hist[n] = hist.get(n, 0) + 1
            if n > best:
                best, best_rec = n, rec
            if record_all:
                allrec.append(rec)
        if free_too:
            free = [(P, gens, m) for P, gens, m in subs
                    if is_principal(fd, steinitz_ideal(fd, G, overlattice_basis(fd, G, gens)))]
            for P, gens, m in free:
                if any(Q != P and _contains(Q, P) for Q, _, _ in free):
                    continue                   # not maximal among free overlattices
                nfree += 1
                n, rec = ev(P, gens, m)
                if n > best_free:
                    best_free, best_free_rec = n, dict(rec, free=True)
        if verbose and nfr % 200 == 0:
            print(f"  D={D}: {nfr} frames, {nlat} lattices, best {best}, "
                  f"{time.time() - t0:.0f}s", flush=True)
    out = dict(D=D, n_reps=len(fd.reps), kappa=len(fd.class_reps), frames=nfr, lattices=nlat,
               S=best, best=best_rec, hist={str(k): v for k, v in sorted(hist.items())},
               sec=round(time.time() - t0, 1))
    if free_too:
        out["S_free"] = best_free
        out["best_free"] = best_free_rec
        out["free_lattices"] = nfree
    if record_all:
        out["all"] = allrec
    return out
