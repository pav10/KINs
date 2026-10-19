#!/usr/bin/env python3
"""R2 probe (HEURISTIC predictor, verified against exact data): distribute the indecomposables of
K = Q(sqrt D) over the ideal classes of E = K(sqrt -d)  (d = 2 for the det-2 family).

For a binary lattice M of determinant d (Gram G in K^2, e1 = (1,0), Q(e1) = alpha) the map
  psi(v) = B(e1, v) + sqrt(-d) * det(e1, v)      satisfies   N_{E/K}(psi(v)) = alpha * Q(v),
so gamma in Q(M)  <=>  gamma*alpha = N(z) for some z in psi(M).  Passing to the O_E-ideal
A = psi(M) O_E (maximal order; PARI bnfinit, GRH), gamma is represented by the O_E-hull only if some
integral ideal B with B*conj(B) = (gamma) * n(A) lies in the class [A]^{-1}, n(A) = (alpha)/N(A).
  n(c) := #{square classes [gamma] of indecomposables : some B with B conj(B) = (gamma) lies in c}
is computed for every class c in Cl(O_E) (a fixed ideal shift n(A) only relabels classes).
For a given lattice the class of A^{-1} * (shift) is reported, so n(c_M) can be compared with the
EXACT s_sq(M).  n(c_M) >= s_sq(M) always (the O_E-hull contains M); equality = model exact.
Usage: python3 scripts/r2_classes.py D [d]          (prints the class distribution)
       import: class_distribution(D, d), lattice_class(...)"""
import os, sys
from fractions import Fraction as F
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from rqf import QF
from complete import FieldData

try:
    import cypari2
    pari = cypari2.Pari()
    pari.allocatemem(2 * 10**9)
except Exception:                       # pragma: no cover
    pari = None


class EData:
    """E = Q(sqrt D, sqrt -d) as an absolute PARI bnf in the variable y, with sqrt D, sqrt -d,
    and complex conjugation rho (fixes sqrt D, sends sqrt -d to -sqrt -d)."""

    def __init__(self, D, d):
        self.D, self.d = D, d
        pol = pari(f"subst(polcompositum(x^2-{D}, x^2+{d})[1], x, y)")
        self.bnf = pari.bnfinit(pol, 1)
        self.sqD = pari(f"(nf)->nfroots(nf, x^2-{D})[1]")(self.bnf)
        self.sqm = pari(f"(nf)->nfroots(nf, x^2+{d})[1]")(self.bnf)
        auts = pari.nfgaloisconj(self.bnf)
        ga = pari("(nf,a,z)->nfgaloisapply(nf,a,z)")
        self.rho = None
        for a in auts:
            if ga(self.bnf, a, self.sqD) == self.sqD and ga(self.bnf, a, self.sqm) == -self.sqm:
                self.rho = a
        assert self.rho is not None
        self._ga = ga
        self.cyc = [int(c) for c in self.bnf.bnf_get_cyc()]
        self._cls = pari("(b,I)->bnfisprincipal(b,I,0)")

    def elt(self, q, sqm_coeff=0):
        """K-element q (QF) + sqm_coeff * sqrt(-d) as an element of E."""
        return pari(str(q.a)) + pari(str(q.b)) * self.sqD + sqm_coeff * self.sqm

    def ecoef(self, x, y):
        """x + y sqrt(-d) for x, y in K (QF)."""
        return self.elt(x) + self.elt(y) * self.sqm

    def cls(self, I):
        return tuple(int(t) for t in self._cls(self.bnf, I))

    def half_norm_ideals(self, g):
        """All integral ideals B of O_E with B * rho(B) = g O_E (g in K)."""
        fa = pari.idealfactor(self.bnf, g)
        primes = [fa[0][i] for i in range(len(fa[0]))]
        exps = [int(fa[1][i]) for i in range(len(fa[1]))]
        done, choices = set(), []
        for i, P in enumerate(primes):
            if i in done:
                continue
            rP = self._ga(self.bnf, self.rho, P)
            j = next(k for k, Q in enumerate(primes)
                     if pari.idealhnf(self.bnf, Q) == pari.idealhnf(self.bnf, rP))
            done |= {i, j}
            e = exps[i]
            if j == i:
                if e % 2:
                    return []
                choices.append([pari.idealpow(self.bnf, P, e // 2)])
            else:
                choices.append([pari.idealmul(self.bnf, pari.idealpow(self.bnf, P, f),
                                              pari.idealpow(self.bnf, primes[j], e - f))
                                for f in range(e + 1)])
        out = [pari.idealhnf(self.bnf, 1)]
        for ch in choices:
            out = [pari.idealmul(self.bnf, B, C) for B in out for C in ch]
        return out


def class_distribution(D, d=2, ed=None, fd=None):
    """{class: set of square-class ids} over all indecomposables mod (O^x)^2, and per rep the
    set of classes of B with B rho(B) = (gamma)."""
    fd = fd or FieldData(D)
    ed = ed or EData(D, d)
    per_rep = []
    dist = {}
    for i, g in enumerate(fd.reps):
        cl = {ed.cls(B) for B in ed.half_norm_ideals(ed.elt(g))}
        per_rep.append(cl)
        for c in cl:
            dist.setdefault(c, set()).add(fd.class_of[i])
    return dist, per_rep, ed, fd


def lattice_class(ed, G, basis):
    """Class in Cl(O_E) of A = psi(M) O_E for the lattice with Z-basis `basis` (K^2 vectors)
    in the frame Gram G (psi uses e1 = (1, 0)), and the norm shift n(A) = (alpha)/N(A) as an ideal
    of E.  Returns the class c such that gamma in Q(M_hull) <=> some B in c has B rho(B) =
    (gamma) * n(A)  -- i.e. c = [A]^{-1}."""
    al, b = G[0][0], G[0][1]
    gens = []
    for (x, y) in basis:
        X = al * x + b * y
        gens.append(ed.ecoef(X, y))
    A = pari.idealhnf(ed.bnf, gens[0])
    for z in gens[1:]:
        A = pari.idealadd(ed.bnf, A, z)
    Ainv = pari.idealinv(ed.bnf, A)
    return ed.cls(Ainv), A


def predicted_for_lattice(D, G, basis, d=2):
    """n(c_M) for the lattice (allowing the shift n(A)): #square classes [gamma] such that some
    integral B with B rho(B) = (gamma) n(A) lies in [A]^{-1}."""
    fd = FieldData(D)
    ed = EData(D, d)
    cM, A = lattice_class(ed, G, basis)
    alpha = ed.elt(G[0][0])
    # n(A) O_E = alpha * A^{-1} * rho(A)^{-1}  (an ideal of E coming from K)
    rA = pari("(nf,a,I)->nfgaloisapply(nf,a,I)")(ed.bnf, ed.rho, A)
    nA = pari.idealmul(ed.bnf, alpha, pari.idealinv(ed.bnf, pari.idealmul(ed.bnf, A, rA)))
    hits = set()
    for i, g in enumerate(fd.reps):
        target = pari.idealmul(ed.bnf, ed.elt(g), nA)
        # ideals B with B rho(B) = target: factor target directly
        Bs = _half_norm_of_ideal(ed, target)
        if any(ed.cls(B) == cM for B in Bs):
            hits.add(fd.class_of[i])
    return len(hits), cM, ed


def _half_norm_of_ideal(ed, I):
    fa = pari.idealfactor(ed.bnf, I)
    primes = [fa[0][i] for i in range(len(fa[0]))]
    exps = [int(fa[1][i]) for i in range(len(fa[1]))]
    done, choices = set(), []
    for i, P in enumerate(primes):
        if i in done:
            continue
        rP = ed._ga(ed.bnf, ed.rho, P)
        j = next(k for k, Q in enumerate(primes)
                 if pari.idealhnf(ed.bnf, Q) == pari.idealhnf(ed.bnf, rP))
        done |= {i, j}
        e = exps[i]
        if j == i:
            if e % 2:
                return []
            choices.append([pari.idealpow(ed.bnf, P, e // 2)])
        else:
            if exps[j] != e:
                return []
            choices.append([pari.idealmul(ed.bnf, pari.idealpow(ed.bnf, P, f),
                                          pari.idealpow(ed.bnf, primes[j], e - f))
                            for f in range(e + 1)])
    out = [pari.idealhnf(ed.bnf, 1)]
    for ch in choices:
        out = [pari.idealmul(ed.bnf, B, C) for B in out for C in ch]
    return out


if __name__ == "__main__":
    D = int(sys.argv[1])
    d = int(sys.argv[2]) if len(sys.argv) > 2 else 2
    dist, per_rep, ed, fd = class_distribution(D, d)
    counts = sorted(((len(v), c) for c, v in dist.items()), reverse=True)
    print(f"D={D} d={d}: Cl(E) = {ed.cyc} (h = {eval('*'.join(map(str, ed.cyc)) or '1')}), "
          f"kappa_sq = {len(fd.class_reps)}, classes hit = {len(dist)}")
    print("  top n(c):", [n for n, _ in counts[:12]])


def e_to_K2(ed, z):
    """z in E -> (X, Y) in K^2 (QF) with z = X + Y sqrt(-d)."""
    D = ed.D
    if not hasattr(ed, "_Binv"):
        M = pari.matconcat([pari.nfalgtobasis(ed.bnf, b) for b in
                            (1, ed.sqD, ed.sqm, ed.sqD * ed.sqm)])
        ed._Binv = pari.matsolve(M, pari.matid(4))
    c = ed._Binv * pari.nfalgtobasis(ed.bnf, z)
    fr = [F(str(c[i])) for i in range(4)]
    return QF(D, fr[0], fr[1]), QF(D, fr[2], fr[3])


def lattice_from_ideal(ed, A, n):
    """The binary lattice (A, Q(z) = N_{E/K}(z)/n) in coordinates z = X + Y sqrt(-d):
    Gram G = diag(1/n, d/n), Z-basis = PARI's Z-basis of A.  n: a totally positive generator
    (QF) of N_{E/K}(A)."""
    D, d = ed.D, ed.d
    H = pari.idealhnf(ed.bnf, A)
    col = pari("(nf,H,j)->nfbasistoalg(nf, H[,j])")
    basis = [e_to_K2(ed, col(ed.bnf, H, j + 1)) for j in range(4)]
    G = [[n.inv(), QF(D, 0)], [QF(D, 0), QF(D, d) * n.inv()]]
    return G, basis


def half_norm_groups(ed, g):
    """(n0, Bs): n0 = HNF key of the product of the self-conjugate primes of E dividing g with odd
    exponent (= primes of K inert in E, odd exponent in g), Bs = all integral B with
    B rho(B) = g * n0^{-1} (n0 viewed in K)."""
    fa = pari.idealfactor(ed.bnf, g)
    primes = [fa[0][i] for i in range(len(fa[0]))]
    exps = [int(fa[1][i]) for i in range(len(fa[1]))]
    done, choices, n0 = set(), [], []
    for i, P in enumerate(primes):
        if i in done:
            continue
        rP = ed._ga(ed.bnf, ed.rho, P)
        j = next(k for k, Q in enumerate(primes)
                 if pari.idealhnf(ed.bnf, Q) == pari.idealhnf(ed.bnf, rP))
        done |= {i, j}
        e = exps[i]
        if j == i:
            if e % 2:
                n0.append(P)
            choices.append([pari.idealpow(ed.bnf, P, e // 2)])
        else:
            choices.append([pari.idealmul(ed.bnf, pari.idealpow(ed.bnf, P, f),
                                          pari.idealpow(ed.bnf, primes[j], e - f))
                            for f in range(e + 1)])
    out = [pari.idealhnf(ed.bnf, 1)]
    for ch in choices:
        out = [pari.idealmul(ed.bnf, B, C) for B in out for C in ch]
    N0 = pari.idealhnf(ed.bnf, 1)
    for P in n0:
        N0 = pari.idealmul(ed.bnf, N0, P)
    return str(N0), out


def class_distribution_groups(D, d=2, ed=None, fd=None):
    """{n0: {class: set of square-class ids}}.  Group n0 = 1 (the unit ideal) is class_distribution.
    max over (n0, c) of the counts = max s_sq over ALL integral lattices with det in d (K^x)^2
    (docs/10_R2.md, Theorem 10.2 / Cor. 10.3), when 2 eps is a square (d = 2)."""
    fd = fd or FieldData(D)
    ed = ed or EData(D, d)
    groups = {}
    for i, g in enumerate(fd.reps):
        key, Bs = half_norm_groups(ed, ed.elt(g))
        dist = groups.setdefault(key, {})
        for B in Bs:
            dist.setdefault(ed.cls(B), set()).add(fd.class_of[i])
    return groups, ed, fd


def classical_sublattices(D, G, basis):
    """Maximal classically integral sublattices M of a norm-integral binary O_K-lattice N
    (Z-basis `basis`, Gram G; Q(N) in O_K, 2B(N,N) in O_K).  Every such M contains 2N, so M/2N is an
    O-stable subspace S of N/2N (F_2^4) with B(S,S) in O_K.  Returns Z-bases of the maximal ones."""
    from itertools import product as iprod
    from lattice import Bform, hnf
    from gram import omega
    w = omega(D)
    n = len(basis)
    Bij = [[Bform(G, basis[i], basis[j]) for j in range(n)] for i in range(n)]
    assert all((2 * Bij[i][j]).is_integral() and Bij[i][i].is_integral()
               for i in range(n) for j in range(n)), "N not norm-integral"
    # omega action in the Z-basis (integer matrix): coordinates via the dual of the trace form
    def coords(v):
        rows = [[b[0].a, b[0].b, b[1].a, b[1].b] for b in basis]
        t = [v[0].a, v[0].b, v[1].a, v[1].b]
        # solve c . rows = t over Q
        A = [[F(x) for x in r] for r in rows]
        M = [A[i][:] + [F(int(i == j)) for j in range(4)] for i in range(4)]
        for c in range(4):
            p = next(r for r in range(c, 4) if M[r][c] != 0)
            M[c], M[p] = M[p], M[c]
            pv = M[c][c]
            M[c] = [x / pv for x in M[c]]
            for r in range(4):
                if r != c and M[r][c] != 0:
                    f_ = M[r][c]
                    M[r] = [M[r][j] - f_ * M[c][j] for j in range(8)]
        inv = [r[4:] for r in M]
        c_ = [sum(F(t[k]) * inv[k][j] for k in range(4)) for j in range(4)]
        assert all(x.denominator == 1 for x in c_)
        return [int(x) for x in c_]
    W = [coords((w * b[0], w * b[1])) for b in basis]          # row i = omega*b_i in basis
    vecs = [c for c in iprod((0, 1), repeat=4)]

    def omega_mod2(c):
        return tuple(sum(c[i] * W[i][j] for i in range(4)) % 2 for j in range(4))

    def add2(c, d):
        return tuple((x + y) % 2 for x, y in zip(c, d))

    def pair_int(c, d):
        s = QF(D, 0)
        for i in range(4):
            if c[i]:
                for j in range(4):
                    if d[j]:
                        s = s + Bij[i][j]
        return s.is_integral()

    # all O-stable subspaces S (closure of spans) that are B-integral
    subs = set()
    frontier = [frozenset([(0, 0, 0, 0)])]
    subs.add(frontier[0])
    while frontier:
        nxt = []
        for S in frontier:
            for v in vecs:
                if v in S:
                    continue
                T = set(S)
                stack = [v]
                while stack:
                    x = stack.pop()
                    if x in T:
                        continue
                    T.add(x)
                    for y in list(T):
                        z = add2(x, y)
                        if z not in T:
                            stack.append(z)
                    stack.append(omega_mod2(x))
                T = frozenset(T)
                if T in subs:
                    continue
                if all(pair_int(x, y) for x in T for y in T):
                    subs.add(T)
                    nxt.append(T)
        frontier = nxt
    maximal = [S for S in subs if not any(S < T for T in subs)]
    out = []
    for S in maximal:
        gens = [(2 * b[0], 2 * b[1]) for b in basis]
        for c in S:
            if any(c):
                gens.append((sum((basis[i][0] * c[i] for i in range(4)), QF(D, 0)),
                             sum((basis[i][1] * c[i] for i in range(4)), QF(D, 0))))
        rows = [[g[0].a, g[0].b, g[1].a, g[1].b] for g in gens]
        den = 1
        for r in rows:
            for x in r:
                den = den * F(x).denominator // __import__("math").gcd(den, F(x).denominator)
        H = hnf([[int(F(x) * den) for x in r] for r in rows])
        out.append([(QF(D, F(h[0], den), F(h[1], den)), QF(D, F(h[2], den), F(h[3], den)))
                    for h in H])
    return out
