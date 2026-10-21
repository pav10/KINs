#!/usr/bin/env python3
"""R3 (docs/11_R3.md): the pigeonhole bound  S(K,2) >= max_c n(c) >= 2 N_E(K) / h(-8D)  for
K = Q(sqrt D), D = 3 (mod 4) squarefree with (H) 2*eps in K^{x2}, computed WITHOUT the class group of E:
  * indecomposables mod <eps> from the CF of sqrt D: beta_{j,t} = alpha_{j-1} + t alpha_j, j even,
    0 <= t < a_{j+1}, alpha_k = p_k + q_k sqrt D, with norm Q_j + 2 t P_{j+1} - t^2 Q_{j+1}
    (x_k = (P_k + sqrt D)/Q_k the complete quotients);
  * square classes mod (K^x)^2: under (H), [beta eps] = [2 beta]; classes are separated by the
    squarefree part of the norm and by quadratic characters at degree-1 primes of norm > 4D (both are
    square-class invariants), and same-signature pairs are compared exactly (rqf.same_square_class),
    so the count is EXACT;
  * gamma O_K is a norm of an integral O_E-ideal  <=>  N(gamma) is represented by x^2 + 2y^2
    (indecomposables are primitive; docs/11_R3.md, Lemma 11.2);
  * h(-8D) = class number of Q(sqrt -2D) (PARI qfbclassno(.,1); reduced-form count as a cross-check).
Usage: python3 scripts/r3_tools.py scan out.jsonl Dmin Dmax [i/m]   (primes D = 3 mod 4: nu*, ...)
       python3 scripts/r3_tools.py proxy out.jsonl Dmin Dmax [i/m]  (CF only: s, iota)
       python3 scripts/r3_tools.py list out.jsonl Dfile [i/m]       (scan the primes listed in Dfile)
       python3 scripts/r3_tools.py field D"""
import os, sys, json, time
from math import isqrt, gcd
from fractions import Fraction as F
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "src"))
from rqf import QF, same_square_class


# ----------------------------------------------------------------------------
# continued fraction of sqrt D and the indecomposables
# ----------------------------------------------------------------------------
def cf_sqrt(D):
    """(a, P, Q): x_k = (P_k + sqrt D)/Q_k, a_k = floor(x_k), k = 0..s (s = period, Q_s = 1)."""
    r = isqrt(D)
    a, P, Q = [r], [0], [1]
    while True:
        p = a[-1] * Q[-1] - P[-1]
        q = (D - p * p) // Q[-1]
        P.append(p)
        Q.append(q)
        a.append((p + r) // q)
        if q == 1:
            return a, P, Q


def has_H(D, cf=None):
    """(H) 2*eps in K^{x2}  <=>  x^2 - D y^2 = +-2 solvable  <=>  some Q_k = 2 (|2| < sqrt D, Lagrange).
    (If mu^2 = 2 eps then N(mu) = +-2; conversely N(mu) = +-2 gives mu^2/2 a totally positive unit eps^j
    with j odd, since 2 is not a square in K.)"""
    a, P, Q = cf or cf_sqrt(D)
    return 2 in Q


def indec_data(D, cf=None):
    """[(j, t, N)] for the indecomposables mod <eps> (D = 3 mod 4, so s is even): iota entries."""
    a, P, Q = cf or cf_sqrt(D)
    s = len(a) - 1
    out = []
    for j in range(0, s, 2):
        Qj, Pn, Qn = Q[j], P[j + 1], Q[j + 1]
        for t in range(a[j + 1]):
            out.append((j, t, Qj + 2 * t * Pn - t * t * Qn))
    return out


def convergents_big(D, cf):
    """p_k, q_k for k = -1..s-1 as lists indexed k+1."""
    a = cf[0]
    p, q = [1, a[0]], [0, 1]
    for k in range(1, len(a) - 1):
        p.append(a[k] * p[-1] + p[-2])
        q.append(a[k] * q[-1] + q[-2])
    return p, q


def element(D, conv, j, t):
    """beta_{j,t} = alpha_{j-1} + t alpha_j as (x, y)."""
    p, q = conv
    return p[j] + t * p[j + 1], q[j] + t * q[j + 1]


# ----------------------------------------------------------------------------
# elementary number theory
# ----------------------------------------------------------------------------
class Sieve:
    """Smallest prime factor up to n (array-backed)."""

    def __init__(self, n):
        from array import array
        self.n = n
        spf = array("I", range(n + 1))
        small = [i for i in range(2, isqrt(n) + 1) if all(i % q for q in range(2, isqrt(i) + 1))]
        for i in reversed(small):                 # smallest prime written last
            cnt = len(range(i * i, n + 1, i))
            spf[i * i::i] = array("I", [i]) * cnt
        self.spf = spf

    def factor(self, m):
        out = {}
        while m > 1:
            p = self.spf[m]
            m //= p
            out[p] = out.get(p, 0) + 1
        return out

    def is_prime(self, m):
        return m > 1 and self.spf[m] == m


def rep_x2_2y2(fac):
    """n > 0 is represented by x^2 + 2y^2  <=>  every prime p = 5, 7 (mod 8) has even exponent."""
    return all(e % 2 == 0 for p, e in fac.items() if p % 8 in (5, 7))


def sqfree(fac):
    out = 1
    for p, e in fac.items():
        if e % 2:
            out *= p
    return out


def _is_prime(n):
    if n < 2:
        return False
    for p in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        if n % p == 0:
            return n == p
    d, r = n - 1, 0
    while d % 2 == 0:
        d //= 2
        r += 1
    for b in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        x = pow(b, d, n)
        if x in (1, n - 1):
            continue
        for _ in range(r - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
    return True


def char_primes(D, m=8):
    """m primes l = 3 (mod 4), l > 4D, (D/l) = 1, with a root r of r^2 = D (mod l)."""
    out, l = [], 4 * D + 1
    while len(out) < m:
        l += 1
        if l % 4 == 3 and _is_prime(l) and pow(D, (l - 1) // 2, l) == 1:
            out.append((l, pow(D, (l + 1) // 4, l)))
    return out


def h_neg(N):
    """Class number h(-N) of primitive forms of discriminant -N < 0 (N = 0, 3 mod 4), by counting
    reduced forms (a, b, c): |b| <= a <= c, b >= 0 if |b| = a or a = c (reference implementation)."""
    h = 0
    for a in range(1, isqrt(N // 3) + 1):
        for b in range(-a + 1, a + 1):
            if (b - N) % 2:
                continue
            num = b * b + N
            if num % (4 * a):
                continue
            c = num // (4 * a)
            if c < a or (c == a and b < 0) or gcd(gcd(a, abs(b)), c) != 1:
                continue
            h += 1
    return h


_pari = None


def pari():
    global _pari
    if _pari is None:
        import cypari2
        _pari = cypari2.Pari()
    return _pari


def h8(D):
    """h(-8D) via PARI qfbclassno(., 1) (Euler products + functional equation)."""
    return int(pari()(f"qfbclassno({-8 * D}, 1)"))


# ----------------------------------------------------------------------------
# exact square classes and N_E
# ----------------------------------------------------------------------------
def field(D, sv=None, nchar=8, need_reps=False):
    """Exact kappa_sq (square classes of indecomposables mod (K^x)^2), N_E (those whose ideal is a norm
    from E = K(sqrt -2)) for D = 3 (mod 4) squarefree with (H).  Classes of the 2*iota elements
    {beta, 2 beta} (= {beta, beta*eps} mod squares under (H))."""
    cf = cf_sqrt(D)
    assert D % 4 == 3 and has_H(D, cf)
    a, P, Q = cf
    s = len(a) - 1
    sv = sv or Sieve(4 * D)
    chars = char_primes(D, nchar)
    # alpha_k mod l for each character prime, k = -1..s-1 (index k+1)
    res = []
    for l, r in chars:
        pm, qm, pc, qc = 1, 0, a[0] % l, 1
        al = [1, (pc + qc * r) % l]
        for k in range(1, s):
            pm, qm, pc, qc = pc, qc, (a[k] * pc + pm) % l, (a[k] * qc + qm) % l
            al.append((pc + qc * r) % l)
        res.append((l, al))
    two_chi = tuple(1 if l % 8 in (1, 7) else -1 for l, _ in chars)
    ind = indec_data(D, cf)
    conv = None
    buckets = {}          # signature -> list of (class id, (j, t, twice))
    classes = []          # per class: (j, t, twice, N, inE) of its first member
    members = []          # per class: all members (j, t, twice, N)
    for (j, t, N) in ind:
        fac = sv.factor(N)
        sq = sqfree(fac)
        inE = rep_x2_2y2(fac)
        chi = tuple(1 if pow((al[j] + t * al[j + 1]) % l, (l - 1) // 2, l) == 1 else -1
                    for (l, al) in res)
        for twice in (0, 1):
            sig = (sq, chi if not twice else tuple(c * d for c, d in zip(chi, two_chi)))
            lst = buckets.setdefault(sig, [])
            hit = False
            if lst:
                if conv is None:
                    conv = convergents_big(D, cf)
                x, y = element(D, conv, j, t)
                me = QF(D, x << twice, y << twice)
                for cid, (j2, t2, tw2) in lst:
                    x2, y2 = element(D, conv, j2, t2)
                    if same_square_class(me, QF(D, x2 << tw2, y2 << tw2)):
                        hit = True
                        break
            if hit:
                members[cid].append((j, t, twice, N))
            else:
                lst.append((len(classes), (j, t, twice)))
                classes.append((j, t, twice, N, inE))
                members.append([(j, t, twice, N)])
    kappa = len(classes)
    nE = sum(1 for c in classes if c[4])
    iota = len(ind)
    iotaE = sum(1 for (j, t, N) in ind if rep_x2_2y2(sv.factor(N)))
    rec = dict(D=D, s=s, iota=iota, iotaE=iotaE, kappa=kappa, nE=nE,
               sigma_odd=sum(a[1:s:2]), sigma_even=sum(a[2:s + 1:2]))
    if need_reps:
        if conv is None:
            conv = convergents_big(D, cf)
        rec["reps"] = [(element(D, conv, j, t), twice, N, inE) for (j, t, twice, N, inE) in classes]
        rec["members"] = members
    return rec


def zagier_T(D, sv):
    """T(D) = sum_{|j| < sqrt D} d(D - j^2)."""
    tot, r = 0, isqrt(D)
    for j in range(-r, r + 1):
        m = D - j * j
        if m > 0:
            dd = 1
            for e in sv.factor(m).values():
                dd *= e + 1
            tot += dd
    return tot


def scan(out, lo, hi, shard=(0, 1), Dlist=None):
    """Exact record per prime p = 3 (mod 4) in [lo, hi] (or in Dlist): nu* = max_f nu(f) = max_c n(c)
    (Theorem 11.3), sum and support of nu, kappa_sq, N_E, h(-8p), h(-p), h(K) (GRH), T(p)."""
    done = set()
    if os.path.exists(out):
        done = {json.loads(l)["D"] for l in open(out)}
    sv = Sieve(4 * (max(Dlist) if Dlist else hi))
    Ds = Dlist or [D for D in range(lo, hi + 1) if D % 4 == 3 and sv.is_prime(D)]
    Ds = [D for k, D in enumerate(Ds) if k % shard[1] == shard[0] and D not in done]
    P = pari()
    with open(out, "a") as f:
        for D in Ds:
            t0 = time.time()
            rec = field(D, sv, need_reps=True)
            nu = nu_distribution(D, rec, sv)
            del rec["reps"], rec["members"]
            rec["h8"] = h8(D)
            rec["nu_max"] = max(nu.values()) if nu else 0
            rec["nu_sum"] = sum(nu.values())
            rec["nu_supp"] = len(nu)
            rec["best_form"] = list(max(nu, key=lambda g: (nu[g], g))) if nu else None
            rec["hm"] = int(P(f"qfbclassno({-D}, 1)")) if D > 3 else 1      # h(-D)
            rec["hK"] = int(P(f"quadclassunit({4 * D})[1]"))               # GRH (informational)
            rec["T"] = zagier_T(D, sv)
            rec["sec"] = round(time.time() - t0, 2)
            f.write(json.dumps(rec) + "\n")
            f.flush()



# ----------------------------------------------------------------------------
# positive definite binary forms of discriminant -8D (pure Python), and the elementary model of n(c)
# ----------------------------------------------------------------------------
def qf_reduce(f):
    a, b, c = f
    while True:
        if c < a:
            a, b, c = c, -b, a
        elif b > a or b <= -a:
            k = (a - b) // (2 * a)            # b + 2ka in (-a, a]
            b2 = b + 2 * k * a
            c = (b2 * b2 - (b * b - 4 * a * c)) // (4 * a)
            b = b2
        else:
            if a == c and b < 0:
                b = -b
            return (a, b, c)


def _xgcd(a, b):
    x0, x1, y0, y1 = 1, 0, 0, 1
    while b:
        q, a, b = a // b, b, a % b
        x0, x1 = x1, x0 - q * x1
        y0, y1 = y1, y0 - q * y1
    return a, x0, y0


def qf_compose(f1, f2):
    """Gauss composition (Cohen, Alg. 5.4.7), reduced."""
    (a1, b1, c1), (a2, b2, c2) = f1, f2
    if a1 > a2:
        (a1, b1, c1), (a2, b2, c2) = (a2, b2, c2), (a1, b1, c1)
    s = (b1 + b2) // 2
    n = b2 - s
    if a2 % a1 == 0:
        y1, d = 0, a1
    else:
        d, u, v = _xgcd(a2, a1)
        y1 = u
    if s % d == 0:
        y2, x2, d1 = -1, 0, d
    else:
        d1, x2, y2 = _xgcd(s, d)
        y2 = -y2
    v1, v2 = a1 // d1, a2 // d1
    r = (y1 * y2 * n - x2 * c2) % v1
    b3 = b2 + 2 * v2 * r
    a3 = v1 * v2
    disc = b1 * b1 - 4 * a1 * c1
    c3 = (b3 * b3 - disc) // (4 * a3)
    return qf_reduce((a3, b3, c3))


def qf_inv(f):
    return qf_reduce((f[0], -f[1], f[2]))


def qf_pow(f, e, disc):
    one = qf_reduce((1, disc % 2, (disc % 2 - disc) // 4))
    if e < 0:
        f, e = qf_inv(f), -e
    out = one
    while e:
        if e & 1:
            out = qf_compose(out, f)
        f = qf_compose(f, f)
        e >>= 1
    return out


def _sqrt_mod(n, p):
    n %= p
    if n == 0:
        return 0
    if p % 4 == 3:
        return pow(n, (p + 1) // 4, p)
    q, s = p - 1, 0
    while q % 2 == 0:
        q //= 2
        s += 1
    z = 2
    while pow(z, (p - 1) // 2, p) != p - 1:
        z += 1
    m, c, t, r = s, pow(z, q, p), pow(n, q, p), pow(n, (q + 1) // 2, p)
    while t != 1:
        i, t2 = 0, t
        while t2 != 1:
            t2 = t2 * t2 % p
            i += 1
        b = pow(c, 1 << (m - i - 1), p)
        m, c, t, r = i, b * b % p, t * b * b % p, r * b % p
    return r


def prime_form(q, disc):
    """A form (q, b, c) of discriminant disc (even, -8D) for a prime q not inert."""
    if disc % q == 0 and q != 2:
        b = 0 if (disc // 4) % q == 0 else q
    elif q == 2:
        b = 0 if (disc // 4) % 2 == 0 else 2
    else:
        b = 2 * _sqrt_mod(disc // 4, q)
    if (b * b - disc) % (4 * q):
        b = q - b if (b % 2) != (disc % 2) else b
    assert (b * b - disc) % (4 * q) == 0, (q, b, disc)
    return qf_reduce((q, b, (b * b - disc) // (4 * q)))


def norm_classes(n, fac, disc, cache):
    """Classes of Cl(disc) containing an ideal of norm n (fac = factorization of n); None if none."""
    opts = [qf_pow((1, 0, -disc // 4), 0, disc)]
    for q, e in fac.items():
        leg = 0 if disc % q == 0 else (1 if q == 2 or pow(disc % q, (q - 1) // 2, q) == 1 else -1)
        if q == 2:
            leg = 0 if disc % 8 == 0 or disc % 4 == 0 else leg
        if leg == -1:
            if e % 2:
                return None
            continue
        if q not in cache:
            cache[q] = prime_form(q, disc)
        P = cache[q]
        exps = [e % 2] if leg == 0 else range(-e, e + 1, 2)
        powers = {qf_pow(P, x, disc) for x in exps}
        opts = list({qf_compose(o, p) for o in opts for p in powers})
    return set(opts)


def nu_distribution(D, rec=None, sv=None):
    """{reduced form f of disc -8D: nu(f)}, nu(f) = #{[gamma] in N_E : f represents N(gamma)}
    (docs/11_R3.md, Theorem 11.4: for D prime the multiset {nu(f)} over the principal genus equals
    {n(c)} over ker(Cl(E) -> Cl^+(K)))."""
    sv = sv or Sieve(4 * D)
    rec = rec or field(D, sv, need_reps=True)
    disc, cache, nu = -8 * D, {}, {}
    for (xy, twice, N, inE), mem in zip(rec["reps"], rec["members"]):
        if not inE:                       # inE is a square-class invariant (Lemma 11.2)
            continue
        cl = set()
        for (j, t, tw, n) in mem:         # a class is represented if ANY member's norm is
            n <<= 2 * tw
            cl |= norm_classes(n, sv.factor(n), disc, cache)
        for f in cl:
            nu[f] = nu.get(f, 0) + 1
    return nu



def proxy_scan(out, lo, hi, shard=(0, 1)):
    """CF-only pass (no factoring): s, iota, sigma_even for primes p = 3 (mod 4) in [lo, hi].  Under the
    heuristic of docs/11_R3.md section 4 the kernel average of nu is ~ C iota / sqrt p, so iota/sqrt p
    preselects fields for the exact computation."""
    with open(out, "a") as f:
        k = 0
        for D in range(lo + (3 - lo) % 4, hi + 1, 4):
            if not _is_prime(D):
                continue
            k += 1
            if k % shard[1] != shard[0]:
                continue
            a, _, _ = cf_sqrt(D)
            s = len(a) - 1
            f.write(json.dumps(dict(D=D, s=s, iota=sum(a[1:s:2]), sigma_even=sum(a[2:s + 1:2]))) + "\n")


if __name__ == "__main__":
    if sys.argv[1] == "list":           # list out.jsonl Dfile [i/m]
        sh = tuple(map(int, sys.argv[4].split("/"))) if len(sys.argv) > 4 else (0, 1)
        scan(sys.argv[2], 0, 0, sh, [int(x) for x in open(sys.argv[3]).read().split()])
    elif sys.argv[1] in ("scan", "proxy"):
        sh = tuple(map(int, sys.argv[5].split("/"))) if len(sys.argv) > 5 else (0, 1)
        (scan if sys.argv[1] == "scan" else proxy_scan)(sys.argv[2], int(sys.argv[3]), int(sys.argv[4]), sh)
    else:
        D = int(sys.argv[2])
        rec = field(D)
        rec["h8"] = h8(D)
        print(rec, "R =", 2 * rec["nE"] / rec["h8"])
