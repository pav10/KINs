#!/usr/bin/env python3
"""R6 (docs/14_R6.md): indecomposable squares and T(K) = max over square classes C of #(I(K) cap C)/(O^x)^2,
for K = Q(sqrt D), D = 2, 3 (mod 4) squarefree, O = Z[sqrt D].

window_test(D, x, y): exact test for alpha = x + y sqrt D >> 0 with gcd(x, y) = 1, y > 0 (Lemma 14.1):
  with m in [1, x], m*y = -1 (mod x), alpha is indecomposable iff m*(x - y sqrt D) < sqrt D, and then
  delta = (b + a/sqrt D)/2 with a = -m, b = (1 + m y)/x is a totally positive codifferent element with
  Tr(delta alpha) = b x + a y = 1 (the witness).  Exact: compare m^2 x^2 < D (1 + m y)^2.
T1(D): T for the square class of 1 = 1 + #{xi mod +-units : xi^2 indecomposable, xi non-unit}; the xi are
  the convergents alpha_k (|N xi| = Q_{k+1} <= sqrt D, Lagrange), so T1 costs one pass over the period.
Usage: python3 scripts/r6_squares.py T1 Dmin Dmax         (records of T1 over squarefree D = 2,3 mod 4)
       python3 scripts/r6_squares.py family nmax           (Theorem 14.2 family D = n^2 + 3)"""
import os, sys, math
from math import gcd, isqrt
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "src"))
sys.path.insert(0, HERE)
from r3_tools import cf_sqrt


def window_test(D, x, y):
    """alpha = x + y sqrt D (x > y sqrt D > 0, gcd(x, y) = 1): (indecomposable?, witness (a, b))."""
    assert y > 0 and x * x > D * y * y and gcd(x, y) == 1
    if x == 1:
        return True, (0, 1)
    m = (-pow(y, -1, x)) % x or x
    if m * m * x * x < D * (1 + m * y) ** 2:     # m sigma_2(alpha) < sqrt D
        a, b = -m, (1 + m * y) // x
        assert b * x + a * y == 1 and b > 0 and D * b * b > a * a
        return True, (a, b)
    return False, None


def square_primitive(D, p, q, n=None):
    """xi = p + q sqrt D (gcd(p, q) = 1); xi^2 = (p^2 + D q^2) + 2pq sqrt D is primitive iff
    gcd(N(xi), 2D) = 1 (Lemma 14.3); checked directly when n = |N(xi)| is not given."""
    if n is not None:
        return gcd(n, 2 * D) == 1
    return gcd(p * p + D * q * q, 2 * p * q) == 1


def T1(D, cf=None, details=False):
    a, P, Q = cf or cf_sqrt(D)
    s = len(a) - 1
    pm, qm, p, q = 1, 0, a[0], 1                 # alpha_{-1} = 1, alpha_0 = a0 + sqrt D
    hits = []
    for k in range(s):                           # alpha_k = p_k + q_k sqrt D, |N| = Q_{k+1}
        n = Q[k + 1]
        if 1 < n and n * n <= D and square_primitive(D, p, q, n):
            x, y = p * p + D * q * q, 2 * p * q
            ok, w = window_test(D, x, y)
            if ok:
                hits.append((k, n))
        pm, qm, p, q = p, q, a[k + 1] * p + pm, a[k + 1] * q + qm
    return (1 + len(hits), hits) if details else 1 + len(hits)


def sqf(D):
    r = isqrt(D)
    return r * r != D and all(D % (p * p) for p in range(2, r + 1))


def family_bound(n):
    """Theorem 14.2: 1 + 2 #{j >= 1 : 9^j < 2n}."""
    j, c = 0, 9
    while c < 2 * n:
        j, c = j + 1, c * 9
    return 1 + 2 * j


if __name__ == "__main__":
    if sys.argv[1] == "T1":
        lo, hi = int(sys.argv[2]), int(sys.argv[3])
        best = 0
        for D in range(lo, hi + 1):
            if D % 4 in (2, 3) and sqf(D):
                t = T1(D)
                if t > best:
                    best = t
                    print(D, "T1 =", t, "s =", len(cf_sqrt(D)[0]) - 1, flush=True)
    elif sys.argv[1] == "family":
        nmax = int(sys.argv[2])
        for n in range(2, nmax + 1, 2):
            D = n * n + 3
            if n % 3 == 0 or not sqf(D):
                continue
            t, hits = T1(D, details=True)
            b = family_bound(n)
            assert t >= b, (n, t, b)
            if n in (10, 100, 1000, 10000) or n == nmax or n == nmax - 1:
                print(f"n={n} D={D}: T1 = {t} >= bound {b}; norms of xi: {sorted(h[1] for h in hits)}", flush=True)
