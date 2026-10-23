#!/usr/bin/env python3
"""R4 checks (docs/12_R4.md):
  prop5 [nmax]      Prop 5 family: D = n^2 + 2 (n even, squarefree), M = O^2 with Gram
                    [[g_t, t mu], [t mu, g_{t+1}]], t = n/2 - 1; exact checks + trace/second-place cosines
  reg Dmax          log eps >= (s - 2) log phi + log(sqrt D - 1) for all squarefree D <= Dmax
  cor213 tmax       proportion of indecomposables of Q(sqrt(t^2+1)) (t odd) whose norm is a sum of two squares
"""
import os, sys, math
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "src"))
sys.path.insert(0, HERE)
from rqf import QF, same_square_class, fundamental_unit, cf_quadratic
from lattice import indecomposable_test

PHI = (1 + 5 ** 0.5) / 2


def sqf(D):
    return all(D % (p * p) for p in range(2, math.isqrt(D) + 1))


def prop5(n):
    """The Prop 5 lattice for even n >= 6 with n^2 + 2 squarefree; returns the checks and cosines."""
    D = n * n + 2
    mu, t = QF(D, n, 1), n // 2 - 1
    g1, g2, b = QF(D, 1) + mu * t, QF(D, 1) + mu * (t + 1), mu * t
    det = g1 * g2 - b * b
    is_ind, _ = indecomposable_test(D)
    ok = dict(integral=all(x.is_integral() for x in (g1, g2, b)), posdef=g1.is_tot_pos() and det.is_tot_pos(),
              indec=is_ind(g1) and is_ind(g2), distinct=not same_square_class(g1, g2),
              det_formula=det == QF(D, 1) + mu * (2 * t + 1) + fundamental_unit(D)[0] * (2 * t))
    cos_tr = float(b.trace()) / math.sqrt(float(g1.trace()) * float(g2.trace()))
    e2 = lambda x: float(x.a) - float(x.b) * math.sqrt(D)
    cos2 = e2(b) / math.sqrt(e2(g1) * e2(g2))
    return ok, cos_tr, cos2


def log_eps(D):
    e, _ = fundamental_unit(D)
    a, b = e.a, e.b                      # e = a + b sqrt D > 1, a ~ b sqrt D
    return math.log(a.numerator) - math.log(a.denominator) + math.log(2) - 1e-12 if a > 10 ** 12 \
        else math.log(float(a) + float(b) * math.sqrt(D))


def reg_check(Dmax):
    bad = []
    for D in range(2, Dmax + 1):
        if math.isqrt(D) ** 2 == D or not sqf(D):
            continue
        P, Q = (1, 2) if D % 4 == 1 else (0, 1)
        s = len(cf_quadratic(P, Q, D)[1])
        if log_eps(D) < (s - 2) * math.log(PHI) + math.log(math.sqrt(D) - 1) - 1e-9:
            bad.append(D)
    return bad


def cor213(tmax):
    from r3_tools import Sieve
    sv = Sieve(tmax * tmax + 2)
    s2 = lambda f: all(e % 2 == 0 for p, e in f.items() if p % 4 == 3)
    out = []
    for t in range(3, tmax, 2):
        D = t * t + 1
        if any(e > 1 for e in sv.factor(D).values()):
            continue
        out.append((t, sum(1 for u in range(-t + 1, t + 1) if s2(sv.factor(D - u * u))) / (2 * t)))
    return out


if __name__ == "__main__":
    what = sys.argv[1]
    if what == "prop5":
        nmax = int(sys.argv[2]) if len(sys.argv) > 2 else 154
        ns = [n for n in range(6, nmax + 1, 2)]
        good = [n for n in ns if sqf(n * n + 2)]
        print(f"even n in [6,{nmax}]: {len(ns)}, with n^2+2 squarefree: {len(good)}")
        for n in good:
            ok, c, c2 = prop5(n)
            if not all(ok.values()) or n in (6, 12, 30, 100) or n == good[-1]:
                print(n, ok, f"1-cos_tr = {1 - c:.4f} <= 2/(n-2) = {2 / (n - 2):.4f}", f"cos at sigma2 = {c2:.4f}")
    elif what == "reg":
        print("violations:", reg_check(int(sys.argv[2])))
    elif what == "cor213":
        rows = cor213(int(sys.argv[2]))
        lo = 3
        for hi in (30, 100, 300, 1000, 3000, 6001):
            b = [p for t, p in rows if lo <= t < hi]
            if b:
                print(f"t in [{lo},{hi}): {len(b)} fields, mean {sum(b) / len(b):.3f}, "
                      f"x sqrt(log hi^2) = {sum(b) / len(b) * math.sqrt(math.log(hi * hi)):.3f}")
            lo = hi
