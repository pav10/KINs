"""Indecomposables (Dress-Scharlau / CF), unit normalization, square classes."""
from fractions import Fraction as F
from math import log
from rqf import QF, cf_quadratic, fundamental_unit, tp_unit, same_square_class


def cf_list(D, nterms):
    """Partial quotients of the CF used in doc3, plus delta."""
    if D % 4 == 1:
        pre, per = cf_quadratic(-1, 2, D)        # (sqrt D - 1)/2
        delta = QF(D, F(1, 2), F(1, 2))          # (sqrt D + 1)/2
    else:
        pre, per = cf_quadratic(0, 1, D)         # sqrt D
        delta = QF(D, 0, 1)
    cf = list(pre)
    L = len(per)
    i = 0
    while len(cf) < nterms:
        cf.append(per[(len(cf) - len(pre)) % L])
        i += 1
    return cf, delta, L


def convergents(cf, n):
    """(p_i, q_i) for i=0..n-1."""
    out = []
    pm1, qm1 = 1, 0
    p0, q0 = cf[0], 1
    out.append((p0, q0))
    for i in range(1, n):
        p = cf[i] * p0 + pm1
        q = cf[i] * q0 + qm1
        out.append((p, q))
        pm1, qm1, p0, q0 = p0, q0, p, q
    return out


def normalize_unit(x, eta):
    """
    Canonical rep of totally positive x in its tp-unit orbit.
    Fundamental domain: ratio emb1/emb2 in [1, eta1^2), tested EXACTLY:
      ratio(y) >= 1     <=>  y.b >= 0
      ratio(y) <  eta1^2 <=> ratio(y*eta.conj()) < 1  <=> (y*eta.conj()).b < 0
    """
    etac = eta.conj()
    for _ in range(100000):
        if x.b < 0:
            x = x * eta
        elif (x * etac).b >= 0:
            x = x * etac
        else:
            return x
    raise RuntimeError("unit normalization did not converge")


def indecomposables(D):
    """Orbit representatives of indecomposables (up to totally positive units)."""
    eta = tp_unit(D)
    nterms = None
    cf, delta, s = cf_list(D, 2 * 50 + 14)
    nterms = 2 * s + 12
    cf, delta, s = cf_list(D, nterms + 2)
    cvs = convergents(cf, nterms + 1)
    # alpha list: alpha[0]=1 ; alpha[k]=p_{k-1}+q_{k-1}*delta for k>=1
    alpha = [QF(D, 1, 0)]
    for (p, q) in cvs:
        alpha.append(QF(D, p, 0) + QF(D, q, 0) * delta)
    cand = []
    for i in range(0, 2 * s + 10, 2):
        for t in range(cf[i + 1] + 1):
            ait = alpha[i] + QF(D, t, 0) * alpha[i + 1]
            cand.append(ait)
    # normalize up to units and dedup (by exact equality of canonical rep)
    reps = []
    seen = set()
    for c in cand:
        if not c.is_tot_pos():
            continue
        cn = normalize_unit(c, eta)
        key = (cn.a, cn.b)
        if key not in seen:
            seen.add(key)
            reps.append(cn)
    return reps


def square_class_reps(D, complete=True):
    """Square classes of indecomposables. complete=True adds alpha*eta (doc).""" 
    eta = tp_unit(D)
    inds = indecomposables(D)
    pool = []
    for a in inds:
        pool.append(a)
        if complete:
            pool.append(a * eta)
    classes = []
    for x in pool:
        if not any(same_square_class(x, r) for r in classes):
            classes.append(x)
    return classes


if __name__ == "__main__":
    # validate square-class counts against doc 6 (complete)
    expect = {2: 2, 3: 2, 5: 1, 6: 4, 7: 4, 10: 5, 11: 6, 13: 3, 14: 4, 15: 2,
              17: 3, 19: 10, 21: 2, 22: 8, 23: 4, 26: 9, 29: 5, 30: 4, 31: 12,
              43: 22, 58: 17, 67: 34}
    ok = True
    for D, e in expect.items():
        n = len(square_class_reps(D, complete=True))
        flag = "" if n == e else f"  <-- expected {e}"
        if n != e:
            ok = False
        print(f"D={D:3d}  #square classes={n}{flag}")
    print("ALL MATCH" if ok else "MISMATCH")
