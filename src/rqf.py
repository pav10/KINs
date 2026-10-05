"""
Pure-Python real quadratic field machinery (no Sage), exact arithmetic.

Element of Q(sqrt D) stored in the sqrt(D)-basis as (a, b) with a,b in Fraction:
    x = a + b*sqrt(D).
O_K = Z[sqrt D]            if D = 2,3 (mod 4)
    = Z[(1+sqrt D)/2]      if D = 1   (mod 4)
"""
from fractions import Fraction as F
from math import isqrt


def _is_rational_square(x: F):
    """If x>=0 is a perfect rational square return its sqrt (Fraction), else None."""
    if x < 0:
        return None
    n, d = x.numerator, x.denominator
    rn, rd = isqrt(n), isqrt(d)
    if rn * rn == n and rd * rd == d:
        return F(rn, rd)
    return None


class QF:
    __slots__ = ("D", "a", "b")

    def __init__(self, D, a, b=0):
        self.D = D
        self.a = F(a)
        self.b = F(b)

    # ---- arithmetic ----
    def __add__(s, o): return QF(s.D, s.a + o.a, s.b + o.b)
    def __sub__(s, o): return QF(s.D, s.a - o.a, s.b - o.b)
    def __neg__(s):    return QF(s.D, -s.a, -s.b)

    def __mul__(s, o):
        if isinstance(o, QF):
            return QF(s.D, s.a * o.a + s.b * o.b * s.D, s.a * o.b + s.b * o.a)
        return QF(s.D, s.a * o, s.b * o)  # scalar

    def __rmul__(s, o): return s.__mul__(o)

    def inv(s):
        n = s.norm()
        if n == 0:
            raise ZeroDivisionError
        return QF(s.D, s.a / n, -s.b / n)

    def __truediv__(s, o):
        o = o if isinstance(o, QF) else QF(s.D, o, 0)
        return s * o.inv()

    def __eq__(s, o):
        o = o if isinstance(o, QF) else QF(s.D, o, 0)
        return s.a == o.a and s.b == o.b

    def __hash__(s): return hash((s.a, s.b))
    def __repr__(s): return f"({s.a}+{s.b}*sqrt{s.D})"

    # ---- field data ----
    def conj(s):  return QF(s.D, s.a, -s.b)
    def norm(s):  return s.a * s.a - s.b * s.b * s.D
    def trace(s): return 2 * s.a
    def emb(s):
        r = s.D ** 0.5
        return (float(s.a) + float(s.b) * r, float(s.a) - float(s.b) * r)

    def is_zero(s): return s.a == 0 and s.b == 0

    def is_tot_pos(s):
        # both embeddings > 0  <=>  a>0 and norm>0
        return s.a > 0 and s.norm() > 0

    def is_tot_nonneg(s):
        if s.is_zero():
            return True
        return s.a >= 0 and s.norm() >= 0

    def is_integral(s):
        if s.D % 4 == 1:
            return (2 * s.b).denominator == 1 and (s.a - s.b).denominator == 1
        return s.a.denominator == 1 and s.b.denominator == 1

    def sqrt_in_K(s):
        """Return r with r*r == s if s is a square in K, else None."""
        N = s.norm()
        rn = _is_rational_square(N)
        if rn is None:
            return None
        for n in (rn, -rn):
            p2 = (s.a + n) / 2
            p = _is_rational_square(p2)
            if p is None:
                continue
            if s.D == 0:
                continue
            q2 = (s.a - n) / (2 * s.D)
            q = _is_rational_square(q2)
            if q is None:
                continue
            for ps in (p, -p):
                for qs in (q, -q):
                    if 2 * ps * qs == s.b and ps * ps + qs * qs * s.D == s.a:
                        return QF(s.D, ps, qs)
        return None

    def is_square_in_K(s):
        return s.sqrt_in_K() is not None


def same_square_class(x: QF, y: QF):
    """x ~ y iff x*y is a square in K (equiv. x/y a square)."""
    return (x * y).is_square_in_K()


# --------------------------------------------------------------------------
# Continued fraction of a quadratic irrational (P + sqrt D)/Q, exact.
# Returns (a0, period) where the expansion is [a0; period-repeated].
# --------------------------------------------------------------------------
def cf_quadratic(P, Q, D):
    """CF of (P + sqrt D)/Q with Q | (D - P^2). Detects the period."""
    assert (D - P * P) % Q == 0
    sqrtD = isqrt(D)
    seen = {}
    seq = []
    p, q = P, Q
    while (p, q) not in seen:
        seen[(p, q)] = len(seq)
        a = (p + sqrtD) // q if q > 0 else -((-p - sqrtD - 1) // (-q))
        # floor((p+sqrtD)/q):
        a = (p + sqrtD) // q
        seq.append(a)
        p = a * q - p
        q = (D - p * p) // q
    start = seen[(p, q)]
    return seq[:start], seq[start:]   # (pre-period, period)


def fundamental_unit_bruteforce(D, ubound=10 ** 7):
    """Reference implementation (slow): minimal solution of t^2 - D u^2 = +-4.
    Kept only for cross-validation of fundamental_unit()."""
    for u in range(1, ubound + 1):
        for sgn in (-4, 4):
            tt = D * u * u + sgn
            if tt <= 0:
                continue
            t = isqrt(tt)
            if t * t != tt:
                continue
            eps = QF(D, F(t, 2), F(u, 2))
            if not eps.is_integral():
                continue
            return eps, int(eps.norm())
    raise RuntimeError(f"no fundamental unit found for D={D} within u<={ubound}")


def fundamental_unit(D):
    """
    Fundamental unit eps > 1 of O_K (K = Q(sqrt D), D squarefree > 1) and N(eps).
    Via the continued fraction of omega_D (omega = sqrt D, or (1+sqrt D)/2 if D=1 mod 4):
    with p_k/q_k the convergents of omega, x_k = p_k - q_k*conj(omega) is a unit for
    k = s-1 (s = period length), and this is the fundamental unit; N = (-1)^s.
    We scan k and return the first x_k with |N(x_k)| = 1 (robust to off-by-one).
    """
    if D % 4 == 1:
        P, Q = 1, 2
        om_conj = QF(D, F(1, 2), F(-1, 2))
    else:
        P, Q = 0, 1
        om_conj = QF(D, 0, -1)
    pre, per = cf_quadratic(P, Q, D)
    cf = list(pre) + list(per) * 3
    pm1, qm1, p, q = 1, 0, cf[0], 1
    k = 0
    while True:
        x = QF(D, p, 0) - QF(D, q, 0) * om_conj
        n = x.norm()
        if abs(n) == 1 and x.emb()[0] > 1:
            return x, int(n)
        k += 1
        if k >= len(cf):
            raise RuntimeError(f"fundamental unit not found for D={D}")
        p, pm1 = cf[k] * p + pm1, p
        q, qm1 = cf[k] * q + qm1, q


def tp_unit(D):
    """Generator eta>1 of totally positive units."""
    eps, Neps = fundamental_unit(D)
    eta = eps if Neps > 0 else eps * eps
    assert eta.norm() == 1 and eta.is_tot_pos() and eta.emb()[0] > 1, (D, eta)
    return eta


if __name__ == "__main__":
    # sanity: fundamental units & their norms
    for D in [2, 3, 5, 6, 7, 11, 13, 43, 67]:
        eps, n = fundamental_unit(D)
        print(f"D={D:3d}  eps={eps}  N(eps)={n:+d}  emb={eps.emb()[0]:.4f}")
