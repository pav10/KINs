"""Experiments: Gram search (exact rank over K), realizability, sums of two squares."""
from fractions import Fraction as F
from math import isqrt
from rqf import QF
from indec import square_class_reps


def rank_over_K(M):
    """Exact rank of a matrix of QF entries over K = Q(sqrt D)."""
    if not M:
        return 0
    D = M[0][0].D
    rows = [row[:] for row in M]
    n, m = len(rows), len(rows[0])
    r = 0
    for col in range(m):
        piv = next((i for i in range(r, n) if not rows[i][col].is_zero()), None)
        if piv is None:
            continue
        rows[r], rows[piv] = rows[piv], rows[r]
        inv = rows[r][col].inv()
        rows[r] = [inv * e for e in rows[r]]
        for i in range(n):
            if i != r and not rows[i][col].is_zero():
                f = rows[i][col]
                rows[i] = [rows[i][j] - f * rows[r][j] for j in range(m)]
        r += 1
        if r == n:
            break
    return r


def omega(D):
    return QF(D, F(1, 2), F(1, 2)) if D % 4 == 1 else QF(D, 0, 1)


def valid_off_diagonal(D, a, b):
    """All c in O_K with a*b - c^2 totally non-negative (classical)."""
    w = omega(D)
    w1, w2 = w.emb()
    a1, a2 = a.emb(); b1, b2 = b.emb()
    B1 = max(0.0, a1 * b1) ** 0.5
    B2 = max(0.0, a2 * b2) ** 0.5
    dw = abs(w1 - w2)
    ymax = (B1 + B2) / dw
    out = []
    for y in range(int(-ymax) - 2, int(ymax) + 3):
        x1lo, x1hi = -B1 - y * w1, B1 - y * w1
        x2lo, x2hi = -B2 - y * w2, B2 - y * w2
        xlo = max(min(x1lo, x1hi), min(x2lo, x2hi))
        xhi = min(max(x1lo, x1hi), max(x2lo, x2hi))
        for x in range(int(xlo) - 2, int(xhi) + 3):
            c = QF(D, x, 0) + QF(D, y, 0) * w
            diff = a * b - c * c
            if diff.is_tot_nonneg():
                out.append(c)
    return out


def largest_gram(D, verbose=False, sizecap=12):
    """Max size of rank<=2 totally PSD Gram matrix on distinct indec. square classes."""
    sq = square_class_reps(D, complete=True)
    cache = {}

    def vc(a, b):
        key = (a.a, a.b, b.a, b.b)
        if key in cache:
            return cache[key]
        r = valid_off_diagonal(D, a, b)
        cache[key] = r
        return r

    best = [0]
    best_M = [None]

    def search(M, start):
        k = len(M)
        if k >= 4 and rank_over_K(M) > 2:
            return
        if k > best[0]:
            best[0] = k
            best_M[0] = [row[:] for row in M]
            if verbose:
                print(f"  size {k}")
        if k >= sizecap:
            return
        for i in range(start, len(sq)):
            vnew = sq[i]
            if k == 0:
                search([[vnew]], i + 1)
                continue
            lists = [vc(M[j][j], vnew) for j in range(k)]

            def bt(idx, cur):
                if idx == k:
                    newM = [row[:] + [cur[j]] for j, row in enumerate(M)]
                    newM.append(list(cur) + [vnew])
                    search(newM, i + 1)
                    return
                for c in lists[idx]:
                    ok = True
                    for j in range(idx):
                        det = (M[j][j] * M[idx][idx] * vnew
                               + QF(D, 2, 0) * M[j][idx] * cur[j] * c
                               - M[j][j] * c * c
                               - M[idx][idx] * cur[j] * cur[j]
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
    return best[0], best_M[0]


def realize(M):
    """
    Test realizability of abstract Gram M as V^T G V with G a 2x2 classical form
    over O_K and integral coordinate vectors. Returns (True, G, coords) or (False, None, None).
    """
    s = len(M)
    D = M[0][0].D
    # find pivot pair with invertible 2x2 minor
    for p in range(s):
        for q in range(p + 1, s):
            det = M[p][p] * M[q][q] - M[p][q] * M[q][p]
            if det.is_zero():
                continue
            G = [[M[p][p], M[p][q]], [M[q][p], M[q][q]]]
            di = det.inv()
            coords = []
            ok = True
            for i in range(s):
                x = (M[q][q] * M[p][i] - M[p][q] * M[q][i]) * di
                y = (M[p][p] * M[q][i] - M[q][p] * M[p][i]) * di
                if not (x.is_integral() and y.is_integral()):
                    ok = False
                    break
                coords.append((x, y))
            if ok:
                # verify Q(v_i) == M[i][i]
                for i in range(s):
                    x, y = coords[i]
                    val = x * x * G[0][0] + QF(D, 2, 0) * x * y * G[0][1] + y * y * G[1][1]
                    assert val == M[i][i]
                return True, G, coords
    return False, None, None


def is_sum_of_two_squares(D, alpha):
    """Is the totally positive alpha = x^2 + y^2 with x,y in O_K?"""
    if not alpha.is_tot_pos():
        return alpha.is_zero()
    w = omega(D)
    w1, w2 = w.emb()
    a1, a2 = alpha.emb()
    B1, B2 = a1 ** 0.5, a2 ** 0.5
    dw = abs(w1 - w2)
    ymax = (B1 + B2) / dw
    for y in range(0, int(ymax) + 3):          # x and -x equivalent: y>=0, and handle x sign in inner
        x1lo, x1hi = -B1 - y * w1, B1 - y * w1
        x2lo, x2hi = -B2 - y * w2, B2 - y * w2
        xlo = max(min(x1lo, x1hi), min(x2lo, x2hi))
        xhi = min(max(x1lo, x1hi), max(x2lo, x2hi))
        for x in range(int(xlo) - 2, int(xhi) + 3):
            c = QF(D, x, 0) + QF(D, y, 0) * w
            rem = alpha - c * c
            if rem.is_tot_nonneg() and rem.is_square_in_K():
                return True
    return False


def count_sos_indec(D):
    """How many indecomposable square classes are represented by x^2+y^2."""
    sq = square_class_reps(D, complete=True)
    return sum(1 for a in sq if is_sum_of_two_squares(D, a)), len(sq)


if __name__ == "__main__":
    print("=== validate max Gram sizes vs doc 6 ===")
    for D, exp in [(19, 4), (22, 4), (31, 4), (58, 6)]:
        sz, _ = largest_gram(D)
        print(f"D={D:3d}  max size={sz}  (doc6={exp})  {'OK' if sz == exp else 'MISMATCH'}")
