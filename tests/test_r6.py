"""R6 (docs/14_R6.md): indecomposable squares, the window criterion, T(K) <= 3s, and the family
D = n^2 + 3 with T(K) -> infinity."""
import os, sys, random
from math import gcd, isqrt
from collections import Counter
import pytest
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "scripts"))
from rqf import QF, same_square_class
from r6_squares import window_test, square_primitive, T1, sqf, family_bound
from r3_tools import cf_sqrt

FIELDS = [D for D in range(2, 300) if D % 4 in (2, 3) and sqf(D)]


def test_window_criterion_matches_cf_test():
    """Lemma 14.1: for primitive alpha = x + y sqrt D >> 0 (y > 0), the window/witness test agrees
    with the CF orbit test of indecomposability (all such alpha with N(alpha) <= D and y < 30)."""
    from lattice import indecomposable_test
    n = 0
    for D in FIELDS[:40]:
        is_ind, _ = indecomposable_test(D)
        for y in range(1, 30):
            x0 = isqrt(D * y * y)
            for x in range(x0 + 1, isqrt(D * y * y + D) + 1):
                if gcd(x, y) == 1:
                    assert window_test(D, x, y)[0] == is_ind(QF(D, x, y)), (D, x, y)
                    n += 1
    assert n > 300


def test_small_norm_primitive_is_indecomposable():
    """Corollary: primitive alpha >> 0 with N(alpha) < 2 sqrt D is indecomposable."""
    from lattice import indecomposable_test
    random.seed(3)
    for D in FIELDS[10:40]:
        is_ind, _ = indecomposable_test(D)
        for y in range(1, 200):
            x = isqrt(D * y * y) + 1
            if gcd(x, y) == 1 and x * x - D * y * y < 2 * D ** 0.5:
                assert is_ind(QF(D, x, y)), (D, x, y)


def test_primitivity_criterion():
    """Lemma 14.3: xi^2 primitive <=> gcd(N(xi), 2D) = 1 (gcd(p, q) = 1)."""
    random.seed(4)
    for _ in range(20000):
        D = random.choice(FIELDS)
        p, q = random.randint(1, 10 ** 5), random.randint(1, 10 ** 5)
        if gcd(p, q) == 1:
            assert square_primitive(D, p, q) == square_primitive(D, p, q, abs(p * p - D * q * q))


def test_T_bounds_and_T1():
    """T1 = fibre of the class of 1 (exact square classes); T(K) <= 2s (R7: classical rank-1 per-face lemma); T1 <= s."""
    from complete import FieldData
    for D in FIELDS:
        fd = FieldData(D)
        c = Counter(fd.class_of)
        cid = next(fd.class_of[i] for i, r in enumerate(fd.reps) if same_square_class(r, QF(D, 1)))
        s = len(cf_sqrt(D)[0]) - 1
        t1 = T1(D)
        assert c[cid] == t1 and t1 <= s and max(c.values()) <= 2 * s, D


def test_family_n2_plus_3():
    """Theorem 14.2: D = n^2 + 3 (n even, 3 not | n, squarefree): xi^{2j} (xi = n + sqrt D) and its
    conjugate are indecomposable for 9^j < 2n, so T1 >= 1 + 2 floor(log_9(2n))."""
    for n in range(2, 400, 2):
        D = n * n + 3
        if n % 3 == 0 or not sqf(D):
            continue
        xi = QF(D, n, 1)
        assert xi.norm() == -3
        z = QF(D, 1)
        for j in range(1, 10):
            z = z * xi * xi
            if 9 ** j >= 2 * n:
                break
            assert window_test(D, int(z.a), int(z.b))[0]   # the conjugate follows by Galois symmetry
        assert T1(D) >= family_bound(n)
