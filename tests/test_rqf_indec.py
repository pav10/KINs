"""Ground-truth regression tests for field arithmetic and indecomposables."""
from math import isqrt
from fractions import Fraction as F
from rqf import QF, fundamental_unit, fundamental_unit_bruteforce, cf_quadratic, tp_unit, same_square_class
from indec import indecomposables, square_class_reps


def squarefree(D):
    r = isqrt(D)
    return r * r != D and all(D % (p * p) for p in range(2, r + 1))


def test_fundamental_unit_cf_matches_bruteforce():
    for D in range(2, 120):
        if squarefree(D):
            try:
                ref = fundamental_unit_bruteforce(D, 10 ** 5)
            except RuntimeError:          # unit too large for the brute-force reference
                continue
            assert fundamental_unit(D) == ref, D


def test_known_units():
    assert fundamental_unit(5)[0] == QF(5, F(1, 2), F(1, 2)) and fundamental_unit(5)[1] == -1
    assert fundamental_unit(19)[0] == QF(19, 170, 39)          # N=+1
    assert tp_unit(31) == QF(31, 1520, 273)


def test_blomer_kala_count():
    """#indecomposables mod totally positive units = u1+u3+..+u_{s-1} (s even), sum(u) (s odd),
    [u1..us] = period of omega_D. Verified D<300 at build time."""
    for D in range(2, 120):
        if not squarefree(D):
            continue
        P, Q = (1, 2) if D % 4 == 1 else (0, 1)
        _, per = cf_quadratic(P, Q, D)
        s = len(per)
        md = sum(per[i] for i in range(0, s - 1, 2)) if s % 2 == 0 else sum(per)
        assert len(indecomposables(D)) == md, D


def test_square_class_counts():
    """kappa_sq(K) = # square classes of indecomposables (seeding alpha and alpha*eta)."""
    expect = {2: 2, 3: 2, 5: 1, 6: 4, 7: 4, 10: 5, 11: 6, 13: 3, 14: 4, 15: 2, 17: 3,
              19: 10, 21: 2, 22: 8, 23: 4, 26: 9, 29: 5, 30: 4, 31: 12, 43: 22, 58: 17, 67: 34}
    for D, e in expect.items():
        assert len(square_class_reps(D)) == e, D


def test_tp_units_are_indecomposable():
    """Every totally positive unit is indecomposable (eps = b+c, b,c >> 0 gives
    b/eps + c/eps = 1 with both totally positive integers of norm < 1)."""
    from lattice import indecomposable_test
    for D in [2, 3, 5, 6, 7, 19, 31, 46]:
        is_indec, eta = indecomposable_test(D)
        assert is_indec(QF(D, 1, 0)) and is_indec(eta) and is_indec(eta * eta)
