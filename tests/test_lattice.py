import random
import pytest
from fractions import Fraction as F
from rqf import QF
from gram import rank_over_K, largest_gram
from lattice import (hnf, in_lattice, det_hnf, enum_short, balanced_classes, find_winner,
                     realized_module, analyze_lattice, trace_gram)


def test_rank_pitfall_matrix():
    """(3/2)I - (1/2)J passes diag>0, 2x2 minors>=0, 3x3 principal minors=0 but has rank 4."""
    D = 2
    one, h = QF(D, 1, 0), QF(D, F(-1, 2), 0)
    M = [[one if i == j else h for j in range(4)] for i in range(4)]
    assert rank_over_K(M) == 4


def test_hnf_membership_random():
    rnd = random.Random(1)
    for _ in range(200):
        gens = [[rnd.randint(-6, 6) for _ in range(4)] for _ in range(rnd.randint(1, 6))]
        H = hnf(gens)
        for g in gens:
            assert in_lattice(H, g)
        # brute check: small combos
        coef = [rnd.randint(-2, 2) for _ in gens]
        v = [sum(c * g[j] for c, g in zip(coef, gens)) for j in range(4)]
        assert in_lattice(H, v)
        # a vector outside: perturb by a non-member if the lattice is not all of Z^4
        if len(H) < 4 or det_hnf(H) > 1:
            assert any(not in_lattice(H, e) for e in
                       ([1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]))
    assert det_hnf(hnf([[2, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1], [1, 0, 0, 0]])) == 1


def test_fincke_pohst_counts():
    I4 = [[int(i == j) for j in range(4)] for i in range(4)]
    assert len(enum_short(I4, 1)) == 4 and len(enum_short(I4, 2)) == 16
    D4 = [[2, -1, 0, 0], [-1, 2, -1, -1], [0, -1, 2, 0], [0, -1, 0, 2]]   # D4 root lattice
    assert len(enum_short(D4, 2)) == 12                                  # 24 roots / +-


@pytest.mark.slow
def test_max_gram_sizes():
    for D, e in [(19, 4), (22, 4), (31, 4)]:
        assert largest_gram(D)[0] == e


@pytest.mark.slow
def test_D31_linear_irreducibility_counterexample():
    """Q(sqrt31): winner binary lattice represents 6 indecomposable square classes
    (value traces 12,34,56,490,902,1314); only [6+sqrt31] has a Z-linearly
    irreducible representative in the trace lattice; shorter vectors span at 20."""
    D = 31
    M = find_winner(D, balanced_classes(D))
    Gp, vs, basis = realized_module(D, M)
    r = analyze_lattice(D, Gp, basis, 1400)
    assert sorted(int(c.trace()) for c, _ in r["classes"]) == [12, 34, 56, 490, 902, 1314]
    assert [c for c, _ in r["lin_irred_classes"]] == [QF(31, 6, 1)]
    assert r["nstar"] == 20


def test_D19_linear_irreducibility_gap():
    D = 19
    M = find_winner(D, balanced_classes(D))
    Gp, vs, basis = realized_module(D, M)
    r = analyze_lattice(D, Gp, basis, 600)
    assert len(r["classes"]) == 4 and len(r["lin_irred_classes"]) == 2


def test_exact_per_edge_method_D57_frame():
    """Delta=57 frame (imported, re-verified): delta = eta, s_sq = 6, 2 values on each of 6 edges."""
    from lattice import exact_sq_classes, edge_deltas
    from rqf import tp_unit
    D = 57
    w, one, z = QF(D, F(1, 2), F(1, 2)), QF(D, 1), QF(D, 0)
    a1, a2, b = one * 10 + w * 3, one * 23 + w * 7, -(one * 13 + w * 4)
    assert a1 * a2 - b * b == tp_unit(D)
    cl, pe = exact_sq_classes(D, [[a1, b], [b, a2]], [(one, z), (w, z), (z, one), (z, w)])
    assert len(cl) == 6 and pe == [2] * 6
    assert len(edge_deltas(19)) == 6 and len(edge_deltas(46)) == 12   # = period s


@pytest.mark.slow
def test_D43_eight_classes():
    """S(Q(sqrt43),2) >= 8: exact, cap-free (L6.2)."""
    from lattice import exact_sq_classes
    D = 43
    M = find_winner(D, balanced_classes(D))
    Gp, vs, basis = realized_module(D, M)
    cl, pe = exact_sq_classes(D, Gp, basis)
    assert len(cl) == 8 and sum(1 for x in pe if x) == 10
