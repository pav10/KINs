"""R7 (docs/15_R7.md): the sharp per-face lemma for classical lattices (<= rd values per functional,
equality only for (M, T_delta) = Z^{rd}), S(K,2) <= 4s, the packing lemma, and the cubic min traces
(Thomas family: t(K) >= a when Z[rho] = O_K)."""
import os, sys, json, gzip, random, importlib
from math import isqrt, sqrt
import pytest
HERE = os.path.dirname(__file__)
sys.path.insert(0, os.path.join(HERE, "..", "scripts"))
from rqf import QF
from r6_squares import sqf

DATA = os.path.join(HERE, "..", "data")


def _records(name):
    op = gzip.open if name.endswith(".gz") else open
    with op(os.path.join(DATA, name), "rt") as f:
        return [json.loads(l) for l in f if l.strip()]


def test_D82_equality_case():
    """Theorem 15.1(b) equality: the R1 extremal lattice of Q(sqrt 82) (s = 1, S = 4) has (M, T_delta) = Z^4
    for its edge functional: det 1, exactly 4 pairs of norm-1 vectors, 4 distinct edge values."""
    pari = pytest.importorskip("cypari2").Pari()
    from complete import FieldData, overlattice_basis
    from lattice import edge_deltas, Bform
    D = 82
    fd = FieldData(D)
    G = [[QF(D, 10, -1), QF(D, -4, 0)], [QF(D, -4, 0), QF(D, 10, 1)]]
    basis = overlattice_basis(fd, G, [[1, 0, 1, 0]])
    (dl,) = edge_deltas(D)
    T = [[(dl * Bform(G, u, v)).trace() for v in basis] for u in basis]
    assert all(x.denominator == 1 for r in T for x in r)
    M = pari.matrix(4, 4, [int(x) for r in T for x in r])
    assert pari.matdet(M) == 1
    V = pari.qfminim(M, 1)[2]
    assert len(V) == 4
    vals = set()
    for k in range(4):
        c = [int(V[i, k]) for i in range(4)]
        x = tuple(sum((ci * u[j] for ci, u in zip(c, basis)), QF(D, 0)) for j in range(2))
        q = Bform(G, x, x)
        assert (dl * q).trace() == 1 and q.is_tot_pos()
        vals.add((q.a, q.b))
    assert len(vals) == 4


def test_edge_functional_norm():
    """Corollary 15.2: N(delta_E * different) = |N(mu_E)| for the intermediate convergent mu_E of every
    edge; it is 1 for some edge iff N(eps) = -1 (s odd)."""
    from indec import cf_list, convergents
    for D in [D for D in range(2, 300) if sqf(D)]:
        _, _, s = cf_list(D, 2)
        cf, dl, s = cf_list(D, 2 * s + 6)
        cvs = convergents(cf, 2 * s + 4)
        al = [QF(D, 1, 0)] + [QF(D, p, 0) + QF(D, q, 0) * dl for (p, q) in cvs]
        Delta = D if D % 4 == 1 else 4 * D
        ones = 0
        for i in range(0, 2 * s, 2):
            a0, a1 = al[i], al[i + 1]
            A = [[2 * a0.a, 2 * a0.b * D], [2 * a1.a, 2 * a1.b * D]]
            det = A[0][0] * A[1][1] - A[0][1] * A[1][0]
            d = QF(D, A[1][1] / det, -A[1][0] / det)          # Tr(d a0) = 1, Tr(d a1) = 0
            assert d.is_tot_pos() and a0.is_tot_pos() and a1.norm() < 0
            assert abs(d.norm()) * Delta == abs(a1.norm()), (D, i)
            ones += abs(a1.norm()) == 1
        assert (ones > 0) == (s % 2 == 1), D


def test_certified_data_respect_sharp_bound():
    """All certified records: S <= 4s (Cor. 15.2), <= 4 values per edge, <= 3 when s is even."""
    for name in ("complete_S.jsonl", "complete_S_long.jsonl"):
        for r in _records(name):
            assert r["S"] <= 4 * r["s"], r["D"]
            if r["best"]:                                     # None when S = 1 (no frame)
                assert max(r["best"]["per_edge"]) <= 4, r["D"]
            if r["s"] % 2 == 0:
                assert r["S"] <= 3 * r["s"]
    for r in _records("edge_stats.jsonl"):
        assert r["max_vals"] <= (4 if r["s"] % 2 else 3), r["D"]
    for name in ("r2_scan.jsonl", "r2_scan_100k.jsonl"):
        for r in _records(name):
            assert r["max_n"] <= 4 * r["s"]
    for r in _records("r3_scan_1e6.jsonl.gz"):
        assert r.get("nu_max", 0) <= 4 * r["s"]


def test_packing_lemma():
    """Lemma 15.4: #{v : T(v) <= t} <= (2 sqrt t + 1)^n for integral lattices (random Gram matrices A^T A)."""
    pari = pytest.importorskip("cypari2").Pari()
    random.seed(7)
    for _ in range(60):
        n = random.randint(1, 5)
        A = pari.matrix(n, n, [random.randint(-2, 2) for _ in range(n * n)])
        if pari.matdet(A) == 0:
            continue
        G = A.mattranspose() * A
        t = random.randint(1, 6)
        N = 1 + int(pari.qfminim(G, t)[0])
        assert N <= (2 * sqrt(t) + 1) ** n


def _r7():
    pytest.importorskip("cypari2")
    try:
        return importlib.import_module("r7_cubic")
    except ImportError:
        pytest.skip("sails not available (SAILS_PATH)")


def test_thomas_family_theorem():
    """Theorem 15.6: alpha = 1 + (a-1)(b - rho) is indecomposable in Z[rho] with m(alpha) = a (exact
    min trace and independent box test, on the order Z[rho]), for small a and b >= 2a-1; and t(O_K) = a
    on the maximal order when Z[rho] = O_K (sails)."""
    R = _r7()
    pari = R.pari
    for a, b in [(3, 5), (3, 6), (4, 7), (4, 9), (5, 9), (6, 11), (6, 12)]:
        pol = pari.Pol(list(reversed(R.thomas(a, b))))
        x = pari("x")
        Z = R.Order(pol, [x ** 0, x, x ** 2])               # the order Z[rho]
        rho = pari.Mod(x, pol)
        alpha = 1 + (a - 1) * (b - rho)
        m, delta = Z.min_trace(alpha)
        assert m == a, (a, b)
        assert R.box_indecomposable(Z, alpha)[0], (a, b)
        O = R.Order(pol, list(pari.nfbasis(pol)))
        assert O.min_trace(alpha)[0] >= a                    # (ii): every order containing Z[rho]
    r = R.field_data(R.thomas(4, 7))
    assert r["index"] == 1 and r["t"] == 4                   # disc 7249: refutes "t_3 <= 3"


def test_min_trace_quadratic_and_simplest_cubic():
    """Cross-checks of the min-trace code: t = 1 for real quadratic fields (t_2 = 1, IMPORTED Kala-Tinkova),
    t = 2 for simplest cubic fields (min traces 1 or 2, IMPORTED Kala-Tinkova)."""
    R = _r7()
    for D in (2, 3, 5, 6, 7, 10, 13, 19, 31, 43, 46, 58, 67, 82, 94):
        c = [-(D - 1) // 4, -1, 1] if D % 4 == 1 else [-D, 0, 1]
        assert R.field_data(c)["t"] == 1, D
    for a in range(-1, 8):
        assert R.field_data(R.simplest(a))["t"] == 2, a


def test_min_trace_independent_box():
    """The qfminim min traces agree with an independent LLL-box enumeration of alpha O^vee (r7_cubic.crosscheck):
    disc 564 (t = 3), 3137 (t = 4), 14229 (t = 5)."""
    R = _r7()
    for c in ([-2, -4, 2, 1], [1, -8, -3, 1], [1, -12, -6, 1]):
        assert max(m for _, m in R.crosscheck(c)) in (3, 4, 5)
