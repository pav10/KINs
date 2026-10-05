"""R1 (complete computation of S(K,2)): cross-checks of every step against independent code."""
import random
from fractions import Fraction as F
import pytest
from rqf import QF
from gram import valid_off_diagonal, omega
from lattice import (hnf, in_lattice, realized_module, exact_sq_classes, analyze_lattice,
                     balanced_classes, find_winner)
from complete import (FieldData, admissible_b, frames, DiscModule, overlattice_basis,
                      complete_S, steinitz_ideal, is_principal, is_principal_box)


def test_admissible_b_matches_box_search():
    for D in [2, 6, 13, 19, 33, 57]:
        fd = FieldData(D)
        R = fd.reps
        for i in range(len(R)):
            for j in range(i, len(R)):
                A = {(b.a, b.b) for b in admissible_b(fd, R[i], R[j])}
                B = {(b.a, b.b) for b in valid_off_diagonal(D, R[i], R[j])
                     if (R[i] * R[j] - b * b).is_tot_pos()}
                assert A == B, (D, i, j)


def test_unit_square_reps():
    """reps = indecomposables mod (O^x)^2: iota (N eps = -1) or 2 iota (N eps = +1);
    square classes = kappa_sq."""
    from indec import indecomposables, square_class_reps
    for D in [2, 3, 5, 6, 7, 13, 19, 21, 31, 33, 43, 57]:
        fd = FieldData(D)
        iota = len(indecomposables(D))
        assert len(fd.reps) == (2 * iota if fd.Neps > 0 else iota), D
        assert len(fd.class_reps) == len(square_class_reps(D)), D


# ---- independent brute-force enumeration of maximal integral overlattices ----------------
def _key(D, basis):
    rows = [[v[0].a, v[0].b, v[1].a, v[1].b] for v in basis]
    den = 1
    for r in rows:
        for c in r:
            den = den * F(c).denominator // __import__("math").gcd(den, F(c).denominator)
    return tuple(tuple(r) for r in hnf([[int(c * den) for c in r] for r in rows])), den


def _inv4(A):
    n = len(A)
    M = [row[:] + [F(int(i == j)) for j in range(n)] for i, row in enumerate(A)]
    for c in range(n):
        p = next(r for r in range(c, n) if M[r][c] != 0)
        M[c], M[p] = M[p], M[c]
        pv = M[c][c]
        M[c] = [x / pv for x in M[c]]
        for r in range(n):
            if r != c and M[r][c] != 0:
                f = M[r][c]
                M[r] = [M[r][j] - f * M[c][j] for j in range(2 * n)]
    return [row[n:] for row in M]


def _brute_maximal(D, G, return_all=False):
    """BFS over O-lattices L <= M <= L^# (L = O^2 with Gram G), extending by x + O x for x in
    a full set of coset reps of L^#/L, testing B-integrality on a Z-basis directly."""
    w, one, z = omega(D), QF(D, 1), QF(D, 0)
    det = G[0][0] * G[1][1] - G[0][1] * G[1][0]
    di = det.inv()
    Gi = [[G[1][1] * di, -G[0][1] * di], [-G[1][0] * di, G[0][0] * di]]
    Lb = [(one, z), (w, z), (z, one), (z, w)]
    Ld = [(Gi[0][0] * x + Gi[0][1] * y, Gi[1][0] * x + Gi[1][1] * y) for x, y in Lb]  # L^# basis

    def B(u, v):
        return G[0][0] * u[0] * v[0] + G[0][1] * (u[0] * v[1] + u[1] * v[0]) + G[1][1] * u[1] * v[1]

    def integral(basis):
        return all(B(u, v).is_integral() for u in basis for v in basis)

    def zbasis(vecs):
        rows = [[v[0].a, v[0].b, v[1].a, v[1].b] for v in vecs]
        den = 1
        for r in rows:
            for c in r:
                den = den * F(c).denominator // __import__("math").gcd(den, F(c).denominator)
        H = hnf([[int(c * den) for c in r] for r in rows])
        return [(QF(D, F(h[0], den), F(h[1], den)), QF(D, F(h[2], den), F(h[3], den))) for h in H]

    # coset reps of L^#/L, independently of complete.py: coordinates of L's basis in the
    # L^#-basis (rational solve), HNF box.
    Pd = [[F(c) for c in (v[0].a, v[0].b, v[1].a, v[1].b)] for v in Ld]
    inv = _inv4(Pd)
    C = []
    for l in Lb:
        q = [F(c) for c in (l[0].a, l[0].b, l[1].a, l[1].b)]
        c = [sum(q[k] * inv[k][j] for k in range(4)) for j in range(4)]
        assert all(t.denominator == 1 for t in c)
        C.append([int(t) for t in c])
    H = hnf(C)
    assert len(H) == 4
    reps = []
    for c in __import__("itertools").product(*[range(H[i][i]) for i in range(4)]):
        reps.append((sum((Ld[i][0] * c[i] for i in range(4)), QF(D, 0)),
                     sum((Ld[i][1] * c[i] for i in range(4)), QF(D, 0))))
    assert len(reps) == abs(det.norm())
    start = zbasis(Lb)
    todo, seen, maxi = [start], {_key(D, start)}, []
    while todo:
        M = todo.pop()
        ext = False
        for x in reps:
            M2 = zbasis(M + [x, (w * x[0], w * x[1])])
            k = _key(D, M2)
            if k == _key(D, M):
                continue
            if not integral(M2):
                continue
            ext = True
            if k not in seen:
                seen.add(k)
                todo.append(M2)
        if not ext:
            maxi.append(_key(D, M))
    if return_all:
        return seen, set(maxi)
    return set(maxi)


def test_maximal_overlattices_against_bruteforce():
    rnd = random.Random(7)
    checked = 0
    for D in [2, 3, 6, 13, 19, 33, 10, 57]:
        fd = FieldData(D)
        fr = list(frames(fd, galois=False))
        rnd.shuffle(fr)
        for key, al, be, b in fr:
            A = DiscModule(fd.ar, fd.ar.from_qf(al), fd.ar.from_qf(b), fd.ar.from_qf(be))
            if A.N > 40:
                continue
            G = [[al, b], [b, be]]
            mine = {_key(D, overlattice_basis(fd, G, gens)) for gens in A.maximal_isotropic()}
            assert mine == _brute_maximal(D, G), (D, key)
            checked += 1
            if checked % 6 == 0:
                break
    assert checked >= 20


def test_all_isotropic_submodules_against_bruteforce():
    rnd = random.Random(11)
    checked = 0
    for D in [6, 10, 15, 19, 26, 30, 33, 57]:
        fd = FieldData(D)
        fr = list(frames(fd, galois=False))
        rnd.shuffle(fr)
        for key, al, be, b in fr[:12]:
            A = DiscModule(fd.ar, fd.ar.from_qf(al), fd.ar.from_qf(b), fd.ar.from_qf(be))
            if A.N > 40:
                continue
            G = [[al, b], [b, be]]
            subs = A.isotropic_submodules()
            mine = {_key(D, overlattice_basis(fd, G, g)) for _, g, _ in subs}
            assert len(mine) == len(subs)
            allb, _ = _brute_maximal(D, G, return_all=True)
            assert mine == allb, (D, key)
            checked += 1
    assert checked >= 20


def test_S_free_against_bruteforce():
    """S_free (free lattices only, h(K) = 2): max of s_sq over ALL free integral overlattices of
    all frames (brute force) == complete_S's free-maximal evaluation."""
    for D in [10, 15, 30, 26]:
        fd = FieldData(D)
        best = 1
        for key, al, be, b in frames(fd):
            A = DiscModule(fd.ar, fd.ar.from_qf(al), fd.ar.from_qf(b), fd.ar.from_qf(be))
            G = [[al, b], [b, be]]
            for _, g, _ in A.isotropic_submodules():
                basis = overlattice_basis(fd, G, g)
                if is_principal(fd, steinitz_ideal(fd, G, basis)):
                    best = max(best, len(exact_sq_classes(D, G, basis)[0]))
        assert complete_S(D)["S_free"] == best, D


def test_galois_reduction_harmless():
    for D in [13, 19, 33]:
        assert complete_S(D, galois=False)["S"] == complete_S(D)["S"]


@pytest.mark.slow
def test_exact_vs_capped_enumeration():
    """s_sq of a few maximal overlattices: exact (per-edge) == capped Fincke-Pohst classes."""
    for D, cap in [(19, 2000), (33, 1500)]:
        fd = FieldData(D)
        n = 0
        for key, al, be, b in frames(fd):
            A = DiscModule(fd.ar, fd.ar.from_qf(al), fd.ar.from_qf(b), fd.ar.from_qf(be))
            G = [[al, b], [b, be]]
            for gens in A.maximal_isotropic():
                basis = overlattice_basis(fd, G, gens)
                cl, _ = exact_sq_classes(D, G, basis)
                r = analyze_lattice(D, G, basis, cap)
                assert len(r["classes"]) <= len(cl)
                if len(cl) >= 4:
                    assert len(r["classes"]) == len(cl), (D, key)
                    n += 1
            if n >= 3:
                break


def test_realized_module_is_OK_module():
    """C17 regression: realized_module must give the O_K-span (it gave Z[sqrt D]-span for
    D = 1 mod 4).  <1, (5+sqrt13)/2> over Q(sqrt13) represents 2 indecomposable classes."""
    D = 13
    one, z = QF(D, 1), QF(D, 0)
    M = [[one, z], [z, QF(D, F(5, 2), F(1, 2))]]
    Gp, vs, basis = realized_module(D, M)
    w = omega(D)
    rows = [[v[0].a, v[0].b, v[1].a, v[1].b] for v in basis]
    H = hnf([[int(4 * c) for c in r] for r in rows])
    for v in basis:
        wv = (w * v[0], w * v[1])
        assert in_lattice(H, [int(4 * c) for c in (wv[0].a, wv[0].b, wv[1].a, wv[1].b)])
    cl, _ = exact_sq_classes(D, Gp, basis)
    assert len(cl) == 2


def _random_ideals(fd, rnd, k):
    w, ar = omega(fd.D), fd.ar
    out = []
    while len(out) < k:
        g = [ar.to_qf((rnd.randint(-12, 12), rnd.randint(-6, 6))) for _ in range(2)]
        H = hnf([list(ar.from_qf(x)) for x in (g[0], g[0] * w, g[1], g[1] * w)])
        if len(H) == 2:
            out.append((H, [ar.to_qf(tuple(H[0])), ar.to_qf(tuple(H[1]))]))
    return out


def test_free_detection():
    """D = 10 (h = 2): the ideal (2, sqrt10) is not principal, (1) is; CF (Serret) test ==
    norm-box search on random ideals."""
    fd = FieldData(10)
    assert is_principal(fd, [QF(10, 1), QF(10, 0, 1)])
    assert not is_principal(fd, [QF(10, 2), QF(10, 0, 1)])
    assert is_principal(fd, [QF(10, 3), QF(10, 0, 3)])
    rnd = random.Random(5)
    for D in [10, 15, 26, 30, 65, 79, 82, 85]:
        fd = FieldData(D)
        for _, I in _random_ideals(fd, rnd, 15):
            assert is_principal(fd, I) == is_principal_box(fd, I), (D, I)


def test_free_detection_against_pari():
    pari = pytest.importorskip("cypari2").Pari()
    isp = pari("(b,i)->bnfisprincipal(b,i,0)==0")
    rnd = random.Random(3)
    for D in [10, 15, 26, 30, 51, 58, 79, 82, 94, 229, 1155, 4171]:
        fd = FieldData(D)
        pol = "x^2-x-%d" % ((D - 1) // 4) if D % 4 == 1 else "x^2-%d" % D
        bnf = pari("bnfinit(%s,1)" % pol)
        for H, I in _random_ideals(fd, rnd, 20):
            M = pari("[%d,%d;%d,%d]" % (H[0][0], H[1][0], H[0][1], H[1][1]))
            assert is_principal(fd, I) == bool(isp(bnf, pari.idealhnf(bnf, M))), (D, H)


def test_complete_small_fields():
    """Exact S(K,2); D = 26 -> 2, 33 -> 4, 19 -> 6 reproduce the imported complete values."""
    expect = {2: 2, 3: 2, 5: 1, 6: 2, 13: 2, 21: 2, 26: 2, 33: 4, 19: 6, 57: 6}
    for D, e in expect.items():
        r = complete_S(D)
        assert r["S"] == e, (D, r["S"])
        assert r["S_free"] == e, (D, r["S_free"])


def test_sails_cross_check():
    """Indecomposables, facet functionals (= edge_deltas), d = 1 and #facets against the
    independent PARI-based `sails` package (skipped if sails/cypari2 are not installed).
    Full range: SAILS_PATH=../sails python3 scripts/check_sails.py 2 3000 (all pass)."""
    import os, sys, importlib
    sp = os.environ.get("SAILS_PATH", os.path.join(os.path.dirname(__file__), "..", "..", "sails"))
    sys.path.insert(0, sp)
    pytest.importorskip("cypari2")
    pytest.importorskip("sails.compute")
    sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "scripts"))
    check = importlib.import_module("check_sails").check
    for D in [2, 3, 5, 13, 19, 31, 43, 46, 57, 67, 94, 181, 331]:
        check(D)
