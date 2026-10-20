"""R3 (docs/11_R3.md): the elementary model of the R2 class distribution -- CF norms, exact square
classes, the x^2 + 2y^2 criterion, binary forms of discriminant -8p -- against the R2/PARI code."""
import os, sys, glob, json, random
from fractions import Fraction as F
import pytest
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "scripts"))
from rqf import QF
from r3_tools import (cf_sqrt, has_H, indec_data, convergents_big, element, Sieve, field,
                      nu_distribution, qf_compose, qf_reduce, h_neg, zagier_T, rep_x2_2y2)


def _sqf(D):
    return all(D % (p * p) for p in range(2, int(D ** .5) + 1))


def test_cf_norms_match_indecomposables():
    """beta_{j,t} = alpha_{j-1} + t alpha_j has norm Q_j + 2t P_{j+1} - t^2 Q_{j+1} = (D - u^2)/Q_{j+1},
    and these are exactly the indecomposables mod <eps> (indec.indecomposables)."""
    from indec import indecomposables
    for D in [D for D in range(3, 800, 4) if _sqf(D)]:
        cf = cf_sqrt(D)
        a, P, Q = cf
        ind = indec_data(D, cf)
        assert sorted(N for _, _, N in ind) == sorted(int(x.norm()) for x in indecomposables(D)), D
        conv = convergents_big(D, cf)
        for j, t, N in ind:
            x, y = element(D, conv, j, t)
            u = Q[j + 1] * t - P[j + 1]
            assert x * x - D * y * y == N and Q[j + 1] * N == D - u * u


def test_H_criterion_and_kappa():
    """(H) <=> some Q_k = 2; exact square-class count kappa_sq = FieldData's."""
    from complete import FieldData
    from r2_scan import two_eps_square
    sv = Sieve(4 * 700)
    for D in [D for D in range(3, 700, 4) if _sqf(D)]:
        H = has_H(D)
        assert H == two_eps_square(D), D
        if H:
            assert field(D, sv)["kappa"] == len(FieldData(D).class_reps), D


def test_forms():
    """Composition is associative/commutative on reduced forms; class numbers by counting."""
    random.seed(0)
    for D in (43, 163, 691):
        disc = -8 * D
        forms = [qf_reduce((a, b, (b * b - disc) // (4 * a))) for a in range(1, 40)
                 for b in range(-a + 1, a + 1, 1) if b % 2 == 0 and (b * b - disc) % (4 * a) == 0]
        forms = sorted(set(forms))
        for _ in range(200):
            f, g, k = (random.choice(forms) for _ in range(3))
            assert qf_compose(f, g) == qf_compose(g, f)
            assert qf_compose(qf_compose(f, g), k) == qf_compose(f, qf_compose(g, k))
    assert [h_neg(8 * D) for D in (3, 7, 43, 691, 5419)] == [2, 4, 10, 38, 82]


def test_nu_equals_pari_class_distribution():
    """Theorem 11.3: for p prime, {n(c)} over Cl(K(sqrt-2)) (PARI bnf) = {nu(f)} over Cl(-8p)^2."""
    pytest.importorskip("cypari2")
    from complete import FieldData
    from r2_classes import class_distribution
    sv = Sieve(4 * 5419)
    for p in (43, 163, 691, 1879, 5419):
        fd = FieldData(p)
        dist, _, _, _ = class_distribution(p, 2, fd=fd)
        nu = nu_distribution(p, sv=sv)
        assert sorted(len(v) for v in dist.values()) == sorted(nu.values()), p


def test_N_E_criterion_against_pari():
    """Lemma 11.2: gamma O_K is a norm of an integral O_E-ideal <=> N(gamma) = x^2 + 2y^2."""
    pytest.importorskip("cypari2")
    from complete import FieldData
    from r2_classes import EData
    for D in (19, 43, 67, 139, 179, 487):
        fd, ed = FieldData(D), EData(D, 2)
        for g in fd.reps:
            n = int(g.norm())
            fac = Sieve(n + 1).factor(n) if n > 1 else {}
            assert bool(ed.half_norm_ideals(ed.elt(g))) == rep_x2_2y2(fac), (D, g)


def test_zagier_hirzebruch():
    """h(Q(sqrt p)) = 1  =>  sum_{i<=s} a_i = T(p)/2 (Zagier) and sum_even - sum_odd = 3 h(-p)
    (Hirzebruch), p = 3 mod 4 prime (IMPORTED identities, checked here)."""
    pari = pytest.importorskip("cypari2").Pari()
    sv = Sieve(4 * 2000)
    n = 0
    for p in range(7, 2000, 4):
        if not sv.is_prime(p) or int(pari(f"quadclassunit({4 * p})[1]")) != 1:
            continue
        a, _, _ = cf_sqrt(p)
        s = len(a) - 1
        odd, even = sum(a[1:s:2]), sum(a[2:s + 1:2])
        assert odd + even == zagier_T(p, sv) // 2
        assert even - odd == 3 * int(pari(f"qfbclassno({-p}, 1)"))
        n += 1
    assert n > 100


@pytest.mark.parametrize("path", sorted(glob.glob(os.path.join(os.path.dirname(__file__), "..", "data",
                                                               "certificates", "D*_form.json"))))
def test_stored_form_certificates(path):
    """data/certificates/D*_form.json (ideal lattices (b O_E, N/a) from a form class of disc -8p), from
    the file alone: integrality; Q(v) = value; indecomposability by the edge-functional witness delta
    (delta >> 0, Tr(delta O) in Z, Tr(delta gamma) = 1); pairwise distinct square classes; and, for
    small p, the CF orbit test as well."""
    from lattice import Bform, indecomposable_test
    from rqf import same_square_class
    c = json.load(open(path))
    D = c["D"]
    q = lambda p: QF(D, F(p[0]), F(p[1]))
    G = [[q(x) for x in row] for row in c["gram"]]
    basis = [(q(v[0]), q(v[1])) for v in c["zbasis"]]
    assert G[0][0].is_tot_pos() and (G[0][0] * G[1][1] - G[0][1] * G[1][0]).is_tot_pos()
    assert all(Bform(G, u, v).is_integral() for u in basis for v in basis)
    is_ind = indecomposable_test(D)[0] if D < 20000 else None
    vals = []
    for vv, val, dl in zip(c["vectors"], c["values"], c["deltas"]):
        x = sum((basis[i][0] * vv[i] for i in range(4)), QF(D, 0))
        y = sum((basis[i][1] * vv[i] for i in range(4)), QF(D, 0))
        Qv, d = Bform(G, (x, y), (x, y)), q(dl)
        assert Qv == q(val) and Qv.is_tot_pos()
        assert d.is_tot_pos() and (d * Qv).trace() == 1
        assert d.trace().denominator == 1 and (d * QF(D, 0, 1)).trace().denominator == 1
        assert is_ind is None or is_ind(Qv)
        assert not any(same_square_class(Qv, w) for w in vals)
        vals.append(Qv)
    assert len(vals) == c["s_sq_lower_bound"] == c["nu_f"]
