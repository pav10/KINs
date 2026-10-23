"""R5 (docs/13_R5.md): the genus of <1, eps> for p prime, against R1's complete enumeration, and the
class count |A| = h(-8D)/2 = h(E)/(2 h(K)) under (H)."""
import os, sys, json
import pytest
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "scripts"))


@pytest.mark.parametrize("p", [7, 19, 31, 43, 59, 83])
def test_genus_sets_equal_r1(p):
    """Represented-class sets of the unimodular det-class-2 lattices found by R1 = the form model sets."""
    pytest.importorskip("cypari2")                 # h(-8p) via PARI (reported only)
    from r5_genus import compare
    r = compare(p)
    assert r["equal"] and r["distinct_r1_sets"] == r["distinct_model_sets"]
    assert 0 < r["pairs_hit"] <= r["A"]        # classes up to inversion with n(c) > 0


def test_genus_class_count():
    """#G = |A| = h(E)/h+(K) = h(-8D)/2 for all (H)-fields D <= 3000 (h(E) from data/r2_scan.jsonl)."""
    pari = pytest.importorskip("cypari2").Pari()
    n = 0
    for l in open(os.path.join(os.path.dirname(__file__), "..", "data", "r2_scan.jsonl")):
        r = json.loads(l)
        if r["D"] > 3000 or not r["two_eps_sq"]:
            continue
        hK = int(pari(f"quadclassunit({4 * r['D']})[1]"))
        assert r["h"] % (2 * hK) == 0 and r["h"] // (2 * hK) * 2 == int(pari(f"qfbclassno({-8 * r['D']}, 1)"))
        n += 1
    assert n > 200
