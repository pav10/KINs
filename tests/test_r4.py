"""R4 (docs/12_R4.md): the imported claims -- Thm 22 table re-run, Prop 5 family, regulator inequality,
incomparability of indecomposables (Thm 4)."""
import os, sys, json
import pytest
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "scripts"))
from rqf import QF


def test_thm22_table():
    """Non-trivial unit-determinant frames reproduce the imported table except D = 22 (4, not 2-3);
    the 'none' fields have only the trivial frame <1, eps> (C26)."""
    from r4_thm22 import unit_frames_S, IMPORTED
    for D, imp in IMPORTED.items():
        r = unit_frames_S(D)
        ntr = r["unit_frames_sq"] + r["unit_frames_eps"]
        val = max(r["S_unit_sq"] or 0, r["S_unit_eps"] or 0)
        if imp == "none":
            assert ntr == 0, D
        elif imp == "2-3":
            assert ntr > 0 and val == (4 if D == 22 else 2), (D, val)
        else:
            assert val == int(imp.lstrip(">=")), (D, val)


def test_prop5_family():
    """Prop 5: integral, positive definite, two indecomposables in distinct classes, 1 - cos <= 2/(n-2)."""
    from r4_checks import prop5, sqf
    for n in (6, 8, 10, 12, 20, 30, 60, 100):
        if not sqf(n * n + 2):
            continue
        ok, c, c2 = prop5(n)
        assert all(ok.values()), (n, ok)
        assert 0 < 1 - c <= 2 / (n - 2) and c2 < 0


def test_regulator_inequality():
    """log eps >= (s - 2) log phi + log(sqrt D - 1), hence s_sq <= 12 s <= 24 + 24.94 (Reg - log(sqrt D - 1))."""
    from r4_checks import reg_check
    assert reg_check(3000) == []


def test_indecomposables_incomparable():
    """Thm 4's input: distinct indecomposables are never comparable (difference not totally positive)."""
    from indec import indecomposables
    from rqf import tp_unit
    for D in (19, 43, 57, 94):
        eta = tp_unit(D)
        I = []
        for g in indecomposables(D):
            I += [g, g * eta, g * eta.conj()]
        for a in I:
            for b in I:
                if a != b:
                    assert not (a - b).is_tot_pos()
