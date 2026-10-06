#!/usr/bin/env python3
"""Independent cross-check of the real quadratic input of L1.4/L1.5/L6.1 and R1 against the
`sails` package (github.com/pav10/sails: PARI-based exact sails of totally real fields).

For every squarefree D in the range, with the defining polynomial of omega (so sails' power basis
is our O-basis 1, omega), checks
  (a) indecomposables mod O^{x,+}: sails == indec.indecomposables(D)            (Dress-Scharlau CF)
  (b) every facet of the totally positive sail has integer distance d = 1     (Kala-Tinkova 3.1)
  (c) facet functionals mod O^{x,+}: sails == lattice.edge_deltas(D)           (our edge functionals)
  (d) the facet lattice points (mod O^{x,+}) are exactly the indecomposables
  (e) #facet orbits = s/2 (s even) or s (s odd)
Usage: SAILS_PATH=/path/to/sails python3 scripts/check_sails.py Dmin Dmax
"""
import os, sys
from fractions import Fraction as F
from math import isqrt
HERE = os.path.dirname(__file__)
sys.path.insert(0, os.path.join(HERE, "..", "src"))
sys.path.insert(0, os.environ.get("SAILS_PATH", os.path.join(HERE, "..", "..", "sails")))
from sails.compute import compute_field
from rqf import QF, tp_unit, cf_quadratic
from indec import indecomposables, normalize_unit
from lattice import edge_deltas


def sqf(D):
    r = isqrt(D)
    return r * r != D and all(D % (p * p) for p in range(2, r + 1))


def check(D):
    if D % 4 == 1:
        coeffs, to_qf = [-(D - 1) // 4, -1, 1], (lambda c: QF(D, c[0] + F(c[1], 2), F(c[1], 2)))
    else:
        coeffs, to_qf = [-D, 0, 1], (lambda c: QF(D, c[0], c[1]))
    rec = compute_field(coeffs, label=f"Q{D}", tp_only=True).record()
    orth = rec["orthants"][0]
    assert orth["signature"] == [0, 0]
    eta = tp_unit(D)
    key = lambda x: (lambda y: (y.a, y.b))(normalize_unit(x, eta))
    ours = {key(a) for a in indecomposables(D)}
    theirs = {key(to_qf(c)) for _, c in orth["indecomposables"]}
    assert ours == theirs, (D, "indecomposables", len(ours), len(theirs))              # (a)
    deltas = set()
    pts = set()
    for f in orth["facets"]:
        assert f["d"] == 1, (D, "integer distance", f["d"])                            # (b)
        c0, c1 = f["c"]                                     # Tr(delta) = c0, Tr(delta*omega) = c1
        x = F(c0, 2)
        y = F(c1, 2 * D) if D % 4 != 1 else (F(c1) - x) / D
        dl = QF(D, x, y)
        assert dl.is_tot_pos()
        deltas.add(key(dl))
        for p in f["lattice_points"]:
            q = to_qf(p)
            assert (dl * q).trace() == 1
            pts.add(key(q))
    assert deltas == {key(d) for d in edge_deltas(D)}, (D, "facet functionals")      # (c)
    assert pts == ours, (D, "facet points")                                          # (d)
    P, Q = (1, 2) if D % 4 == 1 else (0, 1)
    s = len(cf_quadratic(P, Q, D)[1])
    assert len(orth["facets"]) == (s // 2 if s % 2 == 0 else s), (D, "facets", s)    # (e)
    return len(ours), len(orth["facets"]), s


if __name__ == "__main__":
    lo, hi = int(sys.argv[1]), int(sys.argv[2])
    n = 0
    for D in range(lo, hi + 1):
        if sqf(D):
            iota, f, s = check(D)
            n += 1
    print(f"sails cross-check OK for {n} squarefree D in [{lo}, {hi}]")
