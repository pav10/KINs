#!/usr/bin/env python3
"""Independent certificate for a LOWER bound s_sq(M) >= k (no use of the edge machinery for the
checks): M = O^2 + sum O G^{-1} y_g with Gram G.  Checks
  (1) G totally positive definite, M integral (B on a Z-basis of M lies in O_K);
  (2) for each listed vector v (integer coordinates in that Z-basis): Q(v) computed directly;
  (3) Q(v) indecomposable: lattice.indecomposable_test (CF orbit reps + Dress-Scharlau bound);
  (4) the values are pairwise in distinct classes mod (K^x)^2 (rqf.same_square_class).
The vectors are FOUND with the per-edge enumeration, then re-checked as above.
Usage: python3 scripts/certify_lattice.py D a b    (det-2 family frame [[a-sqrt D, b],[b, a+sqrt D]])"""
import os, sys, json
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from rqf import QF, same_square_class
from complete import FieldData, DiscModule, overlattice_basis, steinitz_ideal, is_principal
from lattice import (edge_deltas, enum_short_reduced, vec_value, Bform, indecomposable_test, hnf,
                     in_lattice)


def certify(D, G, gens, verbose=True):
    fd = FieldData(D)
    basis = overlattice_basis(fd, G, gens)
    # (1)
    assert G[0][0].is_tot_pos() and (G[0][0] * G[1][1] - G[0][1] * G[1][0]).is_tot_pos()
    for u in basis:
        for v in basis:
            assert Bform(G, u, v).is_integral()
    # find vectors (one per class) by the per-edge enumeration
    found = []
    for dlt in edge_deltas(D):
        T = [[int((dlt * Bform(G, basis[i], basis[j])).trace()) for j in range(4)] for i in range(4)]
        for _, vv in enum_short_reduced(T, 1):
            val = vec_value(D, G, basis, vv)
            if not any(same_square_class(val, w) for _, w in found):
                found.append((vv, val))
    # (2)-(4) independent re-checks
    is_ind, _ = indecomposable_test(D)
    for vv, val in found:
        x = sum((basis[i][0] * vv[i] for i in range(4)), QF(D, 0))
        y = sum((basis[i][1] * vv[i] for i in range(4)), QF(D, 0))
        q = Bform(G, (x, y), (x, y))
        assert q == val and is_ind(q), (vv, val)
    for i in range(len(found)):
        for j in range(i):
            assert not same_square_class(found[i][1], found[j][1])
    free = is_principal(fd, steinitz_ideal(fd, G, basis))
    if verbose:
        print(f"D={D}: CERTIFIED s_sq >= {len(found)} (M {'free' if free else 'NOT free'})")
        for vv, val in found:
            print(f"   v={vv}  Q(v)={val}  N={val.norm()}")
    return len(found), found, basis, free


if __name__ == "__main__":
    D, a, b = map(int, sys.argv[1:4])
    fd = FieldData(D)
    al, be, bb = QF(D, a, -1), QF(D, a, 1), QF(D, b)
    G = [[al, bb], [bb, be]]
    A = DiscModule(fd.ar, fd.ar.from_qf(al), fd.ar.from_qf(bb), fd.ar.from_qf(be))
    for gens in A.maximal_isotropic():
        print("glue generators (O^2/GO^2 coords):", gens)
        certify(D, G, gens)
