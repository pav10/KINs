#!/usr/bin/env python3
"""R3(b) probe: the Galois-symmetric det-2 frames  G_{a,b} = [[a - sqrt D, b], [b, a + sqrt D]],
D = a^2 - b^2 - 2 squarefree (det G = N(a - sqrt D) - b^2 = 2).  These give the extremal lattices
for D = 19, 22, 38, 43, 46, 58, 73, 82, 862 (ledger L7.8).  For each (a, b) with b^2 + 2 <= Delta/4
(Dress-Scharlau: needed for a - sqrt D to be indecomposable), evaluates EXACTLY s_sq of every
maximal integral overlattice of O^2 with Gram G_{a,b} (CERTIFIED lower bounds for S(K,2)).
Usage: python3 scripts/family_det2.py out.jsonl Dmax [Dmin]"""
import os, sys, json, time
from math import isqrt
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from rqf import QF
from complete import FieldData, DiscModule, evaluate, steinitz_ideal, is_principal, overlattice_basis
from lattice import indecomposable_test


def sqf(D):
    r = isqrt(D)
    return r * r != D and all(D % (p * p) for p in range(2, r + 1))


out, Dmax = sys.argv[1], int(sys.argv[2])
Dmin = int(sys.argv[3]) if len(sys.argv) > 3 else 2
todo = []
for a in range(2, isqrt(Dmax + 2) + 2):
    for b in range(0, a):
        D = a * a - b * b - 2
        if Dmin <= D <= Dmax and D > 1 and sqf(D):
            dl4 = D if D % 4 != 1 else D / 4
            if b * b + 2 <= dl4:
                todo.append((D, a, b))
todo.sort()
for D, a, b in todo:
    t = time.time()
    fd = FieldData(D)
    al, be, bb = QF(D, a, -1), QF(D, a, 1), QF(D, b)
    is_ind, _ = indecomposable_test(D)
    G = [[al, bb], [bb, be]]
    A = DiscModule(fd.ar, fd.ar.from_qf(al), fd.ar.from_qf(bb), fd.ar.from_qf(be))
    res = []
    for gens in A.maximal_isotropic():
        n, cl, pe, basis = evaluate(fd, G, gens)
        free = is_principal(fd, steinitz_ideal(fd, G, basis))
        res.append(dict(s_sq=n, per_edge=pe, free=free, gens=[list(g) for g in gens]))
    best = max(r["s_sq"] for r in res)
    rec = dict(D=D, a=a, b=b, s=len(pe), alpha_indec=is_ind(al), kappa=len(fd.class_reps),
               s_sq=best, lattices=res, sec=round(time.time() - t, 1))
    with open(out, "a") as f:
        f.write(json.dumps(rec) + "\n")
    print(f"D={D} a={a} b={b} s={rec['s']} indec={rec['alpha_indec']} kappa={rec['kappa']} "
          f"s_sq={[r['s_sq'] for r in res]} free={[r['free'] for r in res]} {rec['sec']}s", flush=True)
