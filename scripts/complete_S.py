#!/usr/bin/env python3
"""Roadmap R1: exact S(K,2) (all integral binary lattices) and S_free(K,2) (free lattices),
K = Q(sqrt D), by the provably complete frame + maximal-overlattice search (src/complete.py,
docs/09_R1_complete.md).
Usage: python3 scripts/complete_S.py out.jsonl D1 [D2 ...]   or   ... out.jsonl --range Dmin Dmax"""
import os, sys, json
from math import isqrt
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from complete import complete_S
from rqf import cf_quadratic


def sqf(D):
    r = isqrt(D)
    return r * r != D and all(D % (p * p) for p in range(2, r + 1))


out = sys.argv[1]
if sys.argv[2] == "--range":
    Ds = [D for D in range(int(sys.argv[3]), int(sys.argv[4]) + 1) if sqf(D)]
else:
    Ds = [int(x) for x in sys.argv[2:]]
for D in Ds:
    P, Q = (1, 2) if D % 4 == 1 else (0, 1)
    s = len(cf_quadratic(P, Q, D)[1])
    r = complete_S(D, verbose=True)
    r["s"] = s
    with open(out, "a") as f:
        f.write(json.dumps(r) + "\n")
    print(f"D={D} s={s} kappa={r['kappa']} S={r['S']} S_free={r['S_free']} frames={r['frames']} "
          f"lattices={r['lattices']} hist={r['hist']} {r['sec']}s", flush=True)
