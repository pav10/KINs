#!/usr/bin/env python3
"""Cheap prefilter for the R2 search: for squarefree D = 3 (mod 4) in [Dmin, Dmax] with 2*eps a
square in K, the proxy  r = 2*iota / h(-8D)  (2*iota ~ kappa_sq: indecomposables mod (O^x)^2 via the
Blomer-Kala formula; h(-8D) = class number of Q(sqrt -2D), ~ h(E) when h(K) = 1).  The full class
distribution (scripts/r2_scan.py) is then run on the fields with the largest r.
Usage: python3 scripts/r2_prefilter.py out.jsonl Dmin Dmax [i/m]"""
import os, sys, json
from math import isqrt
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "src"))
from rqf import QF, cf_quadratic, fundamental_unit
import cypari2

pari = cypari2.Pari()


def sqf(D):
    r = isqrt(D)
    return r * r != D and all(D % (p * p) for p in range(2, r + 1))


out, lo, hi = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
shard = tuple(map(int, sys.argv[4].split("/"))) if len(sys.argv) > 4 else (0, 1)
Ds = [D for D in range(lo, hi + 1) if D % 4 == 3 and sqf(D)]
Ds = [D for k, D in enumerate(Ds) if k % shard[1] == shard[0]]
with open(out, "a") as f:
    for D in Ds:
        _, per = cf_quadratic(0, 1, D)
        s = len(per)
        iota = sum(per[i] for i in range(0, s - 1, 2)) if s % 2 == 0 else sum(per)
        e, N = fundamental_unit(D)
        sq = N == 1 and any((QF(D, 2) * x).is_square_in_K() for x in (e, e.conj()))
        if not sq:
            continue
        hK = int(pari(f"quadclassunit({4 * D})[1]"))
        hm = int(pari(f"qfbclassno({-8 * D})"))
        f.write(json.dumps(dict(D=D, s=s, iota=iota, hK=hK, h8=hm, r=2 * iota / hm)) + "\n")
