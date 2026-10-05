#!/usr/bin/env python3
"""Lower bounds for S(K,2) by D: best-found lattice (winner on balanced and on orbit reps),
s_sq computed EXACTLY by the per-edge minimal-vector method (lattice.exact_sq_classes).
Usage: python3 scripts/scan_exact.py Dmin Dmax timeout_per_search out.jsonl"""
import os, sys, json, signal, time
from math import isqrt
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from rqf import cf_quadratic, fundamental_unit
from indec import indecomposables, square_class_reps
from lattice import balanced_classes, find_winner, realized_module, exact_sq_classes, trace_gram

class TO(Exception): pass
def _h(*a): raise TO()
signal.signal(signal.SIGALRM, _h)

def sqf(D):
    r = isqrt(D); return r*r != D and all(D % (p*p) for p in range(2, r+1))

Dmin, Dmax, tl, out = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), sys.argv[4]
for D in range(Dmin, Dmax+1):
    if not sqf(D): continue
    P, Q = (1, 2) if D % 4 == 1 else (0, 1)
    pre, per = cf_quadratic(P, Q, D)
    rec = dict(D=D, s=len(per), iota=len(indecomposables(D)), kappa=len(square_class_reps(D)),
               Neps=fundamental_unit(D)[1])
    for kind in ("balanced", "orbit"):
        R = balanced_classes(D) if kind == "balanced" else indecomposables(D)
        signal.alarm(tl)
        try:
            t = time.time()
            M = find_winner(D, R)
            if len(M) < 2:                      # distinct classes are never collinear
                signal.alarm(0)
                rec[kind] = dict(winner=len(M), s_sq=len(M), edges_touched=None, per_edge=None,
                                 gram=None, diag=[str(M[0][0])], sec=0)
                continue
            Gp, vs, basis = realized_module(D, M)
            cl, pe = exact_sq_classes(D, Gp, basis)
            signal.alarm(0)
            rec[kind] = dict(winner=len(M), s_sq=len(cl), edges_touched=sum(1 for x in pe if x),
                             per_edge=pe, gram=[[str(x) for x in r] for r in Gp],
                             diag=[str(M[i][i]) for i in range(len(M))], sec=round(time.time()-t, 1))
        except TO:
            rec[kind] = "timeout"
        finally:
            signal.alarm(0)
    with open(out, "a") as f:
        f.write(json.dumps(rec) + "\n")
    b = rec["balanced"]; o = rec["orbit"]
    print(D, rec["s"], rec["kappa"], b if isinstance(b, str) else b["s_sq"],
          o if isinstance(o, str) else o["s_sq"], flush=True)
