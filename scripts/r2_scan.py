#!/usr/bin/env python3
"""R2 scan: for squarefree D = 3 (mod 4), E = K(sqrt -2): the class distribution n(c) of the
indecomposables over Cl(O_E) (scripts/r2_classes.py), its maximum over c, and -- for the det-2
family lattices of data/family_det2.jsonl -- the model value n(c_M) next to the EXACT s_sq(M).
max_c n(c) is the exact maximum of s_sq over binary lattices with det in 2 (K^x)^2 when 2*eps is a
square in K (docs/10_R2.md, Theorem 10.2); otherwise it is an upper bound.
Usage: python3 scripts/r2_scan.py out.jsonl Dmin Dmax [i/m]"""
import os, sys, json, time
from math import isqrt
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "src"))
sys.path.insert(0, HERE)
from rqf import QF, fundamental_unit
from complete import FieldData, overlattice_basis
from r2_classes import class_distribution, lattice_class


def sqf(D):
    r = isqrt(D)
    return r * r != D and all(D % (p * p) for p in range(2, r + 1))


def two_eps_square(D):
    e, N = fundamental_unit(D)
    return N == 1 and any((QF(D, 2) * x).is_square_in_K() for x in (e, e.conj()))


if __name__ == "__main__":
    out, lo, hi = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
    shard = tuple(map(int, sys.argv[4].split("/"))) if len(sys.argv) > 4 else (0, 1)
    done = set()
    if os.path.exists(out):
        done = {json.loads(l)["D"] for l in open(out)}
    fam = {}
    fp = os.path.join(HERE, "..", "data", "family_det2.jsonl")
    for l in open(fp):
        r = json.loads(l)
        fam.setdefault(r["D"], []).append(r)
    Ds = [D for D in range(hi, lo - 1, -1) if D % 4 == 3 and sqf(D)]
    Ds = [D for k, D in enumerate(Ds) if k % shard[1] == shard[0] and D not in done]
    for D in Ds:
        t = time.time()
        fd = FieldData(D)
        dist, per_rep, ed, _ = class_distribution(D, 2, fd=fd)
        h = 1
        for c in ed.cyc:
            h *= c
        counts = sorted(((len(v), c) for c, v in dist.items()), reverse=True)
        rec = dict(D=D, two_eps_sq=two_eps_square(D), cyc=ed.cyc, h=h, kappa=len(fd.class_reps),
                   s=None, max_n=counts[0][0] if counts else 0,
                   best_class=list(counts[0][1]) if counts else None,
                   top=[n for n, _ in counts[:8]],
                   avg=round(sum(len(v) for v in dist.values()) / h, 3), fam=[])
        from rqf import cf_quadratic
        rec["s"] = len(cf_quadratic(0, 1, D)[1])
        for r in fam.get(D, []):
            al, be, bb = QF(D, r["a"], -1), QF(D, r["a"], 1), QF(D, r["b"])
            G = [[al, bb], [bb, be]]
            for lat in r["lattices"]:
                basis = overlattice_basis(fd, G, [tuple(g) for g in lat["gens"]])
                cM, _ = lattice_class(ed, G, basis)
                rec["fam"].append(dict(a=r["a"], b=r["b"], exact=lat["s_sq"],
                                       n_cM=len(dist.get(cM, ()))))
        rec["sec"] = round(time.time() - t, 1)
        with open(out, "a") as f:
            f.write(json.dumps(rec) + "\n")
        print(D, rec["two_eps_sq"], rec["cyc"], rec["kappa"], rec["max_n"],
              [(x["exact"], x["n_cM"]) for x in rec["fam"]], rec["sec"], flush=True)
