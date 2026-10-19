#!/usr/bin/env python3
"""Build the ideal lattice M_A = (A, N_{E/K}/n), E = K(sqrt -2), for the class c maximizing n(c)
(scripts/r2_classes.py), and certify s_sq(M_A) >= k exactly (scripts/certify_lattice.py), with
independent re-checks by sails (indecomposables) and PARI nfroots (square classes) when available.
Writes data/certificates/D{D}_ideal.json.   Usage: python3 scripts/certify_ideal.py D"""
import os, sys, json
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "src"))
sys.path.insert(0, HERE)
from rqf import QF, tp_unit
from indec import normalize_unit
from r2_classes import class_distribution, lattice_from_ideal, pari
from certify_lattice import certify_basis


def best_ideal(D, d=2):
    dist, per_rep, ed, fd = class_distribution(D, d)
    c = max(dist, key=lambda c: len(dist[c]))
    for i, g in enumerate(fd.reps):
        for B in ed.half_norm_ideals(ed.elt(g)):
            if ed.cls(B) == c:
                return ed, fd, B, g, c, len(dist[c])
    raise RuntimeError


def main(D):
    ed, fd, B, g, c, n_pred = best_ideal(D)
    G, basis = lattice_from_ideal(ed, B, g)
    k, found, basis, free = certify_basis(D, G, basis, verbose=False)
    vals = [v for _, v in found]
    extra = {}
    issq = pari("(nf,a)->#nfroots(nf, x^2 - a) > 0")
    nf = pari(f"nfinit(y^2-{D})")
    ps = lambda q: pari(f"{q.a} + ({q.b})*y")
    assert not any(issq(nf, ps(vals[i] * vals[j])) for i in range(k) for j in range(i))
    extra["pari_square_classes"] = True
    try:
        sys.path.insert(0, os.environ.get("SAILS_PATH", os.path.join(HERE, "..", "..", "sails")))
        from sails.compute import compute_field
        rec = compute_field([-D, 0, 1] if D % 4 != 1 else [-(D - 1) // 4, -1, 1], tp_only=True).record()
        eta = tp_unit(D)
        key = lambda q: (lambda y: (y.a, y.b))(normalize_unit(q, eta))
        if D % 4 != 1:
            si = {key(QF(D, cc[0], cc[1])) for _, cc in rec["orthants"][0]["indecomposables"]}
            assert all(key(v) in si for v in vals)
            extra["sails_indecomposables"] = True
    except ImportError:
        extra["sails_indecomposables"] = "not available"
    q = lambda x: [str(x.a), str(x.b)]
    cert = dict(D=D, family="ideal lattice (A, N/n), E = K(sqrt-2)", ideal_class=list(c), Cl_E=ed.cyc,
                n_predicted=n_pred, norm_generator=q(g),
                gram=[[q(G[i][j]) for j in range(2)] for i in range(2)], glue=[],
                zbasis=[[q(v[0]), q(v[1])] for v in basis], free=free, s_sq_lower_bound=k,
                vectors=[list(v) for v, _ in found], values=[q(v) for v in vals], checks=extra)
    path = os.path.join(HERE, "..", "data", "certificates", f"D{D}_ideal.json")
    json.dump(cert, open(path, "w"), indent=1)
    print(f"D={D}: predicted n(c)={n_pred}; CERTIFIED s_sq >= {k}; free={free}; checks={extra}; "
          f"norms={sorted(set(int(v.norm()) for v in vals))}")


if __name__ == "__main__":
    main(int(sys.argv[1]))
