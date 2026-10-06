#!/usr/bin/env python3
"""Per-edge statistics over ALL maximal integral overlattices of all frames (R1), D given:
max #values and max #classes on one edge, overall and among lattices with s_sq >= 4.
Usage: python3 scripts/edge_stats.py out.jsonl D1 [D2 ...]"""
import os, sys, json
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from complete import (FieldData, frames, DiscModule, overlattice_basis, edge_values,
                      class_count, same_square_class)

out = sys.argv[1]
for D in map(int, sys.argv[2:]):
    fd = FieldData(D)
    ar = fd.ar
    st = dict(D=D, s=len(__import__("lattice").edge_deltas(D)), max_vals=0, max_cls=0,
              max_vals_ge4=0, max_cls_ge4=0, n=0, witness_cls=None)
    for key, al, be, b in frames(fd):
        A = DiscModule(ar, ar.from_qf(al), ar.from_qf(b), ar.from_qf(be))
        G = [[al, b], [b, be]]
        for _, gens, m in A.isotropic_submodules():
            if not m:
                continue
            st["n"] += 1
            ev = edge_values(fd, G, overlattice_basis(fd, G, gens))
            allv = [v for e in ev for v in e]
            ssq = class_count(allv)
            mv = max(len(e) for e in ev)
            mc = max(class_count(e) for e in ev)
            if mc > st["max_cls"]:
                st["witness_cls"] = dict(frame=[str(al), str(be), str(b)], gens=gens, s_sq=ssq,
                                         per_edge_vals=[[str(v) for v in e] for e in ev])
            st["max_vals"], st["max_cls"] = max(st["max_vals"], mv), max(st["max_cls"], mc)
            if ssq >= 4:
                st["max_vals_ge4"] = max(st["max_vals_ge4"], mv)
                st["max_cls_ge4"] = max(st["max_cls_ge4"], mc)
    with open(out, "a") as f:
        f.write(json.dumps(st) + "\n")
    print(D, {k: v for k, v in st.items() if k != "witness_cls"}, flush=True)
