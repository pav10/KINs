#!/usr/bin/env python3
"""R4 (docs/12_R4.md section 1): re-run of the imported "Thm 22" table.  The imported search was
restricted to UNIT-DETERMINANT frames: Ov + Ow with alpha = Q(v), beta = Q(w) indecomposable in
distinct square classes and det = alpha*beta - b^2 a totally positive unit.  Such a frame is
unimodular, so it has no proper integral overlattice, and the restricted maximum is
  S_unit(K) := max s_sq(frame) over unit-determinant frames   (1 if there is none).
Frames with both alpha, beta units (<1, eps> up to symmetry, s_sq = 2) are counted separately as
"trivial": the imported table evidently excluded them ("no unit-determinant frame" exactly where
only <1, eps> exists).
Every frame up to symmetry comes from complete.frames (the R1 enumeration, provably complete), and
each frame lattice is evaluated exactly (complete.evaluate).  Compared with the complete S(K,2)
of data/complete_S.jsonl.
Usage: python3 scripts/r4_thm22.py out.jsonl D1 D2 ...   (default: the D of the imported table)"""
import os, sys, json, time
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "src"))
from complete import FieldData, frames, evaluate

IMPORTED = {26: "2", 33: "4", 19: "6", 57: "6", 41: "4", 61: "4", 31: ">=4",
            **{D: "2-3" for D in (22, 34, 38, 39, 47, 55, 62, 65, 74, 79, 83)},
            **{D: "none" for D in (35, 51, 53, 66, 69, 70, 77, 78, 82, 85, 87)}}


def unit_frames_S(D):
    """Per determinant class of the unit frames: 'sq' = det in (O^x)^2 (E = K(i)) and 'eps' =
    det in eps^+ (O^x)^2 (E = K(sqrt -eps^+); the same class as 'sq' when N(eps) = -1)."""
    fd = FieldData(D)
    eta = fd.eps if fd.Neps > 0 else fd.u                  # eps^+, generator of O^{x,+}
    res = {"sq": [0, 1, None], "eps": [0, 1, None]}        # [#frames, S, best frame]
    nfr, ntriv, triv_S = 0, {"sq": 0, "eps": 0}, 1
    for key, al, be, b in frames(fd, galois=True):
        nfr += 1
        det = al * be - b * b
        if det.norm() != 1:                      # det >> 0 (admissible b); unit <=> N(det) = 1
            continue
        cls = "sq" if det.is_square_in_K() else "eps"
        if al.norm() == 1 and be.norm() == 1:  # the trivial frame <1, eps> (both units)
            ntriv[cls] += 1
            n, _, _, _ = evaluate(fd, [[al, b], [b, be]], [])
            triv_S = max(triv_S, n)
            continue
        assert cls == "sq" or (det * eta).is_square_in_K()
        n, cl, pe, basis = evaluate(fd, [[al, b], [b, be]], [])
        r = res[cls]
        r[0] += 1
        if n > r[1]:
            r[1], r[2] = n, [str(al), str(be), str(b), str(det)]
    out = dict(D=D, frames=nfr, Neps=fd.Neps, trivial_unit_frames=ntriv, S_trivial=triv_S)
    for cls, (k, n, fr) in res.items():
        out[f"unit_frames_{cls}"] = k
        out[f"S_unit_{cls}"] = n if k else None
        out[f"best_{cls}"] = fr
    k = res["sq"][0] + res["eps"][0]
    out["unit_frames"] = k
    out["S_unit"] = max(res["sq"][1], res["eps"][1]) if k else None
    return out


if __name__ == "__main__":
    out = sys.argv[1]
    Ds = [int(x) for x in sys.argv[2:]] or sorted(IMPORTED)
    S = {}
    p = os.path.join(HERE, "..", "data", "complete_S.jsonl")
    for l in open(p):
        r = json.loads(l)
        S[r["D"]] = r["S"]
    with open(out, "w") as f:
        for D in Ds:
            t = time.time()
            r = unit_frames_S(D)
            r["imported"] = IMPORTED.get(D)
            r["S_complete"] = S.get(D)
            r["sec"] = round(time.time() - t, 1)
            f.write(json.dumps(r) + "\n")
            print(D, "N(eps)=%+d" % r["Neps"], "triv:", sum(r["trivial_unit_frames"].values()), "sq:", r["unit_frames_sq"], r["S_unit_sq"], " eps:",
                  r["unit_frames_eps"], r["S_unit_eps"], " imported:", r["imported"], " S:", r["S_complete"],
                  flush=True)
