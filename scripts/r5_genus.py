#!/usr/bin/env python3
"""R5 (docs/13_R5.md): the genus G of <1, eps> (= (O_E, N), E = K(sqrt -2)) for a prime p = 3 (mod 4).
Theorem 13.1 says G = {M_A : [A] in A_p}, with A_p = ker(Cl(E) -> Cl+(K)) of order h(-8p)/2.
Combined with Theorems 10.2(c) and 11.3, the sets of indecomposable square classes represented by the
lattices of G are exactly the sets
    S(f) = {[gamma] : f represents the norm of a member of [gamma]},   f in Cl(-8p)^2, f ~ f^{-1}.
Independent check (CERTIFIED for the primes p < 100):
  * every lattice of G is unimodular, hence a maximal integral lattice, and it contains a frame as soon
    as s_sq >= 2 (under (H) s_sq is even, so s_sq > 0 already forces a frame);
  * so R1's complete enumeration (complete.frames + maximal isotropic glue, Galois reduction OFF)
    meets every M in G with s_sq > 0.  Collect the represented-class sets of all unimodular lattices of
    det class 2 = eps found by R1 and compare them with the sets S(f) of the form model.
Usage: python3 scripts/r5_genus.py p1 p2 ...      (prints the comparison)"""
import os, sys, time
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "src"))
sys.path.insert(0, HERE)
from rqf import QF, same_square_class
from lattice import Bform, exact_sq_classes
from complete import FieldData, frames, DiscModule, overlattice_basis
from r3_tools import Sieve, field, qf_reduce, h8
from certify_form import class_forms


def model_sets(p):
    """Square-class reps of the indecomposables and the sets S(f) (one per pair f, f^{-1}) with f in
    Cl(-8p)^2 and S(f) nonempty."""
    sv = Sieve(4 * p)
    rec = field(p, sv, need_reps=True)
    reps = [QF(p, x << tw, y << tw) for ((x, y), tw, N, inE) in rec["reps"]]
    byf = {}
    for cid, d in enumerate(class_forms(p, rec, sv)):
        for f in d:
            byf.setdefault(f, set()).add(cid)
    sets = {}
    for f, S in byf.items():
        g = qf_reduce((f[0], -f[1], f[2]))
        assert byf.get(g) == S                      # nu(f) = nu(f^{-1}): conjugate ideals
        sets[min(f, g)] = frozenset(S)
    return reps, sets


def r1_sets(p, reps):
    """Represented-class sets of all unimodular lattices of det class 2 found by the complete R1
    enumeration (all frames, no Galois reduction), with the square classes named by model ids."""
    fd = FieldData(p)
    Delta = 4 * p
    seen, nlat = set(), 0
    for key, al, be, b in frames(fd, galois=False):
        G = [[al, b], [b, be]]
        det = al * be - b * b
        if not (det * 2).is_square_in_K():          # det class 2 (= eps under (H))
            continue
        A = DiscModule(fd.ar, fd.ar.from_qf(al), fd.ar.from_qf(b), fd.ar.from_qf(be))
        for gens in A.maximal_isotropic():
            basis = overlattice_basis(fd, G, gens)
            T = [[(Bform(G, u, v)).trace() for v in basis] for u in basis]
            dT = _det4(T)
            if dT != Delta * Delta:                 # N(vol M) = det(Tr B) / Delta^2 must be 1
                continue
            nlat += 1
            cl, _ = exact_sq_classes(p, G, basis)
            ids = []
            for c in cl:
                m = [i for i, r in enumerate(reps) if same_square_class(c, r)]
                assert len(m) == 1, (p, c)
                ids.append(m[0])
            seen.add(frozenset(ids))
    return seen, nlat


def _det4(M):
    from fractions import Fraction as F
    A = [[F(x) for x in r] for r in M]
    n, d = len(A), F(1)
    for c in range(n):
        r = next((r for r in range(c, n) if A[r][c] != 0), None)
        if r is None:
            return 0
        if r != c:
            A[c], A[r] = A[r], A[c]
            d = -d
        d *= A[c][c]
        for r in range(c + 1, n):
            f = A[r][c] / A[c][c]
            A[r] = [x - f * y for x, y in zip(A[r], A[c])]
    return d


def compare(p):
    t = time.time()
    reps, sets = model_sets(p)
    model = set(sets.values())
    r1, nlat = r1_sets(p, reps)
    return dict(p=p, h8=h8(p), A=h8(p) // 2, pairs_hit=len(sets), distinct_model_sets=len(model),
                r1_unimodular_lattices=nlat, distinct_r1_sets=len(r1), equal=(r1 == model),
                nu_star=max((len(s) for s in model), default=0), sec=round(time.time() - t, 1))


if __name__ == "__main__":
    for p in map(int, sys.argv[1:]):
        print(compare(p), flush=True)
