# 08 — How we got here (short)

1. Code review of Sage modules: rank certification bug (C1), square-class coverage bug (C2).
   Pure-Python reimplementation validated against the Sage outputs.
2. Spherical-code/Ramsey writeup assessed: no uniform bound (C4). O'Meara framework: acute lemma,
   per-field finiteness, mod-2M route dies for class counts (L3.2). Yamamoto assessed.
3. Plan: Tier 0/1 done; Tier 2 = linear-irreducibility hinge; Tier 3 = Gauss/Shyr.
4. Tier 2: after fixing the pivot-sublattice artifact (C5), linear irreducibility REFUTED (L3.1).
5. Parametrization corrected twice: not D, not ι (C6); period s via Kala–Tinková faces.
   Per-face lemma (L1.4) ⇒ s□ ≤ 12s (L1.5).
6. Lower-bound attempt (hug Θ(s) faces via convergents) failed; plateau data led to a lean toward
   boundedness (later weakened, C7, L6).
7. Pareto route: Θ(s), cannot decide (L3.5). Uniformity ⇔ bounded fibres of indecomposable ↦
   ideal class of E (L2.3); analytic wall below Duke/Linnik.
8. Two external registry write-ups merged: value-minimality, representation-number bound, pincer,
   rank-1 core, frame criterion, unit-determinant pattern, "stall at 6".
9. Hand-over session: exact cap-free s□ via per-edge minimal vectors (L6.1); S(Q(√43),2) ≥ 8
   (L6.2) refutes the stall; extremal lattices meet ~all edges (L6.3); my two retractions corrected
   (C9 units ARE indecomposable; C10 ray-class relocation overclaim); C15.
10. R1 session: complete search (frames → maximal integral overlattices = maximal isotropic
   submodules of O²/GO² → exact per-edge evaluation), proof in 09. Exact S(K,2) for D ≤ 100
   (max 8 at D = 43, 67, 86) and long-period fields; S_free differs when h > 1; "≤ 2 classes per
   edge" refuted (D = 82). Bugs fixed: realized_module used Z[√D] (C17), edge_deltas s > 60 (C18).
   Independent cross-check of indecomposables/facets against github.com/pav10/sails (D ≤ 3000).
