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
11. R2 session: the det-2 family lattices are ideal lattices of E = K(√−2). For D ≡ 3 mod 4,
   Tr O_E = 2O_K, so every det-class-2 lattice lies in an integral ideal lattice. Under (H), s□ equals
   the number of indecomposables whose norm-ideals meet one class of Cl(E) (Theorem 10.2), so the best
   lattice comes from a PARI class-group computation instead of a search. Scans: all D ≡ 3 mod 4,
   D ≤ 20000, and all (H)-fields up to 10^5. Certified records 36 → 40 → 48 → 52 (D = 98179), all free and
   all D prime. det class 1 (K(i)) is weak for classical lattices.
12. R3 session.
   - The R2 class distribution is elementary for p prime: n(c) on Cl(K(√−2)) equals ν(f) on the form classes of discriminant −8p (Theorem 11.3). Checked against PARI for all 4808 primes ≤ 10^5.
   - This gives a pure-Python scan of all primes ≤ 10^6 (≈ 1 h; the PARI route took hours to 10^5), class-group-free certificates, and records 56, 60, 64 below 10^6, then 76 at p = 7956859 (top-200 sample of (10^6, 10^7]).
   - Mean ν* grows like 2.8·ln p, and the average over A is ∝ ι/√p.
   - The requested proof of growth along a family was not achieved. It reduces to (G1) a class-number-one-type input and (G2) a twisted-divisor mean value (Proposition 11.5).
13. R4 session. The imported claims were re-checked from their ledger summaries; the source analyses are not in the pack.
   - Thm 22 table: reproduced with R1's complete frame enumeration, except D = 22 (4, not 2–3).
   - Thm 4 and Prop 5 are now PROVED (Prop 5 via an explicit family over Q(√(n²+2))).
   - Thm 17 is unusable as stated; replaced by the proved ceiling s□ ≤ 12s ≤ 24 + 24.94·Reg.
   - The "irreducible core" claim of 01_problem is withdrawn (C27).
14. R5 session.
   - The genus of ⟨1, ε⟩ (D ≡ 3 mod 4, (H)) is the set of ideal lattices over A = ker(Cl E → Cl⁺K): h(−8D)/2 proper classes, all with O⁺ = {±1}.
   - The genus sum of s□ is Σ_c n(c), so the genus average is R3's "avg".
   - The model is checked against R1's complete enumeration for the 30 primes below 200.
   - The Clifford order is never maximal for classical det-class-1 lattices (C29).
   - The maximum exceeds the average by 1–4×.

