# 03 — Corrections, retractions, dead ends (read before restating anything)

C1. **"3×3 principal minors = 0 (plus 2×2 minors ≥ 0, diagonal ≻ 0) certifies rank ≤ 2 / totally
PSD."** FALSE at size ≥ 4: (3/2)I − (1/2)J passes and is rank 4, indefinite. Correct test:
exact rank over K ≤ 2 (sufficient given the 2×2 conditions). [tests: test_rank_pitfall_matrix]

C2. **Square-class coverage.** One rep per O^{×,+}-orbit undercounts when N(ε) = +1: each orbit
can carry two classes [α], [αη]. Seed α and αη, then dedup by class. [test_square_class_counts]

C3. **"Sum of two squares / diagonal forms are good witnesses."** Tautological: diagonal lattices
represent ≤ r classes (L0.5).

C4. **"The spectral/Ramsey spherical-code writeup gives a uniform bound."** Spectral part vacuous
for all totally real fields; Ramsey constant blows up with Δ. Abstract-Gram bounds also ignore
realizability.

C5. **First "D = 31 counterexample" (2 of 4 classes missing).** Artifact: trace lattice was built on
the pivot sublattice Ov_p + Ov_q instead of L = ΣOv_i. Fixed; the corrected computation still
refutes linear irreducibility (L3.1), with different numbers.

C6. **"N_r grows with D" / "parametrize by D".** Wrong parameter (D = m²+1 has period 1). Then
**"parametrize by ι = Σu_odd"** — also wrong (L3.6). The bound parameter is the period s (L1.5).

C7. **"The best form hugs Θ(s) faces (D = 31 touches 3/4 faces)."** Over-read of one
non-certified example; corrected enumeration gives a plateau (2,3,4,4,4 orbits for
s = 2,4,6,8,12). Lean reversed: data favour boundedness. Lower-bound construction via convergents
failed (L5.2).

C8. **"Every indecomposable in a square class is a unit multiple of one balanced rep (norm ceiling
kills non-unit square factors)."** Unproved and contradicted by data: κ□ < 2ι for D = 19, 31
(N(ε) = +1) means some classes contain several unit orbits (αξ², ξ non-unit).

C9. **WRONG RETRACTION (now reversed).** It was claimed that "totally positive units need not be
indecomposable (Maass, Q(√5))" and Cor "units ≤ rank" was withdrawn. This was wrong:
Q(√5) has no totally positive non-square unit (N(ε) = −1), "sum of squares" ≠ "decomposable",
and every totally positive unit IS indecomposable (L0.4, 3-line proof). L1.3 stands.

C10. **"For extremal (δ = ε⁺) lattices the wall moves from Siegel–Linnik to a finite 2-adic
ray-class problem."** Overclaim. Unramified-outside-2 makes the LOCAL (genus/solvability)
conditions 2-adic congruences, but the class condition is in Cl(E) with h⁻(E) = Δ^{1/2+o(1)}
(Brauer–Siegel). Both layers remain (L2.4).

C11. **"Σ_{genus} s□(M') ≤ κ□(K)."** Wrong: a class can be represented by several genus classes;
correct bound carries a divisor-type factor (L5.5).

C12. **"Per-face bound ⇒ s□ ≤ 12·f_K with f_K = faces mod O^{×,+}."** Off by the index
[O^{×,+} : (O^×)²]; correct statement s□ ≤ 12 s (L1.5). Similarly "s = 2 ⇒ ≤ 12" → ≤ 24.

C13. **Statuses downgraded.** "Pareto orbits = Θ(s)" is HEURISTIC (L3.5). "Linear-irreducibility
failure is explained by η-twins of large class-trace" is HEURISTIC. Imported Thm 17 (regulator
ceiling), Thm 22 table, Prop 5, Thm 4 were IMPORTED and (except the Δ = 57 certificate) not
re-verified. [Superseded by R4 (docs/12_R4.md): Thm 4, Prop 5 PROVED; Thm 22 re-run (C26); Thm 17 unusable as stated (C28).]

C14. **Search completeness.** `find_winner` / `largest_gram` on fixed class reps give LOWER bounds
for S(K,2), not S(K,2): classes reachable only through reps αξ² with ξ a non-unit are not explored
(rescaling a vector by a unit is harmless). E.g. our D = 33 search finds 2, imported complete
search finds 4; D = 19 ours 4 (winner on balanced reps), imported 6.

C15. **"At most 2 values per edge in every lattice computed"** (stated at hand-over). False in
general: ⟨1, γ⟩ in odd-period fields gives 3–4 values on one edge (D = 2, 10, 17, 26, 37, 41, 65).
True so far only for the extremal (s□ ≥ 4) lattices, and the right quantity is classes per edge.

C16. **"Stall at 6" (other analysis, 40 fields, Δ ≤ 348).** Their search was restricted to
unit-determinant frames and did not include D = 43, where S ≥ 8 (L6.2).

C17. **`lattice.realized_module` spanned Z[√D]·v_i, not O·v_i** (generators v_i, √D·v_i). For
D ≡ 1 (mod 4) that is an index-4 sublattice which is not an O-module. Fixed (generators v_i, ω v_i).
Affected (all lower bounds, but wrong as values): 07_data table C rows D = 5, 13, 17, 21, 29, 33,
37, 41, 53, 57, 61, 65 (e.g. D = 13, 21, 29, 53: 1 → 2; D = 57 orbit winner 4 → 6); C15's ⟨1, γ⟩
edge counts (D = 65: 4 values, not 3); the D = 33 orbit-plateau entry. Not affected: L3.1 (D = 19,
31; D = 57 used an explicit O-basis), L6.2, frame57, D ≢ 1 data. [test_realized_module_is_OK_module]

C18. **`lattice.edge_deltas` used a fixed 254-term continued fraction** — failed for period s > 60
(first at D > 300). Fixed (length from s). Found by the sails cross-check.

C19. **"Per-edge class count ≤ 2" (R3 working hypothesis, "would give s□ ≤ 2s") and "≤ 2 values per
edge in every lattice with s□ ≥ 4" (L6.3/C15).** Both REFUTED (L7.4): D = 82 has s = 1 and S = 4.

C20. **Do not compute S_free from maximal lattices only.** A free lattice may lie only in non-free
maximal lattices; R1 evaluates lattices maximal among FREE overlattices (L7.1). (Caught before
any value was recorded.)

C21. **"The det-2 family plateaus at 28, so 28 may be a ceiling" (L7.8 discussion).** The plateau is a
property of the Galois-symmetric slice only. The full det-class-2 maximum max_c n(c) (L8.3) exceeds it:
S ≥ 36 at D = 12511, 17539 (CERTIFIED), and it keeps growing for larger D (L8.5, L8.7).

C22. **C7's "lean toward boundedness" (already weakened by L7.8) is not supported.** Under (H),
S(K,2) ≥ max_c n(c) ≥ n̄ (L8.3), and both grow in the data. The dichotomy (L5.2) is still OPEN, but
the evidence now points to unbounded, slow growth.

C23. **Upper bounds from non-classical hulls are not attained.** For E = K(i), the bound 10 at D = 43
exceeds the complete value S = 8, and classical sublattices lose most classes (L8.6). Always state
which E and which integrality the bound refers to.

C24. **"S(K,2) ≥ N_E(K)/h(E), and growth needs N_E/h(E) → ∞" (L8.3, L8.7, roadmap R3) is the wrong normalisation.**
- n(c) lives on a subgroup of order h(−8D)/2, so the bound is 2Σ_c n(c)/h(−8D): h(K) drops out (L9.1).
- The N_E pigeonhole is a poor proxy (corr 0.60 with the truth). The multiplicity Σ_{m|N}χ_{−8}(m) matters, and its mean ∝ L(1, χ_{−8p}) largely cancels h(−8p).
- The governing quantity is ι/√p (L9.6), not "ι against h(−8D)". The L8.7 heuristic "n̄ ≍ (norm-indecomposables)/h(−8D)" is superseded.

C25. **Theorem 11.3 is for D prime only.**
- For composite D ≡ 3 mod 4, h(K) is even, N_{E/k2} is not injective on A, and ν(f) only bounds n(c) from above (12/12 sampled composite fields differ).
- Do not use the form model there; use the PARI scan (r2_scan.py).

C26. **Imported Thm 22 table, D = 22 ("2–3").** Wrong: a unit-determinant frame represents 4 classes, which equals S(Q(√22),2) (R4, L10.1). The table's "no unit-determinant frame" means no frame other than the trivial ⟨1, ε⟩, which exists whenever N(ε) = +1.

C27. **"(a)∩(b) is the irreducible core of question 1, by the Ramsey reduction" (01_problem).** Overstated.
- The Ramsey reduction (L3.4) holds for any set of distinct indecomposables and is empty for d = 2.
- The staircase it produces in degree d is not a real quadratic sail.
- No implication from the real quadratic binary case to degree d follows. The real quadratic binary case remains the natural first case, not a proven core.

C28. **Imported Thm 17 (regulator ceiling) as a usable bound.** Its error term E is unspecified, so as stated it is not a theorem. Use the PROVED s□ ≤ 12s ≤ 24 + 24.94·(Reg − log(√D − 1)) instead (L1.6).

C29. **"Isometry classes in the genus of a free lattice ↔ Shyr's H" (L2.2) without qualification.**
- Shyr's count is for lattices whose Clifford order is the maximal order O_E. A classical lattice has Tr C₀(M) ⊆ 2O_K (Proposition 13.4).
- So for E = K(i), D ≡ 3 mod 4 (det class 1), no classical lattice has Clifford order O_E; Picard groups of non-maximal orders are needed there.
- In det class 2 under (H) the count is h(−8D)/2 (Theorem 13.1).


C30. **"L1.4 (≤ 12 per edge) is the only per-edge bound known" (L7.4) and the constant 12 in s□ ≤ 12s (L1.5, L1.6, 01_problem).** Not sharp for the lattices of this repo.
- They are classical (B(M, M) ⊆ O), so T_δ is an integral lattice and its norm-1 vectors are orthonormal.
- Hence ≤ rd values per functional, so S(K,2) ≤ 4s, and S(K,2) ≤ 3s when N(ε) = +1 (Theorem 15.1, Corollary 15.2). Equality holds at D = 82.
- The kissing-number bound ½τ_{rd} applies only to non-classical lattices (Q(M) ⊆ O, B(M, M) ⊆ ½O).
- Likewise T(K) ≤ 3s (R6) becomes T(K) ≤ 2s.

C31. **"#{v : T_δ(v) ≤ t_d} is not uniform in the lattice" (L5.6).** Wrong. For any integral lattice of rank n, #{v : T(v) ≤ t} ≤ (2√t + 1)^n, by packing (Lemma 15.4). In degree d the field enters only through t(K) and the number of functionals Φ_t(K) (Theorem 15.5).

C32. **"t₃ ≤ 3 (IMPORTED Kala–Tinková)" (L5.6, docs/06_references.md).** Misattributed and false.
- Kala–Tinková treat simplest cubic fields, where min traces are 1 or 2. Tinková (2105.04204) shows min traces are unbounded in families of cubic orders.
- For maximal orders: t(K) = 4 at disc 3137 and t(K) = 5 at disc 14229 (complete scan to disc 20000).
- t(K) ≥ a in Thomas's family whenever Z[ρ] = O_K (Theorem 15.6). CERTIFIED for a ≤ 100, so sup t(K) ≥ 100.
- Whether t₃ = ∞ over maximal orders is HEURISTIC; the gap is infinitely many monogenic parameters.
