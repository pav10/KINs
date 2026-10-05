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
ceiling), Thm 22 table, Prop 5, Thm 4 are IMPORTED and (except the Δ = 57 certificate) not
re-verified.

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
