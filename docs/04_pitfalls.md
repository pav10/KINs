# 04 — Pitfalls (mathematical and computational)

P1. **Rank/PSD certification.** Never accept a Gram matrix from principal-minor tests alone (C1).
Use `gram.rank_over_K` (exact) after the 2×2 total-nonnegativity filter.

P2. **Which lattice?** For an abstract Gram realized by v_1..v_k ∈ K², the honest lattice is
L = Σ O v_i (Z-basis via `lattice.hnf` on generators v_i, √D·v_i). The pivot block
O v_p + O v_q is a sublattice and gives spurious "missing" classes (C5). Pivot coordinates of
other v_i may be non-integral — that is fine.

P3. **Square classes vs unit orbits vs elements.** Three different counts:
elements mod O^{×,+} (ι, Blomer–Kala), classes mod (O^×)² (s_u), classes mod (K^×)² (s□, κ□).
In N(ε) = +1 fields an O^{×,+}-orbit holds 2 classes mod (O^×)²; a class mod (K^×)² may hold
several orbits (C8). Report which one you count.

P4. **Representative choice in searches.** Changing a diagonal rep by a unit square is harmless
(rescales a basis vector by a unit). Changing it by ξ² with ξ a non-unit changes the lattice.
Hence fixed-rep subset searches are lower bounds (C14); complete searches must use the frame
model + overlattice enumeration (roadmap R1).

P5. **"Complete".** Only call a result complete with a written completeness argument.

P6. **Indecomposability test.** `lattice.indecomposable_test(D)`: totally positive, N(α) ≤ Δ/4
(Δ/4 = D for D ≢ 1, D/4 for D ≡ 1 mod 4), and unit-normalized form in the orbit set of
`indec.indecomposables(D)` (exact; no floats). Do not test indecomposability by small searches.

P7. **Floats.** Only allowed for (i) search-box bounds with exact re-check, (ii) diagnostics.
`QF.emb()` is float — never use it for a decision. `normalize_unit` is exact (sign tests);
the earlier float version underflowed.

P8. **Fincke–Pohst caps.** `analyze_lattice` is exhaustive only up to the trace cap; a class with
minimal representing trace above the cap is missed. Always print the cap; increase until the
class list stabilizes AND argue a bound where possible (e.g. via L1.5 the classes on a given edge
appear among minimal vectors of T_δ — a finite, cap-free check per edge; use it).

P9. **Fundamental units.** `rqf.fundamental_unit` is CF-based (fast); `fundamental_unit_bruteforce`
is only a cross-check (fails for D = 94, 151, …). D = 5 is the only purely periodic ω.

P10. **Bilinear normalization.** Q(v) = B(v,v); integral = B(M,M) ⊆ O (classical). Off-diagonal
entries in Gram matrices are B-values (not 2B).

P11. **Container hygiene (if working in ephemeral sandboxes).** Save to the repo immediately;
an earlier session lost fp.py/t2.py to a reset.

P12. **O_K = Z + Zω, not Z + Z√D.** Any Z-basis of an O-lattice built from generators must use
v and ωv (C17). For D ≡ 1 (mod 4) Z[√D] has index 2 in O.

P13. **Huge fundamental units.** Balanced representatives can still be lopsided by a factor up to ε₁
(ε ≈ 3.4·10^16 for D = 958). Never search boxes in embedding coordinates; use exact trace forms
(complete.admissible_b: Tr(b²/αβ) < 2) and exact norm/CF tests (is_principal via Serret). Float
embeddings a − b√D cancel catastrophically for lopsided elements (use N(x)/σ_large).

P14. **Shell hygiene for background runs.** `cmd1 && VAR=... && nohup job &` backgrounds the whole
chain, so VAR is unset afterwards; and `pkill -f pattern` kills the shell whose command line
contains the pattern (use an anchored `^python3 ...`).
