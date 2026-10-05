# 01 — Problem, definitions, scope

## Objects
- K totally real, degree d (=m); O = O_K; σ_1..σ_d real embeddings; α ≻ 0 totally positive;
  α ⪯ β iff β−α ⪰ 0 (totally nonnegative). O⁺ = totally positive integers.
- α ∈ O⁺ is **indecomposable** if α ≠ β+γ with β,γ ∈ O⁺. I(K) = set of indecomposables.
  Units: every totally positive unit is indecomposable (PROVED, ledger L0.4).
- M a totally positive definite O-lattice of rank r with quadratic map Q, bilinear B,
  scale 𝔞 = B(M,M), volume/discriminant 𝔳M. "Classical/integral": B(M,M) ⊆ O.
- s□(M) = #{α(K^×)² : α ∈ Q(M) ∩ I(K)}  — the **primary invariant**.
- s_u(M) = same count mod (O^×)² (finer; conjecturally NOT uniformly bounded, ledger L5.1).
- N_r(K) = sup_M s□(M) (sup over rank-r integral M); S(K,2) := N_2(K) for real quadratic K.
- κ□(K) = #(I(K) modulo (K^×)²) (finite). ι(K) = #(I(K)/O^{×,+}).
- For real quadratic K = Q(√D): s = period length of CF of ω; faces of the sail
  S_K = ∂Conv(O⁺) (Klein sail): vertices = totally positive convergent elements, edges carry
  the semiconvergents; every indecomposable α has a totally positive δ in the codifferent
  O^∨ with Tr(δα) = 1, shared by all indecomposables on one edge (Kala–Tinková 2005.12312,
  Prop. 3.1).

## Questions (priority order, as set by Pavlo)
1. **Uniformity**: is N_r(K) ≤ c(r,d) for all totally real K of degree d? If false, the right
   variant (bound in terms of r, d and a field invariant — which one?).
2. Special cases: (a) binary lattices (Gauss/Clifford composition ↔ ideal classes of the CM
   extension E = K(√−det)); (b) real quadratic fields (all indecomposables explicit via CF).
   **Step 1 of this repo = (a)∩(b).**

## Why (a)∩(b) is the core (not just the easiest case)
Ramsey reduction (ledger L3.4, IMPORTED): a large s□(M) in any degree forces a long
*two-place staircase* of represented indecomposables — the configuration that the real
quadratic sail realizes. So bounding staircases for binary lattices over real quadratic fields
is the irreducible core of question 1.

## Current one-line state
s□(M) ≤ 12·s (s = CF period) is PROVED (per-face lemma, L1.5), and s□(M) is computable exactly
(L6.1). Certified: S(Q(√43),2) ≥ 8 (s = 10); extremal lattices meet (almost) every edge of the
period with ≤ 2 values per edge. Whether S(K,2) is bounded or grows ≍ s is OPEN; the uniform bound
is equivalent to a bounded-fibre statement for indecomposable ↦ ideal class of E = K(√−det M)
(L2.3), which no known technique decides. See 02_ledger.md, 05_roadmap.md.
