# 05 — Roadmap (ordered; each task states what its outcomes would prove)

Status at hand-over: upper bound s□ ≤ 12 s PROVED (L1.5); exact cap-free evaluation of s□(M)
implemented (lattice.exact_sq_classes, L6.1); certified lower bounds S(K,2) ≥ 8 (D = 43),
≥ 6 (D = 19, 31, 57); per-edge value counts ≤ 2 in every extremal lattice; extremal lattices
touch (almost) ALL edges of the period (L6.3). The dichotomy "S(K,2) bounded vs ≍ s" is
genuinely open; the "stall at 6" claim from the other analysis is refuted (L6.2).

## R1 — DONE (09_R1_complete.md, ledger L7). Exact S(K,2), S_free(K,2) for all D ≤ 100 and
## long-period fields; acceptance met (proof written; D = 26/33/19 reproduced; 43 → 8, 31 → 6,
## 46 → 4, 58 → 6, 67 → 8; long periods 478/958/718 → 8, 862 → 12). Original task text below.
## R1 — Complete computation of S(K,2) (decisive computation; do first)
Algorithm (completeness provable):
1. Any M with s□(M) ≥ 2 contains a frame Ov + Ow with Q(v) = α, Q(w) = β indecomposable,
   distinct classes, b = B(v,w) ∈ O admissible (σ(b)² < σ(α)σ(β)); α, β may be taken on fixed
   edge-representatives modulo (O^×)² (rescale v, w by units).
2. M ⊇ Ov+Ow with [M : Ov+Ow]²·𝔳M = (αβ − b²)O (as ideals): finitely many integral overlattices;
   enumerate them (sub-O-modules of (Ov+Ow)^#/(Ov+Ow) that are B-integral).
3. Evaluate each candidate exactly with `exact_sq_classes` (cap-free).
4. Prune by L1.5/L1.2 (norm window ⇒ orthogonality) and symmetry (Galois conjugation, swap).
Acceptance: written completeness proof; reproduces imported complete values D = 26 → 2,
D = 33 → 4, D = 19 → 6; settles D = 43 (≥ 8), D = 31, 46, 58, 67.
Outcome → proof: if S(K,2) tracks s (e.g. ≈ s or ≈ values-per-period), go to R3(b);
if S(K,2) is bounded on long-period families, go to R3(a).

## R2 — Structure of extremal lattices (theory, guided by R1 data) — NOW FIRST PRIORITY
R1 data (L7.6): 20/25 extremal frames (D ≤ 100, S ≥ 4) are Galois-self-conjugate, including all
S = 8 fields (D = 43, 67, 86: det 2, 8, 8, one glue vector). Done for D ≤ 100 (L7.7): S_sym = S in 21/25
fields with S ≥ 4, not in D = 31, 89, 94, 97. Plan: (a) prove a formula/bound for S_sym via genus
theory of Q(√−d), Q(√−dD) (the provable special case); (b) describe the non-symmetric extremal
lattices of D = 31, 89, 94, 97 (what replaces the Galois symmetry — another automorph?).
Observed (EMPIRICAL, L6.4): (i) unit determinant: det = η (D = 19, also imported Δ = 33, 57, 76);
(ii) Galois-self-conjugate frames [[α, b],[b, α']] with b ∈ Z, det = N(α) − b² ∈ Z
(D = 22: det 2; D = 43: det 2; the realized L has unit volume). For (ii), E = K(√−det) =
Q(√D, √−d) is biquadratic: Cl(E) is controlled by Cl(Q(√−d)), Cl(Q(√−dD)), Cl(K) —
classical imaginary-quadratic genus theory. Task: for Galois-symmetric lattices, translate
"which indecomposables are represented" (L2.1) into conditions in Cl(Q(√−dD)) and prove
an explicit formula/bound for s□. This is the most promising PROVABLE special case.
Acceptance: theorem with proof + agreement with D = 22, 43 data.

## R3 — The main theorem
UPDATE (L7.8): a candidate family for (b) exists — det-2 Galois-symmetric frames
[[a − √D, b], [b, a + √D]], D = a² − b² − 2, glued to a unimodular M; certified s□ = 20 at D = 691, 823,
1303, 1579. Next: (1) explain s□(M_{a,b}) — which indecomposables a unimodular Galois-symmetric M
represents — via E = K(√−u) (L2.4) and the sail; (2) find sub-families with s□ → ∞ (prescribed CF,
Friesen) or a uniform bound for them.
(a) If bounded: prove a per-period cancellation — values on different edges collapse into few
    classes (D = 57: 12 values, 6 classes; pairs differ by ξ², ξ non-unit of norm ±1). Find the
    source of ξ (automorphs of the lattice / ambiguous ideal classes of E).
(b) If unbounded: construct a family with s□ → ∞ (Galois-symmetric frames over fields with
    prescribed long periods are the natural candidates; per-edge class count ≤ 2 suggests s□ ≲ 2s),
    and prove the matching lower bound; final answer would be S(K,2) ≍ s, with L1.5 as the upper
    bound. Then the right variant of question 1 is "N_r(K) ≤ c(r,d)·(sail complexity)".
~~Either way, first prove or refute per-edge ≤ 2 classes~~ — REFUTED (L7.4: D = 82, s = 1, S = 4).
s□ ≤ 2s is false; the per-edge bound is L1.4 (≤ 12). A sharper per-edge bound must use the lattice,
not just the edge (e.g. the Galois-symmetric structure L7.6).

## R4 — Verify imported claims
Thm 17 (regulator ceiling), Prop 5 (no minimum angle), Thm 4 (Ramsey staircase) — reprove or
mark unusable. Re-run the imported Thm 22 table with R1 (their search was restricted to
unit-determinant frames; D = 43 was not in it).

## R5 — Genus average via Shyr (L5.5)
Make Σ_{gen} s□ explicit using Shyr's H; compare average vs maximum from R1. Clarify proper vs
improper classes and the order O[√−δ] vs maximal order.

## R6 — Rank-1 core (L5.1)
Decide sup_K T(K) (square-full points on the sail): construct (via prescribed CF) γ with many
γξ² ∈ I(K), or bound. Determines whether s_u (mod unit squares) can ever be uniform.

## R7 — Later
Degree 3 (t_3 ≤ 3; enumerate T_δ vectors of norm ≤ 3), rank 3, Kitaoka connections.

## Do NOT do
- More winner-search plateaus on fixed reps (lower bounds only; superseded by R1).
- Geometry-of-numbers / Pareto / spherical-code bounds (L3.3, L3.5).
- Capped Fincke–Pohst class counts as final answers (use exact_sq_classes).
