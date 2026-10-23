# 05 — Roadmap (ordered; each task states what its outcomes would prove)

Status at hand-over: upper bound s□ ≤ 12 s PROVED (L1.5) [R7: ≤ 4s, sharp]; exact cap-free evaluation of s□(M)
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

## R2 — DONE (docs/10_R2.md, ledger L8): for D ≡ 3 mod 4 with (H), the indecomposables represented by a
## det-class-2 lattice are the norms of ideals in one class of Cl(K(√−2)); s□ = n(c) exactly; S(K,2) ≥
## max_c n(c) ≥ n̄; records 36 (D ≤ 20000) and beyond (D ≤ 100000). Original text kept below.
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

## R3 — WORKED (docs/11_R3.md, ledger L9). Not proved: growth along a family (L9.7 — needs (G1) long
## principal cycles, i.e. class-number-one-type input, and (G2) a twisted-divisor mean value).
## Delivered: h(K)-free reduction (L9.1), exact elementary model ν on Cl(−8p) for p prime (L9.3),
## PARI-free certificates (L9.4), records 56/60/64 up to 10^6 and the 10^7 range (L9.5/L9.8),
## mean ν* ∝ ln p (L9.6). Next for R3: (a) conditional theorem under explicit hypotheses
## (G1)+(G2) with constants; (b) prove (G2) on average over p (large sieve over the family of
## quadratic polynomials (p − u²)/Q); (c) other det classes / composite D (C25).
## R3 — The main theorem (original text)
UPDATE (L8): R3(b) is now an analytic question. Under (H), S(K,2) ≥ N_E(K)/h(E), where N_E(K) is the
number of indecomposable square classes whose ideal is a norm from E = K(√−2), and h(E) ≍ h(−8D).
(1) Prove N_E(K)/h(E) → ∞ along an explicit family (prescribed continued fractions with large partial
quotients against class numbers of Q(√−2D); Friesen-type constructions), or find an obstruction.
(2) Determine max_c n(c) versus n̄ (the "best class" excess, about 3–4×).
(3) Other det classes and D ≢ 3 (mod 4): Lemma 10.1 changes, so redo (a)–(c).
UPDATE (L7.8): a candidate family for (b) exists — det-2 Galois-symmetric frames
[[a − √D, b], [b, a + √D]], D = a² − b² − 2, glued to a unimodular M; certified s□ = 20 at D = 691, 823,
1303, 1579, 24 at D = 3931, 28 at D = 5419; the family record then plateaus at 28 for D ≤ 20000
(L7.8). Next: (1) explain s□(M_{a,b}) — which indecomposables a unimodular Galois-symmetric M
represents — via E = K(√−u) (L2.4) and the sail; (2) given the plateau, first ask whether
s□(M_{a,b}) ≤ 28 is a theorem; otherwise find sub-families with s□ → ∞ (prescribed CF, Friesen).
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

## R4 — DONE (docs/12_R4.md, ledger L10): Thm 22 table re-run (D = 22 corrected, C26); Thm 4 and
## Prop 5 PROVED; Thm 17 unusable as stated (clean ceiling PROVED, C28); Cor 21.3 trivial in lattice form.
## R4 — Verify imported claims (original text)
Thm 17 (regulator ceiling), Prop 5 (no minimum angle), Thm 4 (Ramsey staircase) — reprove or
mark unusable. Re-run the imported Thm 22 table with R1 (their search was restricted to
unit-determinant frames; D = 43 was not in it).

## R5 — DONE for the genus of ⟨1, ε⟩ (docs/13_R5.md, ledger L11): classes ↔ ker(Cl E → Cl⁺K),
## h(−8D)/2 proper classes, O⁺ = {±1}; Σ_gen s□ explicit (Theorem 13.3); trace obstruction for the
## Clifford order (C29); average vs maximum tables. Open: genera of other det classes / non-maximal orders.
## R5 — Genus average via Shyr (L5.5) (original text)
Make Σ_{gen} s□ explicit using Shyr's H; compare average vs maximum from R1. Clarify proper vs
improper classes and the order O[√−δ] vs maximal order.

## R6 — DONE (docs/14_R6.md, ledger L12): sup_K T(K) = ∞ PROVED (D = n² + 3, T ≫ log D); T(K) ≤ 3s (R7: 2s);
## exact witness test for indecomposability; T₁ in O(s); records T₁ ≈ 3.3 D^{1/4} on h = 1 fields.
## Open: unconditional polynomial growth of T; the D ≡ 1 mod 4 analogue of Lemma 14.1.
## R6 — Rank-1 core (L5.1) (original text)
Decide sup_K T(K) (square-full points on the sail): construct (via prescribed CF) γ with many
γξ² ∈ I(K), or bound. Determines whether s_u (mod unit squares) can ever be uniform.

## R7 — DONE (docs/15_R7.md, ledger L13): sharp per-face lemma (≤ rd per functional for classical M, equality
## iff Z^{rd}) ⇒ S(K,2) ≤ 4s, attained at D = 82; degree-d per-functional counts uniform (C31); "t₃ ≤ 3" REFUTED
## (C32), t₃ ≥ 100 CERTIFIED via Thomas's family (Theorem 15.6); Kitaoka: u_E ≤ 2r − 1 for universal classical
## rank r, and indecomposables cannot decide Kitaoka (D = n² − 1).
## R7 — (original text)
Degree 3 (t_3 ≤ 3; enumerate T_δ vectors of norm ≤ 3), rank 3, Kitaoka connections.

## R8 — Next candidates (each states what its outcome would prove)
- (a) **Per-edge ≤ 2 when s is even** (data: yes for every maximal lattice, D ≤ 100). A proof gives S(K,2) ≤ 2s
  for N(ε) = +1, and would locate the extremal fields among N(ε) = −1. A counterexample needs a classical
  binary M with (M, T_δ) ⊇ Z³ — search the Z³ ⊥ ⟨m⟩ shapes directly.
- (b) **t₃ = ∞ unconditionally.** Either a family with provably infinitely many monogenic members (e.g. a
  sub-family of Thomas's where index 1 follows from congruences), or an O_K-version of the box argument of
  Theorem 15.6(i).
- (c) **Degree-3 analogue of R1** for small fields: exact max s□ over binary lattices over Q(ζ₇)⁺ (disc 49), using
  Theorem 15.5 with t = 2 (Φ₂ explicit) as the evaluation step.

## Do NOT do
- More winner-search plateaus on fixed reps (lower bounds only; superseded by R1).
- Geometry-of-numbers / Pareto / spherical-code bounds (L3.3, L3.5).
- Capped Fincke–Pohst class counts as final answers (use exact_sq_classes).
