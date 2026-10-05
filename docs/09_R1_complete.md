# 09 — R1: complete computation of S(K,2) (binary lattices, real quadratic K)

Status: **PROVED** (Theorem 1, written below) + **CERTIFIED** (values, `scripts/complete_S.py`,
cross-checked as in §4).  Code: `src/complete.py`.  Data: `data/complete_S.jsonl`,
`data/edge_stats.jsonl`, `data/complete_S_long.jsonl`.

Conventions as in 01/02: K = Q(√D), O = O_K = Z + Zω, M a totally positive definite O-lattice of
rank 2 (finitely generated O-module, **not necessarily free**) in a binary quadratic space (V, B, Q)
over K, Q(x) = B(x, x); *integral* = classical, B(M, M) ⊆ O.
S(K,2) = max s□(M) over all such M; S_free(K,2) = the same over free M (= binary forms).

## 1. The theorem

Let R be a set of representatives of I(K)/(O^×)² (indecomposables modulo squares of units;
|R| = ι·[O^{×,+} : (O^×)²] = ι or 2ι).  A **frame** is f = (α, β, b) with α, β ∈ R in distinct
classes mod (K^×)², b ∈ O, αβ − b² ≻ 0.  F = set of frames (finite).  L_f := O e1 ⊕ O e2 with
Gram G_f = [[α, b], [b, β]].  Max(f) := maximal integral O-lattices M ⊇ L_f in K L_f;
Max_free(f) := lattices maximal among the free integral O-lattices M ⊇ L_f.

**Theorem 1 (PROVED).** If κ□(K) ≥ 2 then
  S(K,2) = max_{f ∈ F} max_{M ∈ Max(f)} s□(M),   S_free(K,2) = max_{f ∈ F} max_{M ∈ Max_free(f)} s□(M);
if κ□(K) = 1 then S = S_free = 1.  F and every Max(f), Max_free(f) are finite and explicitly
computable, so S(K,2) and S_free(K,2) are computable exactly.

### Proof

**Step 1 (frames).** Let M be integral with s□(M) ≥ 2; pick v, w ∈ M with Q(v), Q(w) ∈ I(K) in
distinct square classes.  I(K) is stable under multiplication by O^{×,+} ⊇ (O^×)² (a totally
positive unit acts as an automorphism of the additive monoid O⁺), so there are units μ, ν with
μ²Q(v) = α ∈ R, ν²Q(w) = β ∈ R; replace v, w by μv, νw ∈ M (Q(μv) = μ²Q(v)).  v, w are
K-independent (w = λv would give β = λ²α, same class).  b := B(v, w) ∈ O (integrality), and
αβ − b² ≻ 0 by the strict Cauchy–Schwarz inequality in the positive definite planes V ⊗_σ R.
So f = (α, β, b) ∈ F and e1 ↦ v, e2 ↦ w is an isometry of L_f onto Ov + Ow ⊆ M.  F is finite:
R is finite (Dress–Scharlau / Blomer–Kala) and b lies in a bounded region of O ⊂ R².

**Step 2 (sandwich).** For L = L_f put L^# = {x ∈ V : B(x, L) ⊆ O} = G_f^{-1} O².  An integral
M ⊇ L satisfies M ⊆ L^#; [L^# : L] = [O² : G_f O²] = N(αβ − b²) < ∞.  Hence the integral lattices
containing L are finitely many, every one lies in a maximal one, and "maximal among integral
lattices containing L" = "maximal integral lattice of V" (any integral lattice containing such
an M contains L).

**Step 3 (monotonicity).** M ⊆ M′ ⇒ Q(M) ⊆ Q(M′) ⇒ s□(M) ≤ s□(M′).  With Steps 1–2 (transported
by the isometry): every M with s□ ≥ 2 has s□(M) ≤ s□(M′) for some M′ ∈ Max(f), f ∈ F; conversely
every M′ ∈ Max(f) is an admissible lattice.  If M is free it lies (Step 2, finiteness) in a lattice
maximal among the free integral lattices containing L, i.e. in Max_free(f).  If κ□ ≥ 2 then
S ≥ S_free ≥ 2 (⟨1⟩ ⊥ ⟨γ⟩ is free), so the lattices with s□ ≤ 1 do not matter.  ∎(reduction)

**Step 4 (overlattices = isotropic submodules).** On A := L^#/L ≅ O²/G_f O² (x = G_f^{-1} y) the
pairing b_A(x, y) = B(x, y) mod O is well defined (B(L^#, L) ⊆ O) and O-bilinear;
b_A(y, y′) = yᵀ G_f^{-1} y′ mod O.  M ↦ M/L is an inclusion-preserving bijection between integral
lattices L ⊆ M ⊆ L^# and O-submodules S ⊆ A with b_A(S, S) = 0 (M integral ⇔ B(x, y) ∈ O on M).

**Step 5 (primary decomposition).** A = ⊕_p A_p (p-primary parts, rational primes p | N(αβ−b²)),
each an O-submodule; for x ∈ A_p, y ∈ A_q, p ≠ q: p^a b_A(x,y) = 0 = q^c b_A(x,y), so b_A(x,y) = 0.
Every submodule is the sum of its p-parts, and is isotropic iff each part is.  So the isotropic
submodules of A are the products of those of the A_p, and the maximal ones are the products of
maximal ones.

**Step 6 (enumeration is complete).** In A_p, depth-first search from 0; the children of an
isotropic S are S + Ox for every x ∈ A_p ∖ S with b_A(x, x) = 0 and b_A(x, g) = 0 for the
O-generators g of S (then b_A(x, S) = 0 by O-bilinearity, and S + Ox is isotropic since
b_A(ox + s, o′x + s′) = oo′ b_A(x, x) + o b_A(x, s′) + o′ b_A(s, x) + b_A(s, s′) = 0).  If T ⊋ S is
isotropic, any x ∈ T ∖ S is a child with S + Ox ⊆ T; by induction on |T/S| every isotropic T ⊇ S
is visited.  Nodes without children are exactly the maximal isotropic submodules.  So the search
lists all isotropic submodules (needed for Max_free) and all maximal ones (Max).

**Step 7 (symmetries; only used to skip frames).** (α, β, b) and (α, β, −b) give the same lattice
(w ↦ −w); (β, α, b) an isometric one.  Galois: x ↦ x^σ (coordinatewise) maps L_f onto L_{f^σ}
(form G_f^σ), integral/maximal/free overlattices onto integral/maximal/free overlattices (the
Steinitz ideal goes to its conjugate, principal ↔ principal), and Q(x^σ) = Q(x)^σ; σ preserves
I(K) and square classes, so s□ is preserved.  f^σ is brought back into F by Step 1's unit
rescaling.  One frame per orbit of ⟨±b, swap, σ⟩ is processed.

**Step 8 (exact evaluation).** s□(M) is computed by L6.1 (PROVED given Kala–Tinková Prop 3.1):
s□(M) = #classes of ∪_E {Q(v) : v ∈ M, Tr(δ_E Q(v)) = 1}, E over the s edge orbits modulo (O^×)².

**Step 9 (freeness).** For M with pseudo-basis 𝔞1 e1 ⊕ 𝔞2 e2 the Z-span 𝔡(M) of {det(x, y) : x, y ∈ M}
is the O-ideal 𝔞1𝔞2 det(e1, e2); M is free ⇔ 𝔡(M) is principal.  An ideal g1(Z + Zθ) is principal
⇔ Z + Zθ is homothetic to O = Z + Zω ⇔ θ ~ ω under GL₂(Z) ⇔ (Serret) the continued fractions of θ
and ω have a common complete quotient in their periods.  ∎

## 2. Implementation (src/complete.py) — what is exact
- `FieldData`: R = balanced reps (minimal trace in the (O^×)²-orbit) of `indec.indecomposables`
  (± seeds αη when N(ε) = +1); exact normalization with exponent tracking.
- `admissible_b`: b with αβ − b² ≻ 0 satisfy Tr(b²/(αβ)) < 2, a positive definite **rational**
  binary form in the coordinates of b; its ellipse is enumerated exactly (Lagrange reduction,
  integer square roots), then αβ − b² ≻ 0 is tested exactly.  Float-free.  (A float box search
  was hopeless for huge ε, e.g. D = 958, ε ≈ 3.4·10^16.)
- `DiscModule`: A = O²/G O² by an integer HNF (|A| = N(det) asserted), pairing via
  det′·adj(G)·y mod N(det) (integers), p-primary parts, DFS of Step 6.
- `overlattice_basis`: Z-basis of M = O² + Σ O G^{-1}y_g by HNF.
- `lattice.exact_sq_classes` (L6.1): exact LLL + Fincke–Pohst with cap 1 per edge; `edge_deltas`
  asserts δ_E ≻ 0, δ_E ∈ O^∨ for every D.
- `steinitz_ideal`, `is_principal` (Serret, exact).
- Floats: none in decisions.  (`emb_stable`, `is_principal_box` are diagnostics/references.)

## 3. Validation (all in `tests/test_complete.py` unless stated)
1. `admissible_b` == independent float-box enumerator (`gram.valid_off_diagonal` + exact filter)
   on all rep pairs, D ∈ {2, 6, 13, 19, 33, 57} (also 67, 94 by hand); ellipse enumerator == brute
   force on 300 random forms.
2. Maximal integral overlattices, and ALL integral overlattices, == brute force (BFS over Z-bases
   in Q⁴, coset reps of L^#/L from an independent rational solve, integrality tested directly) on
   random frames with |A| ≤ 40, D ∈ {2, 3, 6, 10, 13, 15, 19, 26, 30, 33, 57}.
3. Galois reduction harmless (S with/without, D = 13, 19, 33).
4. Exact (per-edge) s□ == capped Fincke–Pohst class count for maximal overlattices with s□ ≥ 4
   (D = 19, 33; slow).
5. Principal-ideal test == norm-box search (8 fields) and == PARI `bnfisprincipal` (12 fields,
   960 random ideals in a separate run, 211 non-principal); skipped without cypari2.
6. S_free == brute force over all free overlattices (D = 10, 15, 26, 30).
7. **sails cross-check** (`scripts/check_sails.py`, independent PARI-based code
   github.com/pav10/sails): for **all 1823 squarefree D ≤ 3000** the indecomposables mod O^{×,+},
   the facet functionals (= `edge_deltas`), integer distance d = 1 of every facet (= KT Prop 3.1)
   and #facets = s/2 or s agree.  So Step 8's imported input is re-verified for D ≤ 3000.
8. Reproduces the IMPORTED complete values D = 26 → 2, 33 → 4, 19 → 6, and the CERTIFIED lower
   bounds 43 ≥ 8, 19/31/57 ≥ 6.

## 4. Results, all squarefree D ≤ 100 (CERTIFIED; `data/complete_S.jsonl`)
`python3 scripts/complete_S.py out.jsonl --range 2 100` (≈ 6 min on 3 cores).

| D | s | κ□ | S | S_free | | D | s | κ□ | S | S_free | | D | s | κ□ | S | S_free |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 2 | 1 | 2 | 2 | 2 | | 35 | 2 | 2 | 2 | 2 | | 66 | 2 | 14 | **4** | 2 |
| 3 | 2 | 2 | 2 | 2 | | 37 | 3 | 5 | 2 | 2 | | 67 | 10 | 34 | **8** | 8 |
| 5 | 1 | 1 | 1 | 1 | | 38 | 2 | 12 | 4 | 4 | | 69 | 4 | 4 | 2 | 2 |
| 6 | 2 | 4 | 2 | 2 | | 39 | 2 | 8 | 2 | 2 | | 70 | 6 | 10 | 3 | 3 |
| 7 | 4 | 4 | 2 | 2 | | 41 | 5 | 7 | 4 | 4 | | 71 | 8 | 12 | 4 | 4 |
| 10 | 1 | 5 | 2 | 2 | | 42 | 2 | 4 | 2 | 2 | | 73 | 9 | 9 | 4 | 4 |
| 11 | 2 | 6 | 2 | 2 | | 43 | 10 | 22 | **8** | 8 | | 74 | 5 | 17 | 3 | 3 |
| 13 | 1 | 3 | 2 | 2 | | 46 | 12 | 12 | 4 | 4 | | 77 | 2 | 2 | 2 | 2 |
| 14 | 4 | 4 | 2 | 2 | | 47 | 4 | 4 | 2 | 2 | | 78 | 4 | 4 | 2 | 2 |
| 15 | 2 | 2 | 2 | 2 | | 51 | 2 | 12 | **4** | 2 | | 79 | 4 | 4 | 2 | 2 |
| 17 | 3 | 3 | 2 | 2 | | 53 | 1 | 7 | 2 | 2 | | 82 | 1 | 17 | **4** | 2 |
| 19 | 6 | 10 | 6 | 6 | | 55 | 4 | 6 | 2 | 2 | | 83 | 2 | 18 | 4 | 4 |
| 21 | 2 | 2 | 2 | 2 | | 57 | 6 | 10 | 6 | 6 | | 85 | 1 | 8 | **3** | 2 |
| 22 | 6 | 8 | 4 | 4 | | 58 | 7 | 17 | **6** | 4 | | 86 | 10 | 28 | **8** | 8 |
| 23 | 4 | 4 | 2 | 2 | | 59 | 6 | 14 | 4 | 4 | | 87 | 2 | 6 | 2 | 2 |
| 26 | 1 | 9 | 2 | 2 | | 61 | 3 | 9 | 4 | 4 | | 89 | 7 | 13 | 5 | 5 |
| 29 | 1 | 5 | 2 | 2 | | 62 | 4 | 4 | 2 | 2 | | 91 | 8 | 18 | **6** | 4 |
| 30 | 2 | 4 | 2 | 2 | | 65 | 3 | 6 | 2 | 2 | | 93 | 2 | 6 | 4 | 4 |
| 31 | 8 | 12 | 6 | 6 | | | | | | | | 94 | 16 | 16 | 6 | 6 |
| 33 | 4 | 4 | 4 | 4 | | | | | | | | 95 | 4 | 4 | 2 | 2 |
| 34 | 4 | 4 | 2 | 2 | | | | | | | | 97 | 9 | 13 | 6 | 6 |

Long-period fields (D ≤ 1000, s = 36–44): §6.

## 5. Findings (statuses as entered in 02_ledger.md, L7)
- **L7.2 CERTIFIED.** max_{D ≤ 100} S(K,2) = 8, attained exactly at D = 43, 67, 86 (all s = 10).
  S(Q(√43),2) = 8 exactly (the L6.2 lower bound is the maximum).
- **L7.3 CERTIFIED.** Non-free lattices beat binary forms: S_free < S for D = 51, 58, 66, 82,
  85, 91 (h(K) = 2, 2, 2, 4, 2, 2).  The invariant is about lattices, and "forms" give a
  different (smaller) number when h(K) > 1.
- **L7.4 CERTIFIED (REFUTES the R3 working hypothesis "≤ 2 classes per edge" and L6.3's "≤ 2
  values per edge when s□ ≥ 4").** D = 82 (s = 1, one edge orbit): a maximal lattice represents
  4 classes, all on the single edge orbit; S(Q(√82),2) = 4 > 2s.  D = 85 (s = 1): S = 3.  D = 58:
  an s□ = 6 lattice with 4 values on one edge.
  Explicit witness (Q(√82), h = 4): Gram [[10−√82, −4], [−4, 10+√82]] on O e1 ⊕ O e2, glued with
  x = (7 + √82/2, 7 − √82/2) (Q(x) = 14): M has Z-basis (1,0), (√82/2, √82/2), (0,1), (0,√82); M is
  NOT free.  Its values on the unique edge orbit: 10−√82, 154−17√82 (norm 18), 46−5√82,
  118−13√82 (norm 66) — four distinct square classes.  Over all maximal lattices (D ≤ 100) the classes on one edge reach 3 (D = 73, 74, 85, 89, 97) and
  4 (D = 58, 82), never more.  Per-edge statistics over every maximal lattice:
  `data/edge_stats.jsonl`.
- **L7.5 EMPIRICAL.** S(K,2) is not monotone in s, κ□ or ι: s = 2 fields reach 4 (D = 38, 51,
  66, 83, 93), s = 12 (D = 46) gives 4, s = 16 (D = 94) gives 6.  S ≥ s happens only for s ≤ 6.
- **L7.6 EMPIRICAL.** In 20 of the 25 fields with S ≥ 4 the first extremal frame found is
  Galois-self-conjugate, [[α, b], [b, α′]], b ∈ Z (all three S = 8 fields; det = N(α) − b² ∈
  {2, 8, 8}).  Supports R2 (E = K(√−d) biquadratic).
- **L7.7 CERTIFIED.** Restricted to Galois-symmetric frames [[α, b], [b, α′]] (b ∈ Z):
  S_sym = S in 21 of 25 fields with S ≥ 4; S_sym < S for D = 31, 89, 94, 97
  (data/symmetric_S.jsonl).
- Hand-over data for D ≡ 1 (mod 4) realized lattices were computed on Z[√D]-sublattices
  (03_corrections C17).

## 6. Long-period fields
(see the table appended below when the runs finish)
