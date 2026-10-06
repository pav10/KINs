# 02 — Results ledger (status-tagged)

Tags: PROVED (proof here) · IMPORTED (source given; "re-verified" if we checked) ·
CERTIFIED (exact computation, command given) · EMPIRICAL · HEURISTIC (gap named) · OPEN · REFUTED.
Notation as in 01_problem.md. Throughout M is integral: Q(M) ⊆ O, B(M,M) ⊆ O unless noted.

---------------------------------------------------------------------------------------------
## L0 Foundations

**L0.1 PROVED.** N_r(Q) = 1: the only indecomposable of Z_{>0} is 1 (n = 1+(n−1)).
The phenomenon is a degree ≥ 2 (units / partial order) effect.

**L0.2 PROVED (value-minimality; the primitive).** If α = Q(v) ∈ I(K), there is no
y ∈ M∖0 with Q(y) ⪯ α, Q(y) ≠ α.
*Proof.* Q(y) ≻ 0. α − Q(y) is a nonzero totally nonnegative element of K, hence totally
positive (a nonzero element of K has no zero conjugate), and lies in O. So
α = Q(y) + (α − Q(y)) is a decomposition. ∎
Also: distinct indecomposables are ⪯-incomparable.

**L0.3 PROVED (acute irreducibility).** If Q(v) ∈ I(K) then v ∉ {y+z : y,z ≠ 0, B(y,z) ⪰ 0}.
*Proof.* Q(v) = Q(y) + (Q(z) + 2B(y,z)) with both summands in O⁺. ∎
(So representing vectors lie in O'Meara's A(M); L0.2 is strictly stronger.)

**L0.4 PROVED (totally positive units are indecomposable).** If ε ∈ O^{×,+} and ε = β+γ,
β,γ ∈ O⁺, then βε⁻¹ + γε⁻¹ = 1 with βε⁻¹, γε⁻¹ ∈ O⁺; each has all conjugates in (0,1), so
norm in (0,1), impossible for a nonzero algebraic integer. ∎
(Supersedes a WRONG retraction issued during the project — see 03_corrections C9.)

**L0.5 PROVED (diagonal forms).** If M = ⊥ O v_i (or any orthogonal sum of rank-1 lattices
𝔞_i v_i), then Q(x) ∈ I(K) forces x in a single summand; hence s□(M) ≤ r.
*Proof.* Q(Σx_i) = ΣQ(x_i), a sum of ≥ 2 elements of O⁺ if two components are nonzero;
on one summand Q(ξ v_i) ∈ Q(v_i)(K^×)². ∎ Consequence: sums of squares are useless
witnesses; all content is in off-diagonal (non-split) lattices.

**L0.6 PROVED (per-field finiteness).** s□(M) ≤ κ□(K) < ∞ (κ□ finite: I(K) is finite modulo
O^{×,+}). Only uniformity is at issue.

**L0.7 PROVED (realizability for h = 1).** If an abstract Gram matrix G = (g_ij) over O of
rank ≤ 2, totally PSD, is realized by v_i ∈ K² (G = Gram(v_i)), then L = Σ O v_i is an
integral O-lattice of rank 2 whose values include all g_ii. If h(K) = 1, L is free.
(No pivot-integrality is needed; see pitfall P2.)

---------------------------------------------------------------------------------------------
## L1 Uniform partial bounds

**L1.1 PROVED (representation numbers; O'Meara 3.7 transplanted).** For α ∈ I(K):
r(α,M) := #{v : Q(v) = α} ≤ 2(2^{rd} − 1).
*Proof.* No representing v lies in 2M (Q(2u) = 4Q(u) = Q(u) + 3Q(u)). Let v, w represent α
with w ≡ v mod 2M, w = v − 2u, u ∉ {0, v}. Then Q(w) − Q(v) = −4B(v−u, u). Since v is
acutely irreducible (L0.3) and v = (v−u) + u with both parts nonzero, B(v−u,u) ⋡ 0, so
σ(B(v−u,u)) < 0 at some σ, giving σ(Q(w)) > σ(Q(v)) — contradicting Q(w) = Q(v). Hence each
nonzero coset of M/2M (2^{rd} − 1 of them) holds at most the pair ±v. ∎
**Scope:** bounds one value; does NOT bound s□ (L3.2).

**L1.2 PROVED (pincer).** v, w independent, α = Q(v), β = Q(w), b = B(v,w). If
N(α)N(β) ≤ N(𝔞)² then b = 0.
*Proof.* Strict Cauchy–Schwarz at each place on the definite plane: σ(b)² < σ(α)σ(β).
Multiply: N(b)² < N(α)N(β) ≤ N(𝔞)². If b ≠ 0, (b) ⊆ 𝔞 gives |N(b)| ≥ N(𝔞). ∎
**Corollary (norm window).** #{α(K^×)² : α ∈ Q(M)∩I(K), N(α) ≤ N(𝔞)} ≤ r, uniformly in K, d.
(Distinct classes ⇒ non-collinear representatives ⇒ pairwise orthogonal ⇒ ≤ r lines;
one class per line.)

**L1.3 PROVED (units).** At most r classes mod (K^×)² of totally positive units are
represented by an integral rank-r M. (L0.4 + L1.2 with N(α)N(β) = 1 ≤ N(𝔞)², 𝔞 ⊆ O.)
IMPORTED source of the regime: arXiv 2409.11082 (recovered mechanism).

**L1.4 PROVED (per-face lemma).** Let δ ∈ O^∨ (codifferent), δ ≻ 0, and
T_δ(v) := Tr_{K/Q}(δ Q(v)) on M (a Z-lattice of rank rd). T_δ is Z-valued (δ ∈ O^∨,
Q(v) ∈ O) and positive definite (δQ totally positive definite), so T_δ ≥ 1 on M∖0.
If Q(v) ∈ I(K) with Tr(δQ(v)) = 1 then v is a minimal vector of T_δ. Hence
#{α ∈ I(K) ∩ Q(M) : Tr(δα) = 1} ≤ ½ τ_{rd}  (τ_n = max kissing number of n-dim lattices;
τ_4 = 24, so ≤ 12 for binary/real quadratic). ∎

**L1.5 PROVED modulo IMPORTED structure (sail bound, real quadratic).** For K = Q(√D):
  s□(M) ≤ ½ τ_{2r} · s       (binary: s□(M) ≤ 12 s),  s = CF period of ω.
*Proof.* IMPORTED (Kala–Tinková 2005.12312 Prop 3.1; Blomer–Kala 1705.03671): indecomposables
lie on the edges of the sail; each edge E has δ_E ≻ 0 in O^∨ with Tr(δ_E α) = 1 for all α ∈ E;
edges modulo O^{×,+} number f = s/2 (s even) or s (s odd) — matching the verified count
ι = Σ_{edges} u_i (tests/test_rqf_indec.py). Since Q(εv) = ε²Q(v), every represented class has
a represented indecomposable on an edge from a fixed set of representatives modulo (O^×)²;
there are f·[O^{×,+} : (O^×)²] = s such edges (N(ε) = (−1)^s: index 2 iff s even). Apply L1.4
per edge. ∎
Note: independent of the partial quotients u_i, hence of ι = #indecomposables.
**Corollary.** Bounded period ⇒ bounded s□ (e.g. s = 2 ⇒ s□ ≤ 24).

**L1.6 IMPORTED, NOT re-verified (regulator ceiling; other analysis "Thm 17").**
s□(M) ≤ V_2^d 2^{d−1} Reg(K) / √N(𝔳M) + E (volume count of vectors with N(Q(u)) ≤ Δ modulo
units). Consequence claimed there ("Cor 21.3"): the proportion of indecomposable classes γ with
αγ a norm from E = K(√−ε⁺) tends to 0 along D = t²+1. Verify before use (task R4).

---------------------------------------------------------------------------------------------
## L2 Exact reformulations (binary)

**L2.1 PROVED (frame criterion).** Let v, w ∈ M, α = Q(v) ≠ 0, b = B(v,w), β = Q(w),
δ = αβ − b² (≻ 0). For γ ∈ K:
  γ ∈ Q(Ov + Ow) ⟺ ∃ X, Y ∈ O: αγ = X² + δY² and X ≡ bY (mod αO).
*Proof.* αQ(xv+yw) = (αx + by)² + δy²; put X = αx+by, Y = y; conversely x = (X − bY)/α ∈ O
iff the congruence holds. ∎ (Applies to the sublattice Ov+Ow; M may be larger.)

**L2.2 IMPORTED (standard; Gauss/Dedekind/Clifford composition; Shyr 1977).** Binary M with
det class δ ↔ an ideal class [𝔞_M] of an order 𝒪 ⊆ E = K(√−δ) (CM, E/K cyclic); Q(M) = values
N_{E/K}(z)/N(𝔞_M), z ∈ 𝔞_M. Isometry classes in the genus of a free lattice ↔ classes counted by
Shyr's H = (h_E/h_K)·2^{1−t}/[U_E^0 : U_K^0] (t = # primes of K ramified in E). Care needed with
non-maximal orders, proper vs improper equivalence and Steinitz class (h(K) > 1).

**L2.3 PROVED given L2.2 (one fibre).** With Φ : I(K) → Cl(𝒪), α ↦ class of an ideal of
relative norm (α) (where solvable),  s□(M) = #{classes [α] with Φ-fibre meeting [𝔞_M]}:
the represented indecomposables of a binary lattice form one fibre of Φ.
⇒ **Uniform bound ⟺ sup over (K, δ, class) of fibre sizes (mod squares) is finite.**

**L2.4 PROVED (unit determinant).** If δ ∈ O^{×,+}, disc(O[√−δ]/O) = (−4δ) = (4), so E/K is
unramified outside 2·∞; local norm conditions away from 2 are automatic/congruence conditions
mod a power of 2. **But** the ideal-class condition lives in Cl(E), with h⁻(E) = Δ_K^{1/2+o(1)}
(Brauer–Siegel, since Δ_E/Δ_K = Δ_K·N(𝔡_{E/K}) and N(𝔡_{E/K}) | 16). So the class-placement
problem is NOT reduced to a finite 2-adic group (see corrections C10).

**L2.5 EMPIRICAL (IMPORTED data, re-verified for Δ=57).** In every certified extremal example
(Δ = 33, 57, 76) the frame determinant is exactly δ = ε⁺ (the totally positive fundamental
unit). Re-verified: Δ=57, α1 = 10+3ω, α2 = 23+7ω, b = −(13+4ω): δ = 151+20√57 = η. Also
N(α1) = N(α2) in all three (4,4 / 4,4 / 9,9). No mechanism known.

---------------------------------------------------------------------------------------------
## L3 Closed routes

**L3.1 REFUTED — linear-irreducibility route. CERTIFIED.**
Hypothesis: every represented indecomposable class has a representative that is Z-linearly
irreducible in the trace lattice (M, Tr∘Q) (v not in the Z-span of trace-shorter vectors).
It would give s□ ≤ 2(2^{rd}−1) via O'Meara 3.18 (linear irreducibility descends to Z).
Counterexamples (h = 1; realized module L = ΣOv_i; exact HNF; `pytest -m slow`,
`scripts/experiments.py analyze 19:600 31:1400`, `frame57`):
- Q(√31): s□ ≥ 6 (value traces 12,34,56,490,902,1314), only [6+√31] linearly irreducible; n* = 20.
- Q(√19): 4 classes, 2 linearly irreducible; n* = 20.
- Q(√57), frame of L2.5: 6 classes, **none** linearly irreducible; n* = 11.
(n* = trace at which shorter vectors already Z-span the lattice.)

**L3.2 PROVED (mod-2M route cannot count classes).** A(M) is closed under multiplication by
units, hence infinite for d ≥ 2 (O'Meara 3.17 analogue fails); some coset of M/2M contains
infinitely many acutely irreducible vectors. L1.1 survives only for fixed values.

**L3.3 REFUTED — metric / spherical-code routes.**
(a) PROVED: the spectral route needs root discriminant < √2 — impossible (min is √5).
(b) IMPORTED (other analysis "Prop 5"), not re-verified: for every δ > 0 there are real
quadratic K, binary M, v, w representing indecomposables in distinct classes with trace-angle
cosine > 1−δ. Separation exists only at a pair-dependent place. Certified extremal vectors
have trace-angle cosines ≈ −0.98 (consistent).

**L3.4 IMPORTED (Ramsey reduction; other analysis "Thm 4").** If s□(M) ≥ R_{d(d−1)}(k), there
are k represented indecomposables forming a two-place staircase. Reduction only.

**L3.5 HEURISTIC (Pareto / value-sail counts are Θ(s)).** Pareto-minimal vectors of the
conjugate pencil, modulo units, ↔ ∪_λ Min(Q_λ) over a torus of covolume Reg; breakpoints
hyperbolically O(1)-spaced ⇒ Θ(s) per period. Gap: the spacing claim is not proved. Moral
(robust): geometry-of-numbers counts of extremal value-vectors ignore the constraint that values
land on the FIELD sail; they cannot beat Θ(s).

**L3.6 CERTIFIED (growth is not governed by D or ι).** Period-2 fields, fan 1 + rα0
(α0 = u0+√D), traces in arithmetic progression, one shared δ: best found lattice represents
exactly 2 classes for ι = 8 (D=66) and ι = 12 (D=146). Consistent with L1.5.

---------------------------------------------------------------------------------------------
## L4 Data
See 07_data.md and 09_R1_complete.md. Headline (R1, CERTIFIED, complete): S(K,2) exactly for all
squarefree D ≤ 100; max 8, attained at D = 43, 67, 86 (s = 10). The IMPORTED complete values
D = 26 (S=2), 33 (S=4=κ□), 19 (S=6) are re-verified. The claim "S ≤ 6 for 40 fields up to
Δ = 348" (restricted to unit-determinant frames) is REFUTED (L6.2; also D = 67, 86 with
Δ = 268, 344).

---------------------------------------------------------------------------------------------
## L5 Open

**L5.1 OPEN (rank-1 / unit-orbit core; IMPORTED "A7").** T(K) = sup_γ #{ξO : γξ² ∈ I(K)}.
Each line contributes 1 class mod (K^×)² but up to T(K) mod (O^×)². Prescribed-CF (Friesen-type)
constructions suggest sup_K T(K) = ∞, which would make s_u non-uniform already in rank 1.
EMPIRICAL support that square classes contain several unit orbits: κ□ < 2ι for D = 19 (10 < 14),
D = 31 (12 < 16) with N(ε) = +1.

**L5.2 OPEN — central dichotomy (reformulated after L6.3; sub-question (i) settled by L7.4).**
Is S(K,2) bounded, or ≍ s? (Exact values: L7.2, 09_R1_complete.md.)
The earlier "faces-touched lemma" (one lattice meets O(1) edges) is contradicted by data:
extremal lattices meet all or almost all edges of a period (L6.3). What stays small is the number
of values per edge (≤ 2 in all data) and the collapse of values into classes.
Sub-questions: (i) per-edge CLASS count ≤ 2 (would give s□ ≤ 2s) — REFUTED, L7.4 (D = 82: 4 classes
on one edge orbit, S = 4 > 2s); (ii) collapse mechanism (L6.4);
(iii) a family with s□ → ∞. Side data: the odd convergents (sail VERTICES) of one period are not
jointly representable (max 2/3, 2/4, 3/5, 2/6 for D = 19, 31, 43, 46; lower-bound searches) —
represented indecomposables are mostly interior edge points, not vertices.

**L5.3 OPEN — bounded fibres of Φ (L2.3)** in the range N(α) ≤ Δ/4 ≪ h(E)^{1+ε}: below
effective Duke / Linnik / Duke–Schulze-Pillot; generic analytic input stops at Siegel.
HEURISTIC caution: κ□(K) can be ≍ Δ^{1/2+o(1)} and #classes in the genus ≍ Δ^{1/2+o(1)};
a *random* placement would give max fibre ≍ log H/log log H → ∞ slowly. So a uniform bound
needs STRUCTURE (sail geometry, L1.4), not just equidistribution. Data (max 6) suggests strong
structure.

**L5.4 OPEN.** Which invariant separates Q(√26) (S = 2) from Q(√19) (S = 6)? Same h = 1,
h⁺ = 2, 2 ramified. (IMPORTED data.)

**L5.5 OPEN.** Genus average: Σ_{M' ∈ gen(M)} s□(M') ≤ κ□(K)·max_α #{ideals of E of norm (α)}
(each class may be hit by several genus classes) — divisor-type factor Δ^{o(1)}. Make precise
with Shyr's H (L2.2) and compare with computed maxima.

**L5.6 OPEN (higher degree).** t_d := sup of min_{δ ∈ O^{∨,+}} Tr(δα) over indecomposables of
degree-d fields (t_2 = 1; t_3 ≤ 3 IMPORTED Kala–Tinková). If t_d < ∞ then L1.4 generalizes with
½τ_{rd} replaced by #{v : T_δ(v) ≤ t_d} — but that count is not uniform in the lattice
(it depends on the successive minima of T_δ); needs work.

---------------------------------------------------------------------------------------------
## L6 Findings of the hand-over session (exact method)

**L6.1 PROVED (exact, cap-free evaluation of s□(M), real quadratic).** For each sail edge E (reps
modulo (O^×)², s of them) with functional δ_E ∈ O^{∨,+}, Tr(δ_E·E) = 1:
  {indecomposables represented by M} ∩ (edge E) = {Q(v) : v ∈ M, T_{δ_E}(v) = 1}.
⊆ by L1.4. ⊇: if T_δ(v) = 1 and Q(v) = β+γ with β, γ ∈ O⁺ then Tr(δβ), Tr(δγ) ∈ Z_{≥1}
(δ ∈ O^{∨,+}), sum ≥ 2 — contradiction; so Q(v) is indecomposable (and lies on E).
Hence s□(M) = #classes of ∪_E {Q(v) : T_{δ_E}(v) = 1}: finitely many norm-1 vectors of s
positive-definite rank-4 Z-lattices. Implemented: `lattice.edge_deltas` (computes δ_E from
Tr(δα_i)=1, Tr(δα_{i+1})=0 and ASSERTS δ_E ≻ 0, δ_E ∈ O^∨ — a computational check of the
imported KT Prop 3.1 for every D run), `lattice.exact_sq_classes` (exact LLL + Fincke–Pohst
with cap 1). Agrees with capped enumeration where the cap suffices; exposes cap failures (D = 43:
cap 3000 finds 6 classes, exact finds 8).

**L6.2 CERTIFIED (S(Q(√43),2) ≥ 8; = 8 exactly by R1, L7.2).** Winner on balanced class reps (`find_winner`), realized
L = ΣOv_i in the ambient form Gp = [[7+√43, −2], [−2, 7−√43]] (det Gp = 2). L is integral,
positive definite, Z-trace-Gram determinant 29584 = 172², i.e. N(𝔳L) = 1. Represented classes
(exact): 7−√43, 1541−235√43, 47207−7199√43, 10731517−1636541√43 (norm 6);
105−16√43, 282−43√43, 730938−111467√43, 1963743−299468√43 (norm 17); pairwise inequivalent mod
(K^×)². Values touch 10/10 edges (per-edge counts 1,1,2,1,1,1,1,2,1,1).
Command: `python3 scripts/scan_exact.py 43 43 600 /tmp/o.jsonl` or see tests (slow).

**L6.3 CERTIFIED (edges touched; per-edge counts).** Exact per-edge value counts:
D = 19 (orbit-rep winner): [1,2,1,1,2,1] → 6 classes, 6/6 edges; D = 31 (balanced):
[1,2,1,0,1,2,1,0] → 6 classes; D = 43: 10/10 edges, 8 classes; D = 57 frame (L2.5):
[2,2,2,2,2,2] → 6 classes, 6/6 edges; D = 46 best found: 8/12 edges, 4 classes.
EMPIRICAL claim "≤ 2 values per edge in every lattice with s□ ≥ 4" — REFUTED by L7.4 (D = 82, D = 58).
Diagonal lattices ⟨1, γ⟩ in odd-period fields carry up to 4 values on one edge (all in ≤ 2
classes): D = 2, 10, 26 (s = 1) give 3–4, D = 17, 37, 41 give 3, D = 65 gives 4 (recomputed after
C17; the stored data/scan_exact.jsonl rows for D ≡ 1 mod 4 are on Z[√D]-sublattices).
These extra values are indecomposable SQUARES ξ² (ξ non-unit) — the rank-1 core L5.1 in action.

**L6.4 EMPIRICAL (shape of extremal lattices).** Two patterns: (i) det = η (tp fundamental
unit): D = 19 winner (det 170+39√19), imported Δ = 33, 57, 76; (ii) Galois-self-conjugate
ambient frames [[α, b], [b, α']], b ∈ Z, det = N(α) − b² = 2 (D = 22, D = 43). In D = 57 the 12
values form 6 classes: values from different edges coincide mod (K^×)² (ratio ξ², ξ ∉ O^×,
equal norms) — the collapse mechanism behind small s□ is "value pairs across edges", not
"few edges".

**L6.5 PROVED (no realizability gap for lattices).** If G = (g_ij) ∈ M_k(O) is totally PSD of
rank 2 with a positive-definite principal 2×2 block, then v_i := coordinates of row i w.r.t. that
block realize G in (K², block form), and L = ΣOv_i is an integral positive-definite binary
O-lattice with Q(v_i) = g_ii. (Rank 2 ⇒ g_ij = g_{i,P} G_P^{-1} g_{P,j}.) Freeness (h(K) = 1) is only
needed if one insists on free lattices / binary FORMS.

---------------------------------------------------------------------------------------------
## L7 R1 — complete computation (docs/09_R1_complete.md, src/complete.py)

**L7.1 PROVED (R1 theorem).** S(K,2) = max over frames f = (α, β, b) (α, β indecomposables mod
(O^×)² in distinct classes, b ∈ O, αβ − b² ≻ 0) of max s□(M) over the maximal integral overlattices
M of L_f; S_free(K,2) (free lattices = binary forms) likewise over lattices maximal among FREE
integral overlattices. Overlattices = isotropic O-submodules of O²/G_fO² (p-primary parts
orthogonal). Proof: 09 §1. Uses L6.1 (hence KT Prop 3.1, IMPORTED; re-verified computationally
for all squarefree D ≤ 3000 against the independent `sails` package, 09 §3.7).

**L7.2 CERTIFIED (exact values).** S(K,2) and S_free(K,2) for all 60 squarefree D ≤ 100 (09 §4,
data/complete_S.jsonl; `python3 scripts/complete_S.py out.jsonl --range 2 100`). max = 8, exactly at
D = 43, 67, 86 (s = 10). Re-verifies the IMPORTED D = 26 → 2, 33 → 4, 19 → 6.

**L7.3 CERTIFIED (lattices ≠ forms).** S_free < S for D = 51 (2 < 4), 58 (4 < 6), 66 (2 < 4),
82 (2 < 4), 85 (2 < 3), 91 (4 < 6); h(K) = 2, 2, 2, 4, 2, 2. For h(K) = 1 they coincide.

**L7.4 CERTIFIED (REFUTES "≤ 2 classes per edge", R3, and L6.3's "≤ 2 values per edge if s□ ≥ 4").**
Q(√82) (s = 1): a maximal integral, non-free lattice (frame [[10−√82, −4], [−4, 10+√82]] glued with
x = (7+√82/2, 7−√82/2); values 10−√82, 154−17√82, 46−5√82, 118−13√82) represents 4 classes, all on
the single edge orbit; S = 4 > 2s. Q(√85) (s = 1): S = 3. Q(√58): an
s□ = 6 lattice with 4 values on one edge. Also for binary FORMS: 3 classes on one edge occur for D = 73, 89, 97 (h = 1, all lattices free).
Hence s□ ≤ 2s is false; the per-face bound L1.4 (≤ 12 per
edge) is the only per-edge bound known. Over every maximal lattice, D ≤ 100: classes on one edge reach 3 (D = 73, 74, 85, 89, 97)
and 4 (D = 58, 82); never more. Per-edge maxima over every maximal lattice:
data/edge_stats.jsonl (`scripts/edge_stats.py`).

**L7.5 EMPIRICAL.** S(K,2) is not monotone in s (s = 2 fields reach 4; s = 12: 4; s = 16: 6) and
S ≥ s occurs only for s ≤ 6 (D ≤ 100). Long periods: 09 §6.

**L7.6 EMPIRICAL (shape of extremal lattices).** In 20 of the 25 fields D ≤ 100 with S ≥ 4 the first
extremal frame found is Galois-self-conjugate [[α, b], [b, α′]] with b ∈ Z (including D = 43, 67, 86,
det = 2, 8, 8); among these first extremal frames, unit determinant (L2.5) occurs only for
D = 33, 41, 61, 71. Supports R2.

**L7.7 CERTIFIED (Galois-symmetric frames, exact over that class).** S_sym(K) := max s□ over the
maximal integral overlattices of frames [[α, b], [b, α′]], α ∈ I(K), b ∈ Z (complete.symmetric_S;
data/symmetric_S.jsonl). D ≤ 100: S_sym = S in 21 of the 25 fields with S ≥ 4 (all of D = 43, 67,
86, 82, 58, 91, 19, 57); S_sym < S for D = 31 (4 < 6), 89 (4 < 5), 94 (4 < 6), 97 (4 < 6). So the
symmetric case (R2) is the typical extremal shape but not the only one.

**L7.8 CERTIFIED (lower bounds well above 8).** The det-2 Galois-symmetric frames
G_{a,b} = [[a − √D, b], [b, a + √D]], D = a² − b² − 2 (the extremal shape for D = 19, 22, 38, 43, 46, 58,
73, 82, 862), glued to their maximal integral overlattices, give
  S(Q(√3931),2) ≥ 24 (a = 63, b = 6, s = 130, free),
  S(Q(√691),2) ≥ 20 (a = 27, b = 6, s = 38),  S(Q(√823),2) ≥ 20 (a = 35, b = 20, s = 44),
  also 20 for D = 1303 (s = 52), 1579 (s = 70); 16 for D = 739; 14 for D = 331; 12 for D = 862.
The D = 691, 823 witnesses are FREE (binary forms). Certificates: data/certificates/*.json
(vectors, values; `scripts/certify_lattice.py D a b`; re-verified by test_stored_certificates and,
independently, indecomposability by `sails` and square classes by PARI nfroots).
Scan: `scripts/family_det2.py out.jsonl 5000` (data/family_det2.jsonl). For D ≡ 2, 3 (mod 4) 2 ramifies
and the index-2 glue makes M unimodular (𝔳M = O) — the unit-determinant pattern L2.5.
Family records (D ≤ 20000, 5354 lattices in 3766 fields): 6, 8, 10, 14, 20, 24, 28 at D = 19, 43, 271,
331, 691, 3931, 5419; **S(Q(√5419),2) ≥ 28** (a = 89, b = 50, s = 122, free; certificate
D5419_det2.json). After D = 5419 no family lattice exceeds 28 (28 recurs at D = 12919, 17431, 17491;
max per D-block of 2500: 20, 24, 28, 22, 24, 28, 28, 20; per s-block of 50: 20, 20, 28, 28, 28, 12, 20)
— EMPIRICAL plateau of the det-2 family at 28 for D ≤ 20000. Says nothing about S(K,2) beyond the family. Complete runs (R1): S = 8 for D = 478, 958 (s = 36), 718 (s = 40); S = 12 for D = 862 (s = 40),
attained by the family lattice — the family gives the exact maximum there.
Consequence for L5.2: the D ≤ 100 maximum 8 is not a ceiling; the "lean toward boundedness" (C7)
is not supported. Within the family the growth in s is irregular (several s ≥ 40 fields give ≤ 8).
