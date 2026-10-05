# CLAUDE.md — Indecomposables represented by binary quadratic forms over real quadratic fields

## Mission
Study which (square classes of) indecomposable totally positive integers a single
totally positive definite quadratic lattice can represent.
**Step 1 (this repo): binary lattices over real quadratic fields K = Q(√D).**
Target theorem: decide whether S(K,2) := max_M s□(M) is bounded by an absolute
constant, and if so prove it; if not, identify the governing invariant and growth rate.
Long-term: rank r, degree d — N_r(K) ≤ c(r,d)? (docs/01_problem.md).

Collaborator: Pavlo Yatsyna (Charles University). Read docs in this order before work:
`docs/01_problem.md` → `docs/02_ledger.md` → `docs/03_corrections.md` → `docs/04_pitfalls.md`
→ `docs/05_roadmap.md`. Data: `docs/07_data.md`. History: `docs/08_history.md`.

## Working standards (non-negotiable)
1. **Proofs over enumeration.** Run a computation only if its outcome changes what we can
   prove (refutes a conjecture, certifies an extremal example, or tests a lemma we are about
   to prove). State beforehand what each possible outcome would establish.
2. **Tag every claim** with exactly one status: `PROVED` (proof written in repo), `IMPORTED`
   (from literature/other analysis, with reference; say whether we re-verified),
   `CERTIFIED` (exact computation, reproducible by a command in this repo), `EMPIRICAL`
   (pattern in data, no proof), `HEURISTIC` (argument with a gap — name the gap),
   `OPEN`, `REFUTED`.
3. **Before restating any result, check `docs/03_corrections.md`.** Several plausible
   statements have already been retracted (some twice). Do not resurrect them.
4. Be concise. No speculative numerology. Distinguish mod (K^×)² from mod (O_K^×)² counts
   at all times (only the former can be uniform).
5. Never call a search "complete" unless it provably covers all lattices (see pitfalls P5).
6. Update `docs/02_ledger.md` and `docs/03_corrections.md` whenever a status changes.

## Code
Pure Python 3, exact arithmetic (`fractions.Fraction`), no Sage/PARI required
(PARI/Sage may be added; keep the pure-Python path as a cross-check).
```
src/rqf.py      Q(√D) arithmetic (√D-basis), CF of quadratic irrationals, fundamental unit (CF), tp unit
src/indec.py    indecomposables (Dress–Scharlau semiconvergents), exact unit-orbit normalization,
                square-class reps (seeds α and αη)
src/gram.py     exact rank over K, admissible off-diagonals, abstract Gram search, realize, SOS
src/lattice.py  balanced reps, find_winner, realized module L=ΣO v_i (HNF), trace lattice,
                exact HNF / membership, Fincke–Pohst, analyze_lattice (classes, orbits,
                Z-linear irreducibility, n*)
scripts/experiments.py   CLI reproducing every number in docs/07_data.md
legacy/real_quadratic_gram.py   Sage version (corrected) — reference only
```
Tests: `python3 -m pytest -q -m "not slow"` (seconds); `python3 -m pytest -q -m slow` (~4 min).
Conventions: ω = √D (D ≡ 2,3 mod 4) or (1+√D)/2 (D ≡ 1 mod 4); Δ = 4D or D;
η = generator > 1 of totally positive units (η = ε if N(ε)=+1 else ε²).
Dress–Scharlau: indecomposables have N(α) ≤ Δ/4.
Blomer–Kala count (verified D<300 in tests): #indecomposables mod O^{×,+} =
u1+u3+…+u_{s−1} (s even) or u1+…+u_s (s odd), [u0; u1..us] the CF of ω.

## Known performance limits
`find_winner` is exponential in the number of reps (fine to ~20 reps); `largest_gram(58)`
is slow. Fincke–Pohst on the rank-4 trace lattice is fine to ~5·10^5 short vectors.
For larger D: use the frame/norm-equation model (docs/05_roadmap.md task R1), not subset search.
