# indec-binary-pack

Representation of (square classes of) indecomposable integers by binary quadratic lattices over
real quadratic fields — research hand-over pack (Yatsyna + Claude).

Start: `CLAUDE.md`, then `docs/01_problem.md` … `docs/08_history.md`.

```bash
python3 -m pytest -q -m "not slow"      # seconds
python3 -m pytest -q -m slow            # ~4 min (D=31, Gram sizes)
python3 scripts/experiments.py frame57  # Δ=57 certificate
python3 scripts/scan_exact.py 43 43 600 /tmp/o.jsonl   # S(Q(√43),2) ≥ 8, exact
```
Requirements: Python ≥ 3.9, numpy not required, pytest for tests. No Sage/PARI needed
(`legacy/real_quadratic_gram.py` is the corrected Sage original, reference only).

Headline status: s□(M) ≤ 4·s proved and attained (s = CF period; R7, was 12s); R1 done: S(K,2) computed exactly
(provably complete search, docs/09_R1_complete.md) for all D ≤ 100 — max 8 (D = 43, 67, 86) —
and long-period fields. R2 (docs/10_R2.md) gives the mechanism behind large s□: for
D ≡ 3 mod 4 the represented indecomposables are norms of ideals in one class of Cl(K(√−2)), so the
best lattice comes from a class-group computation. R3 (docs/11_R3.md) makes this elementary for p prime:
max_c n(c) = max over form classes f of discriminant −8p of #{indecomposable square classes whose norm f
represents} — pure Python, no class group of E. Certified S(K,2) ≥ 64 (p = 814531, complete scan to 10^6)
and **S(K,2) ≥ 76** (p = 7956859), all binary forms of det ε. The mean of max_c n(c) grows like 2.8·ln p.
Growth along a family is NOT proved: it needs class-number-one-type input (long principal Zagier cycles,
ι/√p → ∞) and a twisted-divisor mean value (docs/11_R3.md §5). "≤ 2 classes per edge" refuted (D = 82).
R4 (docs/12_R4.md):
- the imported Thm 22 table was re-run; D = 22 corrected to 4;
- the Ramsey reduction and "no minimum angle" are proved;
- the imported regulator ceiling is replaced by the proved s□ ≤ 12s ≤ 24 + 24.94·Reg (R7: 4s ≤ 8 + 8.32·Reg).

R5 (docs/13_R5.md):
- the genus of ⟨1, ε⟩ has exactly h(−8D)/2 proper classes, the ideal lattices over ker(Cl E → Cl⁺K);
- Σ_gen s□ is explicit;
- the maximum exceeds the genus average by 1–4×.

R6 (docs/14_R6.md):
- sup_K T(K) = ∞ is proved, so the count mod unit squares is not uniform in any rank;
- T(K) ≤ 3s (R7: 2s);
- an exact witness test for indecomposability.

R7 (docs/15_R7.md):
- sharp per-face lemma for classical lattices: ≤ rd values per functional, equality iff (M, Tr δB) ≅ Z^{rd};
  hence S(K,2) ≤ 4s (3s if N(ε) = +1), attained at D = 82;
- degree d: per-functional counts are uniform in M (packing); the field enters through t(K) = max min trace;
- cubic fields: "t₃ ≤ 3" is false (t = 4 at disc 3137, t = 5 at disc 14229; complete scan to 20000), and
  t(K) ≥ a in Thomas's family x(x−a)(x−b)−1 when Z[ρ] = O_K (proved; certified a ≤ 100, so t₃ ≥ 100);
- Kitaoka: classical universal rank-r ⇒ every edge partial quotient ≤ 2r − 1; indecomposables alone cannot
  decide Kitaoka in degree 2 (D = n² − 1).

```bash
python3 scripts/complete_S.py out.jsonl 43 67          # exact S(K,2), S_free(K,2)
python3 scripts/r2_classes.py 5419                      # class distribution over Cl(K(sqrt-2)) (needs cypari2)
python3 scripts/certify_ideal.py 98179                  # build + certify the best ideal lattice
python3 scripts/r3_tools.py scan out.jsonl 3 100000 0/1  # nu* = max_c n(c) for primes p = 3 mod 4 (pure Python)
python3 scripts/certify_form.py 7956859                 # PARI-free certificate S >= 76 (~20 s)
python3 scripts/r3_analysis.py                          # tables of docs/11_R3.md
SAILS_PATH=../sails python3 scripts/check_sails.py 2 3000   # needs cypari2 + pav10/sails
```
