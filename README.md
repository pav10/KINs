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

Headline status: s□(M) ≤ 12·s proved (s = CF period); R1 done: S(K,2) computed exactly
(provably complete search, docs/09_R1_complete.md) for all D ≤ 100 — max 8 (D = 43, 67, 86) —
and long-period fields; "≤ 2 classes per edge" refuted (D = 82); bounded vs ≍ s still open.
Next: roadmap R2 (Galois-symmetric extremal lattices).

```bash
python3 scripts/complete_S.py out.jsonl 43 67          # exact S(K,2), S_free(K,2)
SAILS_PATH=../sails python3 scripts/check_sails.py 2 3000   # needs cypari2 + pav10/sails
```
