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

Headline status: s□(M) ≤ 12·s proved (s = CF period); s□ computable exactly; certified
S(Q(√43),2) ≥ 8; bounded vs ≍ s open. Next: roadmap R1 (provably complete search).
