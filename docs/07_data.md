# 07 — Data (all reproducible; commands in parentheses)

## A. Ground truth for validation (tests/)
- κ□ (square classes of indecomposables, seeding α, αη): D=2:2, 3:2, 5:1, 6:4, 7:4, 10:5, 11:6,
  13:3, 14:4, 15:2, 17:3, 19:10, 21:2, 22:8, 23:4, 26:9, 29:5, 30:4, 31:12, 43:22, 58:17, 67:34.
- Max rank-2 abstract Gram on class reps (`gram.largest_gram`): D=19:4, 22:4, 31:4, 58:6
  (D = 43: 8 on balanced reps).
- Blomer–Kala ι formula: verified all squarefree D < 300.
- Field table D ≤ 100 (s, ι, κ□, N(ε), Reg): data/fields_D100.txt
  (`scripts/experiments.py fields --Dmax 100`).

## B. Certified extremal lattices (exact, `lattice.exact_sq_classes`)
| K | lattice | s□ | per-edge values | edges touched | lin-irred classes |
|---|---|---|---|---|---|
| Q(√43), s=10 | L in Gp=[[7+√43,−2],[−2,7−√43]] (balanced winner), N(𝔳L)=1 | **8** | 1,1,2,1,1,1,1,2,1,1 | 10/10 | 2 (cap 3000) |
| Q(√19), s=6 | orbit-rep winner, det Gp = 170+39√19 = η | 6 | 1,2,1,1,2,1 | 6/6 | 2 |
| Q(√31), s=8 | balanced winner, Gp = diag(6+√31, 17−3√31) ambient | 6 | 1,2,1,0,1,2,1,0 | 6/8 | 1 ([6+√31]) |
| Q(√57), s=6 | frame α1=10+3ω, α2=23+7ω, b=−(13+4ω), δ=151+20√57=η (imported, re-verified) | 6 | 2,2,2,2,2,2 | 6/6 | 0 |
Q(√43) classes: norm 6: 7−√43, 1541−235√43, 47207−7199√43, 10731517−1636541√43;
norm 17: 105−16√43, 282−43√43, 730938−111467√43, 1963743−299468√43.
Q(√31) value traces: 12, 34, 56, 490, 902, 1314. Q(√57) classes: 8±√57, 68±9√57, (23−3√57)/2,
(53−7√57)/2 (each class met on two edges).

## C. Exact scan, lower bounds for S(K,2), D ≤ 66 (`scripts/scan_exact.py 2 66 25 out.jsonl`)
Lattice = winner of `find_winner` on balanced class reps / on orbit reps (25 s limit each), s□
evaluated exactly. LOWER BOUNDS ONLY (fixed reps; e.g. D=57 gives 4 here but 6 by the frame in B).
Raw data with Gram matrices and per-edge counts: data/scan_exact.jsonl.

| D | s | ι | κ□ | N(ε) | s□ (balanced) | s□ (orbit) | edges touched |
|---|---|---|---|---|---|---|---|
| 2 | 1 | 2 | 2 | -1 | 2 | 2 | 1 |
| 3 | 2 | 1 | 2 | +1 | 2 | 1 | - |
| 5 | 1 | 1 | 1 | -1 | 1 | 1 | - |
| 6 | 2 | 2 | 4 | +1 | 2 | 2 | 2 |
| 7 | 4 | 2 | 4 | +1 | 2 | 2 | 3 |
| 10 | 1 | 6 | 5 | -1 | 2 | 2 | 1 |
| 11 | 2 | 3 | 6 | +1 | 2 | 2 | 2 |
| 13 | 1 | 3 | 3 | -1 | 1 | 1 | 1 |
| 14 | 4 | 2 | 4 | +1 | 2 | 2 | 3 |
| 15 | 2 | 1 | 2 | +1 | 2 | 1 | - |
| 17 | 3 | 5 | 3 | -1 | 2 | 2 | 3 |
| 19 | 6 | 7 | 10 | +1 | 4 | 6 | 6 |
| 21 | 2 | 1 | 2 | +1 | 1 | 1 | - |
| 22 | 6 | 6 | 8 | +1 | 4 | 2 | 4 |
| 23 | 4 | 2 | 4 | +1 | 2 | 2 | 3 |
| 26 | 1 | 10 | 9 | -1 | 2 | 2 | 1 |
| 29 | 1 | 5 | 5 | -1 | 1 | 1 | 1 |
| 30 | 2 | 2 | 4 | +1 | 2 | 2 | 2 |
| 31 | 8 | 8 | 12 | +1 | 6 | 2 | 5 |
| 33 | 4 | 4 | 4 | +1 | 2 | 2 | 4 |
| 34 | 4 | 2 | 4 | +1 | 2 | 2 | 3 |
| 35 | 2 | 1 | 2 | +1 | 2 | 1 | - |
| 37 | 3 | 7 | 5 | -1 | 2 | 2 | 3 |
| 38 | 2 | 6 | 12 | +1 | 4 | 2 | 2 |
| 39 | 2 | 4 | 8 | +1 | 2 | 2 | 2 |
| 41 | 5 | 11 | 7 | -1 | 4 | 2 | 5 |
| 42 | 2 | 2 | 4 | +1 | 2 | 2 | 2 |
| 43 | 10 | 13 | 22 | +1 | timeout | 8 | 10 |
| 46 | 12 | 8 | 12 | +1 | 4 | 2 | 5 |
| 47 | 4 | 2 | 4 | +1 | 2 | 2 | 3 |
| 51 | 2 | 7 | 12 | +1 | 4 | 2 | 2 |
| 53 | 1 | 7 | 7 | -1 | 1 | 1 | 1 |
| 55 | 4 | 4 | 6 | +1 | 2 | 2 | 2 |
| 57 | 6 | 7 | 10 | +1 | 4 | 4 | 3 |
| 58 | 7 | 20 | 17 | -1 | timeout | timeout |  |
| 59 | 6 | 9 | 14 | +1 | 4 | 2 | 4 |
| 61 | 3 | 11 | 9 | -1 | 4 | 1 | 3 |
| 62 | 4 | 2 | 4 | +1 | 2 | 2 | 3 |
| 65 | 3 | 9 | 6 | -1 | 2 | 2 | 3 |
| 66 | 2 | 8 | 14 | +1 | 4 | 2 | 2 |

## D. Earlier experiments (superseded but recorded)
- Orbit plateau (winner on orbit reps, capped FP; `experiments.py analyze --reps orbit`):
  orbits represented 2,3,4,4,4 at s = 2,4,6,8,12 (D = 66,33,19,31,46). Classes are the right count:
  e.g. D=19 orbit winner = 6 classes; D=31 orbit winner = 2 classes.
- Period-2 fans (`experiments.py fan 66 146`): best lattice 2 classes for ι = 8, 12.
- Odd convergents jointly rank-2 compatible (`experiments.py convergents 19 31 43 46`):
  2/3, 2/4, 3/5, 2/6.
- Imported (other analysis, Thm 22, unit-determinant frames, not re-run): complete S: D=26:2, 33:4,
  19:6; others: 57:6, 41:4, 61:4, 31:≥4; 22,34,38,39,47,55,62,65,74,79,83: 2–3;
  35,51,53,66,69,70,77,78,82,85,87: no unit-determinant frame.
