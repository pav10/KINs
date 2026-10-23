# 06 — References (what each is used for)

- **Dress–Scharlau**, indecomposables in real quadratic fields: indecomposables = semiconvergents
  of ω; N(α) ≤ Δ/4. (indec.py; P6.)
- **Blomer–Kala, arXiv:1705.03671**: number of indecomposables mod O^{×,+} via partial quotients;
  N(ε) = (−1)^s. Verified in tests for D < 300 (formula in CLAUDE.md).
- **Kala–Tinková, arXiv:2005.12312, Prop. 3.1**: every indecomposable of a real quadratic field has
  δ ∈ O^{∨,+} with Tr(δα) = 1; indecomposables on one edge share δ. Basis of L1.4–L1.5, L6.1;
  computationally checked by assertions in `edge_deltas`. Their cubic results are for simplest cubic
  fields (min traces 1 or 2), NOT "t₃ ≤ 3" (retracted, C32).
- **Kala–Man, arXiv:2403.18390**: sails S_K = ∂Conv(O⁺) in all signatures; indecomposables lie on
  the sail (converse iff d = 2); ι(K) ≫ Δ^{1/12} for some biquadratics but ≍ log Δ for a family.
- **O'Meara 1980, "On indecomposable quadratic forms"**: H(M) ⊆ A(M) ⊆ G(M); 3.7 (2M-cosets),
  3.17 (#A(M_Z) ≤ 2(2^n−1)), 3.18 (linear irreducibility and the trace lattice). Used in L1.1, L3.1,
  L3.2.
- **Shyr, Bull. AMS 83 (1977)**: class number H of proper classes in the genus of a free totally
  positive binary lattice: H = (h_K/h_k)·2^{1−r}/[U_K^0:U_k^0] (k base field, K = k(√−δ), r = #
  ramified primes); explicit Shintani-type formula. Used in L2.2, R5.
- **Yamamoto, arXiv:2304.02299**: angles realizable by integer vectors; dimension 3 exceptional
  (Hasse–Minkowski). Assessed: pairwise only, does not bound counts.
- **arXiv:2409.11082**: unit-representation regime (recovered by L1.3).
- **Duke–Schulze-Pillot; Cogdell et al. arXiv:1402.1332; Blomer–Harcos GAFA 2010; Harcos**:
  representation by ternary/Hilbert forms; effective ranges stop at Siegel — the analytic wall
  for L5.3.
- **Kala–Yatsyna arXiv:2402.03850**: sums of squares; same below-effective-bound class-number
  (exceptional element) phenomenon as L5.3.
- **Friesen**: prescribed continued fractions (constructions for R6 / L5.1 and for long-period
  families in R3(b)).
- Imported "registry" analyses (rounds 1 and E-family): Thms 1–6, 15–17, 21, 22, Prop 5, A7 —
  see ledger for which were re-verified.
- **R3 (docs/11_R3.md):**
  - **Herglotz, Math. Z. 12 (1922)**; **Lemmermeyer, "Kuroda's class number formula", Acta Arith. 66 (1994)**: h(E) = ½·Q·h(K)h(k1)h(k2) for imaginary biquadratic E with w_E = 2 (Theorem 11.1(a)). Re-verified against PARI for 4055 + 4808 fields.
  - **Genus theory** (Gauss; e.g. Cox, *Primes of the form x² + ny²*, §3, §6): h(Q(√p)) is odd for p prime; disc −8p has 2 genera, with principal genus = Cl² (Theorem 11.3). Class field theory: a norm map is onto the (narrow) class group when the extension is ramified at a finite prime.
  - **Zagier, "Nombres de classes et fractions continues", Astérisque 24–25 (1975)**: reduced forms and partial quotients. With h(K) = 1, Σ_{i≤s} a_i = T(D)/2 = ½Σ_{|u|<√D} d(D − u²).
  - **Hirzebruch, Enseign. Math. 19 (1973)** (proof: Zagier, loc. cit.): for p ≡ 3 mod 4 prime with h(Q(√p)) = 1, Σ(−1)^i a_i = 3h(−p).
  - Both identities are used only in §4–§5 and are re-verified for all h = 1 primes ≤ 10^6 in the scan.
- **R7 (docs/15_R7.md):**
  - **Tinková, arXiv:2105.04204 (Ramanujan J. 61, 2023)**: min traces of indecomposables are unbounded in families of
    cubic orders (abstract only; the paper could not be fetched in the R7 session). Context for Theorem 15.6.
  - **E. Thomas, J. reine angew. Math. 310 (1979)**: the family x(x − a)(x − b) − 1 and its units. Source of the family only.
  - **Chan–Kim–Raghavan 1996**: the real quadratic fields with a classical universal ternary form are Q(√2),
    Q(√3), Q(√5). IMPORTED, not re-verified; the list was confirmed from secondary sources (e.g. arXiv:2501.19371). Used in Prop. 15.9.
  - **Blomer–Kala arXiv:1503.05736; Kala–Tinková arXiv:2005.12312**: trace pairing with codifferent elements as the
    tool for lower bounds on ranks of universal forms (context for Cor. 15.8).
  - **PARI/GP `nflist`** (cubic fields by discriminant): IMPORTED as complete for the scan of §3.3.
  - Standard facts used without proof: Euler's lemma Z[ρ]^∨ = f′(ρ)⁻¹Z[ρ]; the norm-2 vectors of an integral lattice
    without norm-1 vectors form an ADE root system; disc of the trace form = N(δ^r 𝔳M)|d_K|^r.
